#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uwoc_lut_gen.py — Generates the ROM for uwoc_channel.v from uwoc_channel_model.py
═══════════════════════════════════════════════════════════════════════════════

Why generate it offline?
  σ²_ho is a double integral over the Nikishov spectrum — it has NO closed form and
  cannot be evaluated on the FPGA. Nor can the inverse CDF of Gamma / Weibull / Lognormal.
  → compute them in Python here and load them into the FPGA as lookup tables.

Produces: verilog/uwoc_channel_rom.vh  with 6 tables

  ┌ hs_rom(cls, addr)      256×12b × 8 classes  inverse CDF of h_s (Gamma)
  ├ ho_rom(row, addr)      256×12b × 8 rows     inverse CDF of h_o (LN/Weibull)
  │                                             rows 0-5: turbulence levels L0-L5
  │                                             rows 6-7: matched-variance family
  │                                                       controls [rev B]
  ├ psig_rom(λ, w, dist)   24b × 3 × 3 × 16     nominal P_sig (μ_ref, h=1)
  ├ epol_rom(λ, w, dist)   24b × 3 × 3 × 16     e_pol(d) — polarisation error from scattering
  ├ hscls_rom(water, dist) 3b  × 3 × 16         (water,dist) → σ_s² class
  ├ hooff_rom(water, dist) 3b  × 3 × 16         range-based turbulence-level offset
  └ nexp_inv_rom(λ, w, d)  16b × 3 × 3 × 16     32768/E[clicks per window]
                                                (so channel_monitor can normalise SNR
                                                 with NO divider)

Fixed conventions (the RTL must follow them exactly):
  • h_s, h_o : 12-bit, value 256 = mean (E[h] = 1.0), saturating at 4095 (h < 16)
  • probability : 24-bit unsigned, value p ↔ round(p · 2^24), saturating at 2^24-1
  • ROM addressed by inverse CDF: addr u ∈ [0,255] ↔ quantile (u + 0.5)/256

[rev B, 2026-09-15]
  • ho_rom has 8 rows. Row 6 is the Weibull law at the L3 scintillation index and
    row 7 the log-normal law at the L4 index, so the pairs (3, 6) and (7, 4) share
    σ²_ho and differ only in the distribution family. uwoc_channel.v selects them
    with turb_level = 6 / 7, bypassing the range offset.
  • _rtl_exact_sim follows the RTL operation order ((psig·m)<<g >>3)·h_f >>8 and the
    RTL level clamp (1..5 for a nonzero command).

Run:
    python python/uwoc_lut_gen.py
    python python/uwoc_lut_gen.py --verify      # check the ROM against float
"""

from __future__ import annotations

import argparse
import datetime
import os
import sys

import numpy as np
from scipy.stats import gamma as gamma_dist, weibull_min, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uwoc_channel_model import (            # noqa: E402
    WATER_TYPES, OCEAN_TURB_LEVELS, LAMBDA_OPTIONS, LinkConfig, UWOCChannel,
    sigma2_s, sigma2_ho_clamped, lognormal_params, weibull_params,
    h_ell, e_pol, sample_h_s, sample_h_o, n_bar, p_click,
)

# ═══════════════════════════════════════════════════════════════════════════
# QUANTISATION CONVENTIONS  (the RTL uwoc_channel.v must match EXACTLY)
# ═══════════════════════════════════════════════════════════════════════════
#
# BIT WIDTHS — chosen by numerical experiment, NOT arbitrarily:
#
#  • h_s, h_o: 12-bit, H_MEAN = 256 ↔ h = 1.0, representing h ∈ [0, 16).
#    Why not 8-bit/128 like the old gamma_gamma_final.v? Because the
#    atmospheric model has light tails, while Gamma (large σ_s²) and Weibull
#    (large σ²_ho) underwater have VERY HEAVY TAILS: h_max reaches 16–20.
#    Clipping at h = 2.0 (8-bit/128) drops E[h] from 1.00 to 0.56 and σ² from
#    4.50 to 1.74 — badly wrong. At 12-bit/256 the E[h] error is ≤ 2.4% at σ² = 5.
#
#  • Probability: 24-bit, p ↔ round(p·2²⁴), resolution 6·10⁻⁸.
#    Why not 16-bit? Because underwater P_sig drops to ~10⁻⁷ at long range;
#    at 16 bits round(2.9e-6 · 65536) = 0 — the signal VANISHES entirely and
#    the emulator is left with dark counts only. This is a fundamental
#    difference from atmospheric FSO, where P_click ~ O(1).
#
H_BITS = 12
H_MEAN = 256            # 12-bit: 256 ↔ h = 1.0
H_MAX = (1 << H_BITS) - 1
ROM_DEPTH = 256         # inverse-CDF points
PROB_BITS = 24
PROB_SCALE = 1 << PROB_BITS
PROB_MAX = PROB_SCALE - 1
N_DIST = 16             # range points per water type
N_WATER = 3
N_HS_CLASS = 8          # quantisation classes for σ_s²
MU_REF_STEP = 8         # mu_level = 8 ↔ the nominal μ of LinkConfig

WATER_ORDER = ("clear_ocean", "coastal", "harbor")

# σ_s² class grid. Log-spaced because σ_s² varies exponentially with range.
HS_CLASS_SIGMA2 = (0.0, 0.02, 0.05, 0.12, 0.30, 0.75, 1.80, 4.50)


def q_prob(p: float) -> int:
    """Quantise a probability to 24 bits."""
    return int(np.clip(round(float(p) * PROB_SCALE), 0, PROB_MAX))


def q_h(x: float) -> int:
    """Quantise a fading coefficient to 12 bits (256 = 1.0)."""
    return int(np.clip(round(float(x) * H_MEAN), 0, H_MAX))


# ═══════════════════════════════════════════════════════════════════════════
# 1. INVERSE-CDF ROM
# ═══════════════════════════════════════════════════════════════════════════

def inv_cdf_gamma(sigma2: float) -> list[int]:
    """h_s ~ Gamma(shape=1/σ², scale=σ²), E=1  → a 256-entry 12-bit table."""
    u = (np.arange(ROM_DEPTH) + 0.5) / ROM_DEPTH
    if sigma2 <= 1e-6:
        return [H_MEAN] * ROM_DEPTH
    v = gamma_dist.ppf(u, a=1.0 / sigma2, scale=sigma2)
    return [q_h(x) for x in v]


def inv_cdf_ho(sigma2: float, family: str | None = None) -> list[int]:
    """h_o inverse CDF, E=1 → a 256-entry 12-bit table.

    family=None applies the modelling rule of the paper (log-normal below
    σ² = 1, Weibull at and above). "lognormal" / "weibull" force the family,
    which is how the matched-variance ablation rows are built."""
    u = (np.arange(ROM_DEPTH) + 0.5) / ROM_DEPTH
    if sigma2 <= 1e-6:
        return [H_MEAN] * ROM_DEPTH
    if family is None:
        family = "lognormal" if sigma2 < 1.0 else "weibull"
    if family == "lognormal":
        mu_X, sX2 = lognormal_params(sigma2)
        v = np.exp(2.0 * norm.ppf(u, loc=mu_X, scale=np.sqrt(sX2)))
    elif family == "weibull":
        b1, b2 = weibull_params(sigma2)
        v = weibull_min.ppf(u, c=b1, scale=b2)
    else:
        raise ValueError("unknown family %r" % family)
    return [q_h(x) for x in v]


def ho_row_specs(cfg: LinkConfig) -> list[tuple[str, float, str]]:
    """(name, σ²_ho at 20 m, family) for the 8 ho_rom rows.

    Rows 0-5 follow the modelling rule. Rows 6 and 7 are the matched-variance
    family controls: Weibull at the L3 index and log-normal at the L4 index."""
    s2 = [sigma2_ho_clamped(20.0, t["eps"], t["chi_T"], t["w"],
                            cfg.lam_nm, cfg.wave_type)
          for t in (OCEAN_TURB_LEVELS[lv] for lv in range(6))]
    rows = []
    for lv in range(6):
        fam = ("none" if s2[lv] <= 1e-6
               else ("lognormal" if s2[lv] < 1.0 else "weibull"))
        rows.append((OCEAN_TURB_LEVELS[lv]["name"], s2[lv], fam))
    rows.append(("ABL_WB@L3", s2[3], "weibull"))
    rows.append(("ABL_LN@L4", s2[4], "lognormal"))
    return rows


def hs_class_of(sigma2: float) -> int:
    """Assign σ_s² to the nearest class in HS_CLASS_SIGMA2 (by log distance)."""
    ref = np.array(HS_CLASS_SIGMA2)
    d = np.abs(np.log(np.maximum(ref, 1e-4)) - np.log(max(sigma2, 1e-4)))
    if sigma2 < 1e-3:
        return 0
    return int(np.argmin(d))


# ═══════════════════════════════════════════════════════════════════════════
# 2. BUILD ALL THE TABLES
# ═══════════════════════════════════════════════════════════════════════════

def build_tables(cfg: LinkConfig, window: int) -> dict:
    tabs: dict = {}

    # --- inverse-CDF ---
    tabs["hs_rom"] = [inv_cdf_gamma(s2) for s2 in HS_CLASS_SIGMA2]
    specs = ho_row_specs(cfg)
    tabs["ho_specs"] = specs
    tabs["ho_rom"] = [inv_cdf_ho(s2, None if fam == "none" else fam)
                      for _, s2, fam in specs]
    # NOTE: ho_rom takes σ²_ho at the 20 m REFERENCE range. The range dependence
    # of rows 1-5 is handled separately by hooff_rom below (an integer level
    # offset) — this keeps the ROM small while still tracking σ²_ho(d).

    # --- WAVELENGTH-dependent tables: address = {lam_idx[1:0], water[1:0], dist[3:0]}
    # This λ dimension is what makes the wavelength-adaptation knob in
    # adaptive_controller.v actually affect the emulator (otherwise λ would be a
    # decorative variable). Cost: 2 × 3 λ × 3 waters × 16 ranges × 24 bit ≈ 6.9 kbit — cheap.
    psig, epol, nexp_inv = [], [], []
    for lam in LAMBDA_OPTIONS:
        c_lam = LinkConfig(**{**cfg.__dict__, "lam_nm": lam})
        for wi, water in enumerate(WATER_ORDER):
            r_psig, r_epol, r_ninv = [], [], []
            for d in WATER_TYPES[water].d_grid:
                d = float(d)
                p_sig = float(1.0 - np.exp(-n_bar(1.0, d, water, c_lam)))
                r_psig.append(q_prob(p_sig))
                r_epol.append(q_prob(e_pol(d, water, c_lam)))
                exp_clicks = (p_sig + c_lam.Y0) * window
                r_ninv.append(int(np.clip(round(32768.0 / max(exp_clicks, 1e-9)),
                                          1, 65535)))
            psig.append(r_psig); epol.append(r_epol); nexp_inv.append(r_ninv)

    # --- λ-independent tables: address = {water[1:0], dist[3:0]} ---
    # σ_s² is fitted from Table IV of [2] at the reference wavelength; that paper
    # does not characterise the λ dependence of scattering fading, so we keep it λ-independent.
    hscls, ho_lvl, meta_rows = [], [], []
    for wi, water in enumerate(WATER_ORDER):
        r_cls, r_ho = [], []
        for di, d in enumerate(WATER_TYPES[water].d_grid):
            d = float(d)
            s2s = sigma2_s(d, water)
            r_cls.append(hs_class_of(s2s))
            r_ho.append(_ho_level_offset(d, cfg))
            lam0 = LAMBDA_OPTIONS.index(cfg.lam_nm) if cfg.lam_nm in LAMBDA_OPTIONS else 0
            p_sig = psig[lam0 * N_WATER + wi][di] / PROB_SCALE
            meta_rows.append((water, di, d, s2s, r_cls[-1], p_sig,
                              e_pol(d, water, cfg), (p_sig + cfg.Y0) * window))
        hscls.append(r_cls); ho_lvl.append(r_ho)

    tabs.update(psig_rom=psig, epol_rom=epol, hscls_rom=hscls,
                nexp_inv_rom=nexp_inv, hooff_rom=ho_lvl, meta=meta_rows)
    return tabs


def _ho_level_offset(d: float, cfg: LinkConfig) -> int:
    """
    Range-based turbulence-level offset.

    ho_rom is generated at the 20 m reference range. σ²_ho grows almost linearly
    with d, and levels L1..L5 are ~4× apart in σ². So the level offset needed is

        offset = clip( round( log₄( σ²_ho(d) / σ²_ho(20 m) ) ), -3, 3 )

    evaluated with the L3 turbulence parameters at the configured wavelength. It
    depends on the range only, not on the water type. The RTL uses
        eff_level = clamp(turb_level + offset, 1, 5)   for turb_level = 1..5,
    and ignores it for the ablation rows 6/7. That way only 6 turbulence ROMs are
    needed while still reflecting σ²_ho(d) across the whole range grid.
    """
    t = OCEAN_TURB_LEVELS[3]
    ref = sigma2_ho_clamped(20.0, t["eps"], t["chi_T"], t["w"], cfg.lam_nm)
    cur = sigma2_ho_clamped(d, t["eps"], t["chi_T"], t["w"], cfg.lam_nm)
    if ref <= 0 or cur <= 0:
        return 0
    return int(np.clip(round(np.log(cur / ref) / np.log(4.0)), -3, 3))


# ═══════════════════════════════════════════════════════════════════════════
# 3. VERILOG EMISSION
# ═══════════════════════════════════════════════════════════════════════════

def _emit_rom_mem(name: str, out_bits: int, data: list[list[int]],
                  per_line: int = 8) -> str:
    """
    Emit one ROM as a MEMORY ARRAY + initial block.

    Why not a `function` + nested `case` like the old gamma_gamma_final.v?
    Quartus synthesises a nested case into LOGIC (LEs), not into M4K. With 2048
    12-bit entries the LE count explodes. The form `reg [W-1:0] mem [0:N-1]` +
    initial is reliably inferred by Quartus II as a ROM in M4K blocks.

    Flat address: addr = {sel, idx}  (sel in the high bits).
    """
    rows, cols = len(data), len(data[0])
    L = [f"reg [{out_bits-1}:0] {name} [0:{rows*cols-1}];",
         "initial begin"]
    flat = [v for row in data for v in row]
    for i in range(0, len(flat), per_line):
        body = " ".join(f"{name}[{i+j}]={out_bits}'d{v};"
                        for j, v in enumerate(flat[i:i + per_line]))
        L.append("    " + body)
    L.append("end")
    return "\n".join(L)


def emit_verilog(tabs: dict, cfg: LinkConfig, window: int, path: str):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out = [
        "// " + "=" * 74,
        "// AUTO-GENERATED — DO NOT EDIT BY HAND",
        "// Generated by: python/uwoc_lut_gen.py",
        f"// Timestamp: {now}",
        "// " + "=" * 74,
        "// Model: h = h_l(d,λ;w) · h_f,  h_f = h_s · h_o   (Salcedo-Serrano ICC'22 eq.2)",
        f"// λ = {cfg.lam_nm} nm   μ_ref = {cfg.mu}   η_det = {cfg.eta_det}",
        f"// Y₀ = {cfg.Y0:.6e}   f_rep = {cfg.f_rep:.3e} Hz",
        f"// Monitoring window used to compute nexp_inv_rom: {window:,} attempts",
        "//",
        "// Conventions:",
        f"//   h_s, h_o : {H_BITS}-bit, {H_MEAN} = mean (h = 1.0), saturating at {H_MAX}",
        f"//   probability: {PROB_BITS}-bit, p ↔ round(p·{PROB_SCALE}), saturating at {PROB_MAX}",
        f"//   mu_level : 4-bit, {MU_REF_STEP} = nominal μ",
        "//",
        "// water_type order: 0=clear ocean, 1=coastal, 2=harbor",
        "// " + "=" * 74,
        "",
        "// ---- Range grid [m] per water type (reference / UART only) ----",
    ]
    for wi, water in enumerate(WATER_ORDER):
        g = ", ".join(f"{d:g}" for d in WATER_TYPES[water].d_grid)
        out.append(f"//   water {wi} ({WATER_TYPES[water].name}): {g}")
    out.append("")

    out.append("// ---- σ_s² for each h_s class ----")
    for i, s2 in enumerate(HS_CLASS_SIGMA2):
        out.append(f"//   class {i}: σ_s² = {s2:g}")
    out.append("")

    out.append("// ---- σ²_ho and family of each ho_rom row (at d_ref = 20 m) ----")
    for row, (name, s2, fam) in enumerate(tabs["ho_specs"]):
        kind = {"none": "—", "lognormal": "Lognormal", "weibull": "Weibull"}[fam]
        tag = "   [ablation row, turb_level %d, no range offset]" % row if row >= 6 else ""
        out.append(f"//   row {row} {name:<9}: σ²_ho = {s2:7.4f}  ({kind}){tag}")
    out.append("")

    out.append(f"// hs_rom  : address = {{hs_class[2:0], addr[7:0]}}  (8×256)")
    out.append(_emit_rom_mem("hs_rom", H_BITS, tabs["hs_rom"]))
    out.append("")
    out.append(f"// ho_rom  : address = {{eff_level[2:0], addr[7:0]}}  "
               f"({len(tabs['ho_rom'])}×256)")
    out.append(_emit_rom_mem("ho_rom", H_BITS, tabs["ho_rom"]))
    out.append("")
    out.append("// ---- λ-DEPENDENT tables: address = {lam_sel[1:0], water[1:0], dist[3:0]}")
    out.append("//      lam_sel: 0=" + " 1=".join(str(l) for l in LAMBDA_OPTIONS) + " nm")
    out.append(_emit_rom_mem("psig_rom", PROB_BITS, tabs["psig_rom"], 4))
    out.append("")
    out.append(_emit_rom_mem("epol_rom", PROB_BITS, tabs["epol_rom"], 4))
    out.append("")
    out.append(f"// nexp_inv_rom = round(32768 / E[clicks per {window:,}-attempt window])")
    out.append("// → lets channel_monitor normalise SNR with NO divider.")
    out.append("// ⚠ These values are tied to the window size: the ATTEMPT_LOG2 parameter")
    out.append(f"//   of channel_monitor.v MUST be {int(np.log2(window))}.")
    out.append(_emit_rom_mem("nexp_inv_rom", 16, tabs["nexp_inv_rom"], 4))
    out.append("")
    out.append("// ---- λ-INDEPENDENT tables: address = {water[1:0], dist[3:0]} ----")
    out.append(_emit_rom_mem("hscls_rom", 3, tabs["hscls_rom"], 8))
    out.append("")
    # hooff_rom returns a 3-bit signed value (−3..3) → encoded as two's complement
    hooff = [[v & 0x7 for v in row] for row in tabs["hooff_rom"]]
    out.append("// hooff_rom: range-based turbulence offset, 3-bit two's complement (−3..+3)")
    out.append(_emit_rom_mem("hooff_rom", 3, hooff, 8))
    out.append("")

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")

    n_bits = (len(HS_CLASS_SIGMA2) * ROM_DEPTH * H_BITS
              + len(tabs["ho_rom"]) * ROM_DEPTH * H_BITS
              + len(LAMBDA_OPTIONS) * N_WATER * N_DIST
              * (PROB_BITS + PROB_BITS + 16)
              + N_WATER * N_DIST * (3 + 3))
    return n_bits


# ═══════════════════════════════════════════════════════════════════════════
# 4. VERIFICATION: does the ROM reproduce the float distribution?
# ═══════════════════════════════════════════════════════════════════════════

def verify(tabs: dict, cfg: LinkConfig, window: int) -> bool:
    print("\n" + "═" * 78)
    print("  ROM VERIFICATION — 12-bit tables vs floating-point model")
    print("═" * 78)
    rng = np.random.default_rng(2024)
    ok = True

    print(f"\n  (a) inverse-CDF h_s (Gamma) — {ROM_DEPTH} levels, {H_BITS}-bit"
          f" (maximum representable h = {H_MAX/H_MEAN:.2f})")
    print(f"      {'class':>6}{'σ_s²':>9}{'theoretical h_max':>18}{'E[h] ROM':>11}"
          f"{'σ² ROM':>10}{'E[h] error':>13}{'σ² error':>12}")
    print("      " + "-" * 80)
    for i, s2 in enumerate(HS_CLASS_SIGMA2):
        v = np.array(tabs["hs_rom"][i], dtype=float) / H_MEAN
        m, var = v.mean(), v.var() / max(v.mean() ** 2, 1e-12)
        hmax = (gamma_dist.ppf((ROM_DEPTH - 0.5) / ROM_DEPTH,
                               a=1.0 / s2, scale=s2) if s2 > 0 else 1.0)
        e_m = abs(m - 1.0) * 100
        e_v = abs(var - s2) / max(s2, 1e-9) * 100 if s2 > 0 else 0.0
        ok &= (e_m < 5.0 and (e_v < 10.0 or s2 == 0))
        print(f"      {i:>6}{s2:>9.3f}{hmax:>18.2f}{m:>11.4f}{var:>10.4f}"
              f"{e_m:>12.2f}%{e_v:>11.2f}%")

    print(f"\n  (b) inverse-CDF h_o — {ROM_DEPTH} levels, {H_BITS}-bit"
          " (rows 6-7: matched-variance distribution-family controls)")
    print(f"      {'row':>6}{'σ²_ho':>10}{'distribution':>13}{'E[h] ROM':>11}"
          f"{'σ² ROM':>10}{'E[h] error':>13}{'σ² error':>12}")
    print("      " + "-" * 74)
    for row, (name, s2, fam) in enumerate(tabs["ho_specs"]):
        v = np.array(tabs["ho_rom"][row], dtype=float) / H_MEAN
        m, var = v.mean(), v.var() / max(v.mean() ** 2, 1e-12)
        e_m = abs(m - 1.0) * 100
        e_v = abs(var - s2) / max(s2, 1e-9) * 100 if s2 > 1e-6 else 0.0
        kind = {"none": "—", "lognormal": "LogNormal", "weibull": "Weibull"}[fam]
        ok &= (e_m < 5.0 and (e_v < 12.0 or s2 <= 1e-6))
        print(f"      {row:>6}{s2:>10.4f}{kind:>13}{m:>11.4f}{var:>10.4f}"
              f"{e_m:>12.2f}%{e_v:>11.2f}%")

    print(f"\n  (c) RTL-exact golden model (integers + ROM only) vs floating-point model")
    print("      Consider only VIABLE operating points (P_click ≥ 1e-4); at longer")
    print("      ranges, the link is dead due to attenuation; deviations are only shot noise.")
    print(f"      {'water':<13}{'d':>6}{'lv':>4}{'P_click float':>15}"
          f"{'P_click ROM':>14}{'#click':>9}{'QBER float':>12}"
          f"{'QBER ROM':>11}{'Δ':>8}")
    print("      " + "-" * 84)
    N_SIM = 2_000_000
    n_checked = 0
    for wi, water in enumerate(WATER_ORDER):
        wt = WATER_TYPES[water]
        for di in range(0, N_DIST, 2):
            d = float(wt.d_grid[di])
            for lv in (1, 4):
                fl = UWOCChannel(water, d, lv, cfg).mean_metrics()
                if fl["p_click"] < 1e-4:
                    continue
                rt = _rtl_exact_sim(tabs, wi, di, lv, cfg, n=N_SIM, rng=rng,
                                    lam_idx=LAMBDA_OPTIONS.index(cfg.lam_nm))
                dq = abs(rt["qber"] - fl["qber"]) * 100
                dp = abs(rt["p_click"] - fl["p_click"]) / fl["p_click"] * 100
                ok &= (dq < 3.0 and dp < 25.0)
                n_checked += 1
                print(f"      {water:<13}{d:>6.1f}{lv:>4}"
                      f"{fl['p_click']:>15.4e}{rt['p_click']:>14.4e}"
                      f"{rt['n_click']:>9d}{fl['qber']*100:>11.2f}%"
                      f"{rt['qber']*100:>10.2f}%{dq:>7.2f}%")
    print(f"      ({n_checked} operating points checked)")

    print(f"\n  {'PASS' if ok else 'REVIEW NEEDED'} — ROM is sufficiently accurate for RTL")
    print("  Note: P_click ROM < float under strong turbulence is CORRECT — because")
    print("  fading is heavy-tailed; most h samples are BELOW the mean (median < mean).")
    return ok


def _rtl_exact_sim(tabs: dict, wi: int, di: int, turb_level: int,
                   cfg: LinkConfig, n: int, rng, mu_level: int = MU_REF_STEP,
                   lam_idx: int = 0, mu_gain_log2: int = 0):
    """
    Simulate EXACTLY AS THE RTL DOES: integers only, using the generated ROMs.
    This is the golden model for cross-checking uwoc_channel.v.
    """
    cls = tabs["hscls_rom"][wi][di]
    off = tabs["hooff_rom"][wi][di]
    if turb_level == 0:
        eff_lv = 0
    elif turb_level >= 6:                       # ablation rows: verbatim
        eff_lv = turb_level
    else:                                       # RTL clamp is 1..5
        eff_lv = int(np.clip(turb_level + off, 1, 5))

    hs = np.array(tabs["hs_rom"][cls], dtype=np.int64)[
        rng.integers(0, ROM_DEPTH, n)]
    ho = np.array(tabs["ho_rom"][eff_lv], dtype=np.int64)[
        rng.integers(0, ROM_DEPTH, n)]

    # h_f = (h_s · h_o) >> 8 , saturated to 12 bits          (RTL §4 step 1)
    hf = np.minimum((hs * ho) >> 8, H_MAX)

    # p_sig = ((((psig_ref · mu_level) << g) >> 3) · h_f) >> 8  (RTL §4 step 2)
    lw = lam_idx * N_WATER + wi                 # address {lam, water}
    psig_ref = int(tabs["psig_rom"][lw][di])
    psig_mu = min(((psig_ref * mu_level) << mu_gain_log2) >> 3, PROB_MAX)
    psig = np.minimum((psig_mu * hf) >> 8, PROB_MAX)

    p_noise = q_prob(cfg.Y0)
    epol = tabs["epol_rom"][lw][di]

    sig_det = rng.integers(0, PROB_SCALE, n) < psig
    noise_det = rng.integers(0, PROB_SCALE, n) < p_noise
    click = sig_det | noise_det

    r = rng.integers(0, PROB_SCALE, n)
    err = np.where(sig_det, r < epol, r < (PROB_SCALE // 2))

    nc = int(click.sum())
    return dict(p_click=nc / n, n_click=nc,
                qber=float((err & click).sum() / nc) if nc else 0.5)


# ═══════════════════════════════════════════════════════════════════════════
# 5. MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description="Generate UWOC ROM for FPGA")
    ap.add_argument("--out", default="verilog/uwoc_channel_rom.vh")
    ap.add_argument("--lam", type=int, default=450, choices=LAMBDA_OPTIONS)
    ap.add_argument("--mu", type=float, default=0.1, help="nominal μ")
    ap.add_argument("--window", type=int, default=65536,
                    help="monitoring window [attempts]; MUST be a power of 2 and "
                         "match 2^ATTEMPT_LOG2 in channel_monitor.v")
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()

    cfg = LinkConfig(lam_nm=args.lam, mu=args.mu)
    window = args.window or cfg.coherence_pulses

    print("█" * 78)
    print("  UWOC LUT GENERATOR")
    print(f"  λ={cfg.lam_nm} nm  μ_ref={cfg.mu}  η={cfg.eta_det}  "
          f"Y₀={cfg.Y0:.3e}")
    print(f"  monitoring window = {window:,} attempts")
    print("█" * 78)

    tabs = build_tables(cfg, window)

    print("\n  Lookup table by (water type, range):")
    print(f"  {'water':<13}{'idx':>4}{'d [m]':>8}{'σ_s²':>10}{'class':>5}"
          f"{'P_sig':>12}{'e_pol':>9}{'E[click]/window':>13}")
    print("  " + "-" * 76)
    for (water, di, d, s2s, cls, psig, ep, exp_c) in tabs["meta"]:
        if di % 3 == 0 or di == N_DIST - 1:
            print(f"  {water:<13}{di:>4}{d:>8.2f}{s2s:>10.4f}{cls:>5}"
                  f"{psig:>12.3e}{ep*100:>8.2f}%{exp_c:>13.1f}")

    n_bits = emit_verilog(tabs, cfg, window, args.out)
    print(f"\n  ✔ Wrote {args.out}")
    print(f"    Total ROM capacity ≈ {n_bits:,} bits "
          f"= {n_bits/4096:.1f} M4K blocks (Cyclone II EP2C20 has 52 blocks)")

    good = True
    if args.verify:
        good = verify(tabs, cfg, window)

    print("\n  Next step: compile verilog/uwoc_channel.v in Quartus")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
