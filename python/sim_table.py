#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sim_table.py — Full-matrix SIMULATION table + COMPARISON table against FPGA measurements
═══════════════════════════════════════════════════════════════════════════════

DOES NOT access the COM port: this script only reads data/fpga_points.csv
(+ data/clicks_*.csv) and recomputes the model. It is safe to run in parallel
with an fpga_collect.py measurement session.

Three tables, each with a different purpose:

  [1] --matrix  SIMULATION TABLE (data/sim_table.csv)
      The complete water × turbulence level × distance × wavelength matrix,
      calculated using the analytical model + window statistics. This is the
      "theoretical results" table for the paper, and it also indicates WHICH
      POINTS ARE WORTH MEASURING in hardware: rows that differ only in their
      std/outage columns but are identical in their P_click/QBER columns will
      also produce exactly the same curve when measured again on the FPGA.

  [2] --compare  COMPARISON TABLE (data/compare_table.csv + .md)
      Matches every measured point with its model prediction: P_click ratio,
      Poisson z-score, whether the measured QBER's Clopper-Pearson interval
      contains the model QBER, and the fading field's σ² inferred from the irrad
      column in clicks_*.csv.

  [3] --blocks  COHERENCE-BLOCK TABLE (data/block_table.csv)
      Only works with data acquired using the [v12] bitstream and --coh >= 1:
      each block of 2^(coh+5) qubits is then ONE frozen-fading sample, so each
      block's QBER is a genuinely measured quantity. The table compares the
      blockwise mean / standard deviation / outage against
      UWOCChannel.window_statistics(window = 2^(coh+5)). THIS is where L1…L5
      separate; POOLED QBER does not — see the note below.

Run:
    python python/sim_table.py                      # all three tables
    python python/sim_table.py --matrix             # simulation table only
    python python/sim_table.py --compare            # comparison table only
    python python/sim_table.py --blocks             # blockwise table only
    python python/sim_table.py --matrix --lams 450,532,650

WHY POOLED QBER DOES NOT REVEAL TURBULENCE:
    fpga_collect writes total errors / total sifted bits, i.e. an average
    WEIGHTED by the number of clicks. Deeply faded blocks contribute very few
    clicks and therefore receive little weight, so the final number barely
    changes with turbulence level (the model at 25 m gives 3.71% at L1 and
    3.70% at L5). The UNWEIGHTED mean across individual blocks instead rises
    from 3.72% to 13.50% — and that is the security-determining quantity,
    because BB84 distills keys block by block rather than pooling an entire
    session.

SECONDARY INDICATOR — CV(irrad):
    h = h_s·h_o with E[h_s] = E[h_o] = 1 and the two variables independent, so
        Var(h) = (1 + σ²_s)(1 + σ²_ho) − 1
    The irrad field is reported PER CLICK (top_module.v, scale 128 = 1.0), so
    its histogram correctly samples the marginal distribution of h. This
    quantity also works with data acquired using an older bitstream, for which
    blockwise statistics were not yet measurable.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np                                          # noqa: E402
from scipy.stats import beta as _beta, norm as _norm        # noqa: E402

from uwoc_channel_model import (                            # noqa: E402
    OCEAN_TURB_LEVELS, WATER_TYPES, LinkConfig, UWOCChannel,
    h2, max_secure_distance, sigma2_ho_clamped, skr,
)
from uwoc_lut_gen import _ho_level_offset                   # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data")

WATER_ORDER = ("clear_ocean", "coastal", "harbor")
QBER_LIMIT = 0.11          # BB84 / Shor-Preskill
OUTAGE_LIMIT = 0.10
CONF = 0.95
MIN_SIFT = 16

SIM_FIELDS = [
    "water", "water_name", "distance_m", "turb", "turb_name", "lam_nm",
    "sigma2_s", "sigma2_ho", "sigma2_h", "cv_h", "L_dB", "nbar", "e_pol_pct",
    "p_click", "sift_eff", "qber_pct", "skr_bps",
    "qber_win_std_pct", "qber_win_p90_pct", "outage", "skr_outage_bps",
    "secure_mean", "secure_outage",
]

CMP_FIELDS = [
    "tag", "phase", "water", "distance_m", "turb", "lam_nm",
    "n_qubit", "n_click", "n_sift", "n_err",
    "p_click_meas", "p_click_sim", "p_click_ratio", "z_click",
    "qber_meas_pct", "qber_lo_pct", "qber_hi_pct", "qber_sim_pct",
    "qber_in_CI", "skr_meas_pp", "skr_sim_pp",
    "cv_h_meas", "cv_h_sim", "sigma2_h_meas", "sigma2_h_sim", "irrad_sat",
    "n_blocks", "var_click", "seconds", "secure_meas",
]


# ═══════════════════════════════════════════════════════════════════════════
# Statistics helpers
# ═══════════════════════════════════════════════════════════════════════════
def cp_interval(k: int, n: int, conf: float = CONF) -> tuple[float, float]:
    """Two-sided Clopper-Pearson interval — exact, so it stays honest at n_err = 0."""
    if n <= 0:
        return 0.0, 1.0
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else float(_beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k >= n else float(_beta.ppf(1 - a, k + 1, n - k))
    return lo, hi


def sigma2_h(s2_s: float, s2_ho: float) -> float:
    """Var(h_s·h_o) for two INDEPENDENT unit-mean factors."""
    return (1.0 + s2_s) * (1.0 + s2_ho) - 1.0


def rom_channel(water: str, d: float, turb: int, cfg: LinkConfig) -> UWOCChannel:
    """The channel the FPGA ACTUALLY implements, not the continuous model.

    ho_rom holds one inverse-CDF per turbulence level, generated at the 20 m
    reference range only; range dependence is carried by an INTEGER level offset
    round(log₄(σ²_ho(d)/σ²_ho(20 m))) — see _ho_level_offset in uwoc_lut_gen.py and
    uwoc_channel.v §3. So σ²_ho on hardware is quantised onto the ladder
    {0.02, 0.08, 0.30, 1.00, 3.00}, in steps of 4×.

    At 25 m the offset is 0, so the board draws L4 from σ²_ho = 1.00 while the
    continuous model at that range says 1.55 — a 55 % difference in variance, which
    lands squarely on the deep-fade tail that sets the outage. Comparing measured
    block statistics against the continuous model there charges the RTL for a
    discrepancy that is really the ROM's range quantisation, so this is the
    reference the block table uses.
    """
    ch = UWOCChannel(water, d, turb, cfg)
    eff = max(1, min(5, turb + _ho_level_offset(d, cfg)))
    t = OCEAN_TURB_LEVELS[eff]
    ch.sigma2_ho = sigma2_ho_clamped(20.0, t["eps"], t["chi_T"], t["w"],
                                     cfg.lam_nm, cfg.wave_type)
    return ch


def irrad_stats(tag: str) -> tuple[int, float, float]:
    """(n, cv, saturated fraction) of the per-click irrad field.

    irrad is h on a 128 = 1.0 scale, logged once per click. Since clicks are spread
    over many coherence blocks, its histogram samples the marginal law of h, so
    CV(irrad)² estimates Var(h).

    ⚠ THE FIELD CLIPS AT 255 (top_module.v: h_f[8:1] saturated), i.e. at h = 2.0.
    At weak turbulence nothing reaches that, but σ²_ho = 3 at L5 puts a large part
    of the distribution above it and the clipped CV then UNDER-estimates the true
    spread — badly. The saturated fraction is returned so the caller can say so
    instead of quietly reporting a wrong number; the per-block QBER statistics
    (--blocks) are immune, which is another reason to prefer them.
    """
    p = os.path.join(DATA_DIR, "clicks_%s.csv" % tag)
    if not os.path.exists(p):
        return 0, float("nan"), float("nan")
    with open(p, newline="", encoding="utf-8") as f:
        v = np.array([int(r["irrad"]) for r in csv.DictReader(f)
                      if r.get("irrad", "").strip().isdigit()], dtype=float)
    if v.size < 8:
        return int(v.size), float("nan"), float("nan")
    m = float(v.mean())
    cv = float(v.std() / m) if m > 0 else float("nan")
    return int(v.size), cv, float((v >= 255).mean())


def _short(path: str) -> str:
    """Path relative to the repo when possible — relpath raises across drives."""
    try:
        return os.path.relpath(path, ROOT)
    except ValueError:
        return path


def _write_csv(path: str, fields: list[str], rows: list[dict]):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("  → %s  (%d rows)" % (_short(path), len(rows)))


# ═══════════════════════════════════════════════════════════════════════════
# [1] Simulation matrix
# ═══════════════════════════════════════════════════════════════════════════
def build_sim_row(water: str, d: float, turb: int, lam_nm: int,
                  cfg0: LinkConfig, n_windows: int, mc: int) -> dict:
    cfg = LinkConfig(**{**cfg0.__dict__, "lam_nm": lam_nm})
    ch = UWOCChannel(water, d, turb, cfg)
    m = ch.mean_metrics()
    w = ch.window_statistics(n_windows=n_windows, qber_limit=QBER_LIMIT)

    s2h = sigma2_h(m["sigma2_s"], m["sigma2_ho"])
    q_basis = 0.5                                   # p_z = 0.5, both scripts use this
    row = dict(
        water=water, water_name=WATER_TYPES[water].name, distance_m=d,
        turb=turb, turb_name=OCEAN_TURB_LEVELS[turb]["name"], lam_nm=lam_nm,
        sigma2_s=m["sigma2_s"], sigma2_ho=m["sigma2_ho"], sigma2_h=s2h,
        cv_h=math.sqrt(s2h), L_dB=m["L_dB"], nbar=m["nbar"],
        e_pol_pct=100 * m["e_det"], p_click=m["p_click"],
        sift_eff=q_basis * m["p_click"], qber_pct=100 * m["qber"],
        skr_bps=skr(m["p_click"], m["qber"], cfg, q=q_basis),
        qber_win_std_pct=100 * w["qber_std"],
        qber_win_p90_pct=100 * w["qber_p90"],
        outage=w["outage"],
        skr_outage_bps=w["skr_per_pulse"] * cfg.f_rep,
        secure_mean=m["qber"] < QBER_LIMIT,
        secure_outage=(m["qber"] < QBER_LIMIT) and (w["outage"] < OUTAGE_LIMIT),
    )
    if mc:
        r = ch.simulate(mc, seed=1234 + int(d * 7) + turb)
        row["mc_p_click"] = r["p_click"]
        row["mc_qber_pct"] = 100 * r["qber"]
        row["mc_n"] = mc
    return row


def cmd_matrix(args, cfg0: LinkConfig):
    waters = [w.strip() for w in args.waters.split(",") if w.strip()]
    turbs = [int(t) for t in args.turbs.split(",") if t.strip()]
    lams = [int(l) for l in args.lams.split(",") if l.strip()]

    rows = []
    for water in waters:
        grid = WATER_TYPES[water].d_grid
        for lam in lams:
            for turb in turbs:
                for d in grid:
                    rows.append(build_sim_row(water, float(d), turb, lam,
                                              cfg0, args.windows, args.mc))

    fields = list(SIM_FIELDS) + (["mc_p_click", "mc_qber_pct", "mc_n"]
                                 if args.mc else [])
    print("\n" + "═" * 108)
    print("  SIMULATION TABLE — %d configurations" % len(rows))
    print("═" * 108)

    # Print one compact block per (water, λ); turbulence levels side by side so the
    # invariance of the mean is visible at a glance instead of buried in the CSV.
    for water in waters:
        for lam in lams:
            sub = [r for r in rows if r["water"] == water and r["lam_nm"] == lam]
            if not sub:
                continue
            print("\n  ▸ %s, λ = %d nm, μ = %g"
                  % (WATER_TYPES[water].name, lam, cfg0.mu))
            print("  %6s %4s %10s %8s %8s %9s %9s %8s %8s %6s"
                  % ("d[m]", "L", "P_click", "QBER%", "SKR bps", "σ²_ho",
                     "CV(h)", "σQBER%", "outage", "secure"))
            print("  " + "-" * 96)
            for r in sub:
                print("  %6g L%-3d %10.3e %8.2f %8.3g %9.3f %9.3f %8.2f %8.3f %6s"
                      % (r["distance_m"], r["turb"], r["p_click"], r["qber_pct"],
                         r["skr_bps"], r["sigma2_ho"], r["cv_h"],
                         r["qber_win_std_pct"], r["outage"],
                         "✓" if r["secure_outage"] else
                         ("mean" if r["secure_mean"] else "✗")))

    # Headline table for the paper: d_max per (water, turbulence level).
    print("\n" + "═" * 108)
    print("  MAXIMUM SECURE DISTANCE d_max [m] — outage criterion < %.0f%% "
          "(mean criterion in parentheses)" % (100 * OUTAGE_LIMIT))
    print("═" * 108)
    for lam in lams:
        cfg = LinkConfig(**{**cfg0.__dict__, "lam_nm": lam})
        print("\n  λ = %d nm" % lam)
        print("  %-18s" % "Water" + "".join("%16s" % ("L%d" % l) for l in turbs))
        print("  " + "-" * (18 + 16 * len(turbs)))
        for water in waters:
            cells = []
            for lv in turbs:
                d_out = max_secure_distance(water, lv, cfg, QBER_LIMIT,
                                            criterion="outage",
                                            outage_limit=OUTAGE_LIMIT)
                d_mean = max_secure_distance(water, lv, cfg, QBER_LIMIT,
                                             criterion="mean")
                cells.append("%7.1f (%5.1f)" % (d_out, d_mean))
            print("  %-18s" % WATER_TYPES[water].name
                  + "".join("%16s" % c for c in cells))
    print("\n  The parenthesized column barely changes with L: mean QBER is")
    print("  independent of turbulence (P_click is nearly linear in h, so fading")
    print("  cancels in expectation). Only the outage criterion separates the turbulence levels.\n")

    _write_csv(os.path.join(DATA_DIR, args.sim_out), fields, rows)
    return rows


# ═══════════════════════════════════════════════════════════════════════════
# [2] Measured vs simulated
# ═══════════════════════════════════════════════════════════════════════════
def cmd_compare(args, cfg0: LinkConfig):
    path = os.path.join(DATA_DIR, args.points)
    if not os.path.exists(path):
        print("  ⚠ %s does not exist yet — skipping the comparison table." % _short(path))
        return []

    with open(path, newline="", encoding="utf-8") as f:
        meas = list(csv.DictReader(f))
    if not meas:
        print("  ⚠ %s is empty." % _short(path))
        return []

    rows = []
    for r in meas:
        water = r["water"]
        d = float(r["distance_m"])
        turb = int(r["turb"])
        lam = int(r["lam_nm"])
        n_q = int(r["n_qubit"])
        n_c = int(r["n_click"])
        n_s = int(r["n_sift"])
        n_e = int(r["n_err"])

        cfg = LinkConfig(**{**cfg0.__dict__, "lam_nm": lam})
        ch = UWOCChannel(water, d, turb, cfg)
        m = ch.mean_metrics()

        # Var(n_click) has TWO terms once the fading block is long:
        #     n_click = Σ_b Bin(N/B, p·h_b)   over B coherence blocks
        #     Var = N·p                    (Poisson / counting)
        #         + (N·p)²·σ²_h / B        (only B independent fading draws)
        # At the legacy coh_sel = 0 the second term vanishes (B ~ 10^5). At
        # coh_sel = 11 a 5 m point spans only ~25 blocks and the fading term
        # dominates: scoring it as pure Poisson would flag a perfectly good point
        # as a 3σ discrepancy.
        s2h_sim = sigma2_h(m["sigma2_s"], m["sigma2_ho"])
        coh_q = int(float(r.get("coh_qubits") or 0))
        n_blocks = (n_q / coh_q) if coh_q > 0 else 0.0
        exp_click = n_q * m["p_click"]
        var_click = exp_click + (exp_click ** 2 * s2h_sim / n_blocks
                                 if n_blocks > 0 else 0.0)
        z = (n_c - exp_click) / math.sqrt(max(var_click, 1e-9))
        lo, hi = cp_interval(n_e, n_s)
        q_meas = n_e / n_s if n_s else float("nan")
        q_basis = 0.5

        n_ir, cv_meas, sat = irrad_stats(r["tag"])
        s2h_sim = sigma2_h(m["sigma2_s"], m["sigma2_ho"])

        rows.append(dict(
            tag=r["tag"], phase=r["phase"], water=water, distance_m=d,
            turb=turb, lam_nm=lam, n_qubit=n_q, n_click=n_c, n_sift=n_s,
            n_err=n_e,
            p_click_meas=n_c / max(n_q, 1), p_click_sim=m["p_click"],
            p_click_ratio=(n_c / max(n_q, 1)) / max(m["p_click"], 1e-300),
            z_click=z,
            qber_meas_pct=100 * q_meas, qber_lo_pct=100 * lo,
            qber_hi_pct=100 * hi, qber_sim_pct=100 * m["qber"],
            qber_in_CI=bool(n_s >= MIN_SIFT and lo <= m["qber"] <= hi),
            skr_meas_pp=float(r["skr_per_pulse"]),
            skr_sim_pp=q_basis * m["p_click"] * max(0.0, 1 - 2 * h2(m["qber"])),
            cv_h_meas=cv_meas, cv_h_sim=math.sqrt(s2h_sim),
            sigma2_h_meas=cv_meas ** 2 if cv_meas == cv_meas else float("nan"),
            sigma2_h_sim=s2h_sim, n_blocks=n_blocks, var_click=var_click,
            irrad_sat=sat,
            seconds=float(r["seconds"]), secure_meas=r["secure"],
        ))

    # ---- console -------------------------------------------------------
    print("\n" + "═" * 118)
    print("  MEASUREMENT (FPGA) ↔ SIMULATION COMPARISON — %d points" % len(rows))
    print("═" * 118)
    print("  %-22s %5s %10s %10s %6s %7s %8s %-18s %8s %5s"
          % ("configuration", "d[m]", "P_meas", "P_model", "ratio", "z",
             "QBER meas", "CP 95%", "QBER mdl", "match"))
    print("  " + "-" * 114)
    for r in rows:
        tag = r["tag"].replace("fixed_", "").replace("adaptive_", "")
        verdict = ("match" if r["qber_in_CI"] else
                   ("MISMATCH" if r["n_sift"] >= MIN_SIFT else "low N"))
        print("  %-22s %5g %10.3e %10.3e %6.3f %+7.2f %7.2f%% [%5.2f%%,%6.2f%%] "
              "%7.2f%% %5s"
              % (tag, r["distance_m"], r["p_click_meas"], r["p_click_sim"],
                 r["p_click_ratio"], r["z_click"], r["qber_meas_pct"],
                 r["qber_lo_pct"], r["qber_hi_pct"], r["qber_sim_pct"], verdict))

    tot_c = sum(r["n_click"] for r in rows)
    tot_ce = sum(r["n_qubit"] * r["p_click_sim"] for r in rows)
    tot_s = sum(r["n_sift"] for r in rows if r["n_sift"] >= MIN_SIFT)
    tot_e = sum(r["n_err"] for r in rows if r["n_sift"] >= MIN_SIFT)
    tot_ee = sum(r["n_sift"] * r["qber_sim_pct"] / 100
                 for r in rows if r["n_sift"] >= MIN_SIFT)
    n_ok = sum(1 for r in rows if r["qber_in_CI"])
    n_val = sum(1 for r in rows if r["n_sift"] >= MIN_SIFT)
    # Variances add across points; each already carries its own fading term.
    tot_var = sum(r["var_click"] for r in rows)
    z_tot = (tot_c - tot_ce) / math.sqrt(max(tot_var, 1e-9))
    print("  " + "-" * 114)
    print("  POOLED: clicks %d / expected %.1f = %.3f (z = %+.2f, p = %.3f)"
          % (tot_c, tot_ce, tot_c / max(tot_ce, 1e-9), z_tot,
             2 * float(_norm.sf(abs(z_tot)))))
    if any(r["n_blocks"] and r["n_blocks"] < 100 for r in rows):
        print("  ⚠ Some points span < 100 coherence blocks: the z-score above ALREADY")
        print("    includes fading dispersion; do not compare it with the old pure-Poisson z-score.")
    if tot_s:
        lo, hi = cp_interval(tot_e, tot_s)
        print("  POOLED: QBER %.3f%% [%.3f%%, %.3f%%] vs model %.3f%% → %s"
              % (100 * tot_e / tot_s, 100 * lo, 100 * hi, 100 * tot_ee / tot_s,
                 "match" if lo <= tot_ee / tot_s <= hi else "MISMATCH"))
        print("  %d/%d sufficiently sampled points contain the model within the 95%% CP interval "
              "(expected ~95%%)" % (n_ok, n_val))

    # ---- turbulence signature -----------------------------------------
    have_cv = [r for r in rows if r["cv_h_meas"] == r["cv_h_meas"]]
    if have_cv:
        print("\n" + "═" * 118)
        print("  TURBULENCE INDICATOR — σ² of h inferred from the per-click irrad column")
        print("═" * 118)
        print("  %-22s %5s %4s %10s %10s %10s %10s %8s"
              % ("configuration", "d[m]", "L", "CV meas", "CV model",
                 "σ²_h meas", "σ²_h mdl", "clipped"))
        print("  " + "-" * 94)
        n_sat = 0
        for r in have_cv:
            sat = r.get("irrad_sat", float("nan"))
            flag = " ←" if sat == sat and sat > 0.01 else ""
            if flag:
                n_sat += 1
            print("  %-22s %5g L%-3d %10.3f %10.3f %10.3f %10.3f %7.1f%%%s"
                  % (r["tag"].replace("fixed_", "").replace("adaptive_", ""),
                     r["distance_m"], r["turb"], r["cv_h_meas"], r["cv_h_sim"],
                     r["sigma2_h_meas"], r["sigma2_h_sim"],
                     100 * sat if sat == sat else float("nan"), flag))
        if n_sat:
            print("\n  ← %d points have >1%% of clicks CLIPPED at irrad = 255 (h = 2.0)."
                  % n_sat)
            print("    The measured CV in those rows is a LOWER BOUND, not a correct estimate.")
            print("    Use the --blocks table for these points: blockwise QBER is not clipped.")
        print("\n  For data acquired using the OLD bitstream, this is the only quantity")
        print("  that separates L1…L5 at the same distance. From [v12] onward, the --blocks table is far better.")

    _write_csv(os.path.join(DATA_DIR, args.cmp_out), CMP_FIELDS, rows)
    _write_markdown(os.path.join(DATA_DIR, os.path.splitext(args.cmp_out)[0]
                                 + ".md"), rows)
    return rows


# ═══════════════════════════════════════════════════════════════════════════
# [3] Per-coherence-block statistics — the turbulence result [v12]
# ═══════════════════════════════════════════════════════════════════════════
WIN_FIELDS = [
    "tag", "water", "distance_m", "turb", "lam_nm", "coh_qubits", "n_qubit",
    "seg_mode", "n_blocks", "n_blocks_expected", "n_blocks_valid",
    "sift_per_block",
    "qber_pooled_pct", "qber_win_mean_pct", "qber_win_std_pct",
    "qber_win_p90_pct", "outage", "n_dead",
    "sim_qber_win_mean_pct", "sim_qber_win_std_pct", "sim_qber_win_p90_pct",
    "sim_outage", "sim_qber_pooled_pct", "sim_sigma2_ho_rom",
    "cont_qber_win_mean_pct", "cont_qber_win_std_pct", "cont_outage",
]


def blocks_from_clicks(tag: str, coh_q: int, n_qubit: int, n_click: int):
    """Per-coherence-block (n_sift, n_err) reconstructed from the per-click log.

    Two segmentation modes, picked automatically:

    'attempt' — the correct one. Block = attempt // coh_qubits, with the TOTAL block
        count taken from n_qubit rather than from the log, so that blocks producing
        NO click still appear. Those are the deep fades that drive the outage
        probability; dropping them would bias every dispersion statistic downwards.

    'irrad' — the salvage path. Bitstreams before the qub_index fix reported
        top_module's total_cnt in that field, and total_cnt only advances in
        FSM_PROCESS, a state a lost photon never reaches — so the column counts
        CLICKS and every click lands in block 0. h is constant inside a block, so a
        run of consecutive clicks sharing one irrad value IS one block, and the data
        is still usable. Two limits, both stated by the caller: a block that
        produced zero clicks is invisible, and two adjacent blocks that happened to
        draw the same irrad byte merge into one.

    Returns (sift, errs, mode, n_blk_expected) or None.
    """
    p = os.path.join(DATA_DIR, "clicks_%s.csv" % tag)
    if not os.path.exists(p) or coh_q <= 0 or n_qubit <= 0:
        return None
    with open(p, newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r.get("bmatch", "").strip()]
    if not rows:
        return None

    n_blk = int(math.ceil(n_qubit / float(coh_q)))
    att = [r.get("attempt", "").strip() for r in rows]
    have_att = all(a for a in att)

    # A click counter maxes out at n_click; a real attempt counter reaches n_qubit.
    # Underwater the two differ by ~3 decades, so the test is not delicate.
    if have_att and max(int(a) for a in att) <= 1.05 * max(n_click, 1) \
            and n_click < 0.5 * n_qubit:
        have_att = False

    if have_att:
        sift = np.zeros(n_blk, dtype=np.int64)
        errs = np.zeros(n_blk, dtype=np.int64)
        for r, a in zip(rows, att):
            b = int(a) // coh_q
            if b >= n_blk:                  # line that arrived after the last chunk
                continue
            if int(r["bmatch"]):
                sift[b] += 1
                errs[b] += int(r["err"])
        return sift, errs, "attempt", n_blk

    # ---- salvage: one block per run of constant irrad ----
    sift_l, errs_l = [], []
    prev = None
    for r in rows:
        ir = int(r["irrad"])
        if ir != prev:
            sift_l.append(0)
            errs_l.append(0)
            prev = ir
        if int(r["bmatch"]):
            sift_l[-1] += 1
            errs_l[-1] += int(r["err"])
    return (np.array(sift_l, dtype=np.int64), np.array(errs_l, dtype=np.int64),
            "irrad", n_blk)


def cmd_windows(args, cfg0: LinkConfig):
    path = os.path.join(DATA_DIR, args.points)
    if not os.path.exists(path):
        print("  ⚠ %s does not exist yet — skipping the blockwise table." % _short(path))
        return []
    with open(path, newline="", encoding="utf-8") as f:
        meas = list(csv.DictReader(f))

    rows, skipped = [], 0
    for r in meas:
        coh_q = int(float(r.get("coh_qubits") or 0))
        n_qubit = int(r["n_qubit"])
        blk = blocks_from_clicks(r["tag"], coh_q, n_qubit, int(r["n_click"]))
        if blk is None:
            skipped += 1
            continue
        sift, errs, mode, n_blk_exp = blk

        # Match window_statistics() exactly: a block with no sifted click carries no
        # information and is scored at QBER = 0.5, so it counts as an outage rather
        # than quietly vanishing from the average.
        valid = sift > 0
        qw = np.where(valid, errs / np.maximum(sift, 1), 0.5)
        pooled = errs.sum() / max(sift.sum(), 1)

        water, d, turb, lam = (r["water"], float(r["distance_m"]),
                               int(r["turb"]), int(r["lam_nm"]))
        cfg = LinkConfig(**{**cfg0.__dict__, "lam_nm": lam})
        ch = rom_channel(water, d, turb, cfg)
        w = ch.window_statistics(n_windows=40000, window=coh_q,
                                 qber_limit=QBER_LIMIT)
        sim_pooled = float((w["qber_win"] * w["clicks_win"]).sum()
                           / max(w["clicks_win"].sum(), 1))
        # Continuous model kept alongside: the gap between the two IS the ROM's
        # range quantisation, and it belongs in the paper's error budget.
        wc = UWOCChannel(water, d, turb, cfg).window_statistics(
            n_windows=40000, window=coh_q, qber_limit=QBER_LIMIT)

        # How many blocks are genuinely dead rather than merely noisy. The mean and
        # the std are both driven by this handful of blocks, so its Poisson error
        # is the real uncertainty of the whole row.
        n_dead = int((sift == 0).sum() + ((sift > 0) & (errs / np.maximum(sift, 1)
                                                        > 0.35)).sum())

        rows.append(dict(
            tag=r["tag"], water=water, distance_m=d, turb=turb, lam_nm=lam,
            coh_qubits=coh_q, n_qubit=n_qubit, n_blocks=int(sift.size),
            n_blocks_expected=n_blk_exp, seg_mode=mode,
            n_blocks_valid=int(valid.sum()),
            sift_per_block=float(sift.mean()),
            qber_pooled_pct=100 * pooled,
            qber_win_mean_pct=100 * float(qw.mean()),
            qber_win_std_pct=100 * float(qw.std()),
            qber_win_p90_pct=100 * float(np.percentile(qw, 90)),
            outage=float((qw > QBER_LIMIT).mean()),
            n_dead=n_dead,
            sim_qber_win_mean_pct=100 * w["qber_mean"],
            sim_qber_win_std_pct=100 * w["qber_std"],
            sim_qber_win_p90_pct=100 * w["qber_p90"],
            sim_outage=w["outage"],
            sim_qber_pooled_pct=100 * sim_pooled,
            sim_sigma2_ho_rom=ch.sigma2_ho,
            cont_qber_win_mean_pct=100 * wc["qber_mean"],
            cont_qber_win_std_pct=100 * wc["qber_std"],
            cont_outage=wc["outage"],
        ))

    if not rows:
        print("\n  ⚠ No point has both the `attempt` and `coh_qubits` columns "
              "(%d points skipped)." % skipped)
        print("    The [v12] bitstream and fpga_collect.py running with --coh >= 1 are required.")
        return []

    print("\n" + "═" * 118)
    print("  COHERENCE-BLOCK STATISTICS — where turbulence is actually revealed")
    print("═" * 118)
    print("  %-20s %5s %4s %13s %7s %9s %9s %9s %9s"
          % ("configuration", "d[m]", "L", "blocks", "sift/blk", "pooled QBER",
             "block QBER", "block std", "outage"))
    print("  " + "-" * 112)
    for r in rows:
        nb = ("%d" % r["n_blocks"] if r["seg_mode"] == "attempt"
              else "%d/%d*" % (r["n_blocks"], r["n_blocks_expected"]))
        print("  %-20s %5g L%-3d %13s %7.1f %8.2f%% %8.2f%% %8.2f%% %9.3f"
              % (r["tag"].replace("fixed_", "").replace("adaptive_", ""),
                 r["distance_m"], r["turb"], nb, r["sift_per_block"],
                 r["qber_pooled_pct"], r["qber_win_mean_pct"],
                 r["qber_win_std_pct"], r["outage"]))
        print("  %-20s %5s %4s %13s %7s %8.2f%% %8.2f%% %8.2f%% %9.3f   ← ROM (σ²_ho=%.2f)"
              % ("", "", "", "", "", r["sim_qber_pooled_pct"],
                 r["sim_qber_win_mean_pct"], r["sim_qber_win_std_pct"],
                 r["sim_outage"], r["sim_sigma2_ho_rom"]))
        print("  %-20s %5s %4s %13s %7s %8s %8.2f%% %8.2f%% %9.3f   ← continuous model"
              % ("", "", "", "", "", "", r["cont_qber_win_mean_pct"],
                 r["cont_qber_win_std_pct"], r["cont_outage"]))
        # The mean and the std both ride on this handful of blocks; its Poisson
        # error is the honest uncertainty of the row.
        nd = r["n_dead"]
        print("  %-20s dead blocks (QBER>35%%): %d ± %.1f  →  contribute %.2f%% ± %.2f%% "
              "to block QBER"
              % ("", nd, max(nd, 1) ** 0.5,
                 100 * nd * 0.5 / max(r["n_blocks"], 1),
                 100 * max(nd, 1) ** 0.5 * 0.5 / max(r["n_blocks"], 1)))
        if r["sift_per_block"] < 32:
            print("  %-20s ⚠ only %.1f sifted bits/block — blockwise QBER is dominated "
                  "by counting noise" % ("", r["sift_per_block"]))

    n_salv = sum(1 for r in rows if r["seg_mode"] == "irrad")
    if n_salv:
        print("\n  * %d points had to be SALVAGED by segmenting runs of constant irrad:"
              % n_salv)
        print("    the bitstream used for acquisition did not yet have qub_index, so the "
              "`attempt` field counts CLICKS, not attempts. Each consecutive irrad run")
        print("    is treated as one block. This has two consequences:")
        print("      · blocks with NO clicks are invisible ⇒ std and outage are LOWER BOUNDS")
        print("      · two adjacent blocks with the same irrad value are merged into one")
        print("    The 'blocks' column reports <measured>/<theoretical = n_qubit/coh_qubits>.")
    if skipped:
        print("\n  (%d points skipped because their logs do not yet have the `attempt` column.)" % skipped)
    print("\n  'Pooled QBER' is the value written by fpga_collect: total errors / total sifted bits.")
    print("  'Block QBER' is the UNWEIGHTED mean over individual coherence blocks —")
    print("  the security-determining quantity, because BB84 distills keys block by block.")

    _write_csv(os.path.join(DATA_DIR, args.win_out), WIN_FIELDS, rows)
    return rows


def _write_markdown(path: str, rows: list[dict]):
    """Paper-ready table: measurement and model in one row, CI attached."""
    lines = [
        "| Configuration | d [m] | L | λ [nm] | N qubits | Measured P_click | Model P_click "
        "| Ratio | Measured QBER (95% CP) | Model QBER | Match |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--:|",
    ]
    for r in rows:
        lines.append(
            "| %s | %g | %d | %d | %s | %.3e | %.3e | %.3f | %.2f %% [%.2f, %.2f] "
            "| %.2f %% | %s |"
            % (r["tag"].replace("fixed_", "").replace("adaptive_", ""),
               r["distance_m"], r["turb"], r["lam_nm"],
               "{:,}".format(r["n_qubit"]).replace(",", " "),
               r["p_click_meas"], r["p_click_sim"], r["p_click_ratio"],
               r["qber_meas_pct"], r["qber_lo_pct"], r["qber_hi_pct"],
               r["qber_sim_pct"], "✅" if r["qber_in_CI"] else
               ("—" if r["n_sift"] < MIN_SIFT else "❌")))
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("  → %s" % _short(path))


# ═══════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser(
        description="UWOC-BB84 simulation tables and comparison against FPGA measurements "
                    "(DOES NOT open the COM port)")
    ap.add_argument("--matrix", action="store_true", help="output only the simulation table")
    ap.add_argument("--compare", action="store_true", help="output only the comparison table")
    ap.add_argument("--blocks", dest="do_blocks", action="store_true",
                    help="output only the COHERENCE-BLOCK statistics table (requires data "
                         "acquired using the [v12] bitstream with --coh >= 1)")
    ap.add_argument("--waters", default=",".join(WATER_ORDER))
    ap.add_argument("--turbs", default="1,2,3,4,5")
    ap.add_argument("--lams", default="450")
    ap.add_argument("--mu", type=float, default=0.1)
    ap.add_argument("--windows", type=int, default=4000,
                    help="number of Monte Carlo windows for σQBER/outage statistics")
    ap.add_argument("--mc", type=int, default=0,
                    help="add qubit-level Monte Carlo columns, N pulses per point "
                         "(0 = disabled; only for checking the analytical calculation)")
    ap.add_argument("--points", default="fpga_points.csv",
                    help="measurement file inside the --data directory")
    ap.add_argument("--data", default=None,
                    help="data directory (default: <repo>/data)")
    ap.add_argument("--sim-out", dest="sim_out", default="sim_table.csv")
    ap.add_argument("--cmp-out", dest="cmp_out", default="compare_table.csv")
    ap.add_argument("--win-out", dest="win_out", default="block_table.csv")
    args = ap.parse_args()

    global DATA_DIR
    if args.data:
        DATA_DIR = os.path.abspath(args.data)

    cfg0 = LinkConfig(mu=args.mu)
    both = not (args.matrix or args.compare or args.do_blocks)

    print("█" * 108)
    print("  SIMULATION / COMPARISON TABLES — μ = %g, η_det = %g, Y₀ = %.2e, "
          "f_rep = %.0e Hz" % (cfg0.mu, cfg0.eta_det, cfg0.Y0, cfg0.f_rep))
    print("  QBER limit = %.0f%%, outage limit = %.0f%%, CP %.0f%%"
          % (100 * QBER_LIMIT, 100 * OUTAGE_LIMIT, 100 * CONF))
    print("█" * 108)

    if args.matrix or both:
        cmd_matrix(args, cfg0)
    if args.compare or both:
        cmd_compare(args, cfg0)
    if args.do_blocks or both:
        cmd_windows(args, cfg0)


if __name__ == "__main__":
    main()
