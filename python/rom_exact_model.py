#!/usr/bin/env python3
"""
rom_exact_model.py -- population model of the emulator at the RTL level.

Everything here uses the tables that are actually in verilog/uwoc_channel_rom.vh
and the integer arithmetic of uwoc_channel.v (rev B):

    h_f   = min((h_s * h_o) >> 8, 4095)
    P_sig = min((((psig_ref * m) << g) >> 3) * h_f >> 8, 2^24 - 1)
    decoy: P_sig >> shift          vacuum: 0
    click = sig | noise,  err = sig ? (u < e_pol) : coin

so a prediction made here is the expectation of what the board measures, not of
the continuous channel model. Two uses:

  * the equal-resolution comparison between model and measurement (the "L4
    discrepancy" of the first draft was an 8-bin measurement compared with a
    full-resolution model), and
  * planning: matched-variance ablation, background dose-response and the
    decoy-state campaign.

Nothing in this module talks to the board.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
ROM_PATH = ROOT / "verilog" / "uwoc_channel_rom.vh"

S24 = 1 << 24
P_MAX = S24 - 1
H_MAX = 4095
Y0_CODE = 218                 # round(1.3e-5 * 2^24)
MU_NOM = 0.1                  # mu_level = 8, gain 1
ATTEMPTS_PER_BLOCK = 1 << 16
LN2 = math.log(2.0)


# ===========================================================================
# ROM access
# ===========================================================================

@dataclass
class Rom:
    hs: np.ndarray        # [8, 256]
    ho: np.ndarray        # [rows, 256]
    psig: np.ndarray      # [9, 16]  index lam*3 + water
    epol: np.ndarray      # [9, 16]
    hscls: np.ndarray     # [3, 16]
    hooff: np.ndarray     # [3, 16]  signed

    def hs_class(self, water: int, dist: int) -> int:
        return int(self.hscls[water, dist])

    def offset(self, water: int, dist: int) -> int:
        return int(self.hooff[water, dist])

    def psig_ref(self, lam: int, water: int, dist: int) -> int:
        return int(self.psig[lam * 3 + water, dist])

    def epol_code(self, lam: int, water: int, dist: int) -> int:
        return int(self.epol[lam * 3 + water, dist])


def _table(text: str, name: str) -> np.ndarray:
    pairs = re.findall(r"(?<![A-Za-z_])" + name + r"\[(\d+)\]=\d+'d(\d+);", text)
    if not pairs:
        raise ValueError("table %s not found" % name)
    arr = np.zeros(max(int(i) for i, _ in pairs) + 1, dtype=np.int64)
    for i, v in pairs:
        arr[int(i)] = int(v)
    return arr


def load_rom(path: Path | str = ROM_PATH) -> Rom:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    ho = _table(text, "ho_rom")
    off = _table(text, "hooff_rom").reshape(3, 16)
    off = np.where(off >= 4, off - 8, off)             # 3-bit two's complement
    return Rom(hs=_table(text, "hs_rom").reshape(8, 256),
               ho=ho.reshape(ho.size // 256, 256),
               psig=_table(text, "psig_rom").reshape(9, 16),
               epol=_table(text, "epol_rom").reshape(9, 16),
               hscls=_table(text, "hscls_rom").reshape(3, 16),
               hooff=off)


def effective_row(turb: int, offset: int) -> int:
    if turb == 0:
        return 0
    if turb >= 6:
        return turb
    return int(min(5, max(1, turb + offset)))


def fading_values(rom: Rom, water: int, dist: int, turb: int) -> np.ndarray:
    """All equally likely 12-bit h_f values of one operating point."""
    ho = rom.ho[effective_row(turb, rom.offset(water, dist))]
    cls = rom.hs_class(water, dist)
    if cls == 0:
        return ho.copy()                                # h_s == 256
    hs = rom.hs[cls]
    return np.minimum((hs[:, None] * ho[None, :]) >> 8, H_MAX).ravel()


def psig_codes(rom: Rom, lam: int, water: int, dist: int, h_f: np.ndarray,
               mu_level: int = 8, mu_gain_log2: int = 0,
               shift: int | None = None) -> np.ndarray:
    ref = rom.psig_ref(lam, water, dist)
    psig_mu = min(((ref * mu_level) << mu_gain_log2) >> 3, P_MAX)
    p = np.minimum((psig_mu * np.asarray(h_f, dtype=np.int64)) >> 8, P_MAX)
    if shift is not None:
        p = p >> shift
    return p


def event_probs(p_sig_code, p_noise_code: int, epol_code: int):
    """Per-attempt click probability and click-and-error probability."""
    ps = np.asarray(p_sig_code, dtype=float) / S24
    pn = p_noise_code / S24
    ep = epol_code / S24
    click = 1.0 - (1.0 - ps) * (1.0 - pn)
    err = ps * ep + (1.0 - ps) * pn * 0.5
    return click, err


def mu_of(mu_level: int = 8, mu_gain_log2: int = 0, shift: int | None = None) -> float:
    mu = MU_NOM * mu_level / 8.0 * (1 << mu_gain_log2)
    return mu / (1 << shift) if shift is not None else mu


# ===========================================================================
# Security functional
# ===========================================================================

def h2(e):
    e = np.asarray(e, dtype=float)
    out = np.zeros_like(e)
    m = (e > 0.0) & (e < 1.0)
    x = e[m]
    out[m] = -(x * np.log2(x) + (1.0 - x) * np.log2(1.0 - x))
    return out


def V(G, T):
    """Reduced privacy functional G[1 - H2(T/G)], unclipped on [0, 1]."""
    G = np.asarray(G, dtype=float)
    T = np.asarray(T, dtype=float)
    e = np.divide(T, G, out=np.zeros_like(T), where=G > 0)
    return G * (1.0 - h2(e))


def defect(G, T, w):
    w = np.asarray(w, dtype=float)
    w = w / w.sum()
    G0, T0 = float((w * G).sum()), float((w * T).sum())
    V0 = float(V(np.array([G0]), np.array([T0]))[0])
    return float((w * V(G, T)).sum()) - V0, V0, G0, T0


def bkm(G, T, w):
    w = np.asarray(w, dtype=float)
    w = w / w.sum()
    G0, T0 = float((w * G).sum()), float((w * T).sum())
    e0 = T0 / G0
    k = 1.0 / (LN2 * e0 * (1.0 - e0) * G0)
    dG, dT = np.asarray(G) - G0, np.asarray(T) - T0
    return 0.5 * float(np.sum(w * k * (e0 * e0 * dG * dG - 2 * e0 * dG * dT + dT * dT)))


def quantile_labels(values, q) -> np.ndarray:
    return pd.qcut(pd.Series(np.asarray(values, dtype=float)), q=q, labels=False,
                   duplicates="drop").to_numpy().astype(int)


def binned_moments(values, G, T, q):
    lab = quantile_labels(values, q)
    nb = lab.max() + 1
    cnt = np.bincount(lab, minlength=nb).astype(float)
    return (np.bincount(lab, weights=G, minlength=nb) / cnt,
            np.bincount(lab, weights=T, minlength=nb) / cnt, cnt)


def point_defects(rom: Rom, water: int = 0, dist: int = 4, turb: int = 5,
                  lam: int = 0, mu_level: int = 8, mu_gain_log2: int = 0,
                  p_noise: int = Y0_CODE, bins: int = 8, q_sift: float = 0.5) -> dict:
    """Population defect of one operating point at 8 bins and full resolution."""
    hf = fading_values(rom, water, dist, turb)
    ps = psig_codes(rom, lam, water, dist, hf, mu_level, mu_gain_log2)
    click, err = event_probs(ps, p_noise, rom.epol_code(lam, water, dist))
    G, T = q_sift * click, q_sift * err
    w = np.ones(len(hf))
    dfull, V0, G0, T0 = defect(G, T, w)
    # Full resolution groups equal h_f values: identical population statistics.
    Gb, Tb, cnt = binned_moments(hf, G, T, bins)
    d8 = defect(Gb, Tb, cnt)[0]
    b8 = bkm(Gb, Tb, cnt)
    hv = hf / 256.0
    return dict(sigma2_realized=float(hv.var() / hv.mean() ** 2),
                P_h_below_hstar=float(np.mean(ps < p_noise)) if p_noise else 0.0,
                d_bins_over_V=d8 / V0, d_full_over_V=dfull / V0,
                bkm_bins_over_V=b8 / V0,
                bkm_ratio_bins=(b8 / d8) if d8 > 0 else float("nan"),
                d_bins_bits=d8 / G0, d_full_bits=dfull / G0,
                QBER=T0 / G0, G0=G0, V0=V0, bins=bins)


# ===========================================================================
# Block-level Monte Carlo of the counting statistics (planning, coverage)
# ===========================================================================

def simulate_blocks(rom: Rom, n_blocks: int, rng, water: int = 0, dist: int = 4,
                    turb: int = 5, lam: int = 0, mu_level: int = 8,
                    mu_gain_log2: int = 0, p_noise: int = Y0_CODE,
                    classes: tuple = ((1.0, None),)) -> pd.DataFrame:
    """Blocks with binomial counts, one fading draw per block.

    classes: sequence of (probability, shift) with shift None = signal and
    shift "vac" = vacuum. Returns per-block counts per class, basis-resolved,
    with h_mean in 12-bit units like the capture CSV."""
    hf_all = fading_values(rom, water, dist, turb)
    hf = hf_all[rng.integers(0, len(hf_all), n_blocks)]
    epol = rom.epol_code(lam, water, dist)
    out = {"h_mean": hf, "attempts": np.full(n_blocks, ATTEMPTS_PER_BLOCK)}
    remaining = np.full(n_blocks, ATTEMPTS_PER_BLOCK)
    probs = [p for p, _ in classes]
    for k, (p, shift) in enumerate(classes):
        cond = p / max(1e-300, 1.0 - sum(probs[:k]))
        n_k = remaining if k == len(classes) - 1 else rng.binomial(remaining, min(1.0, cond))
        remaining = remaining - n_k
        if shift == "vac":
            ps = np.zeros(n_blocks, dtype=np.int64)
        else:
            ps = psig_codes(rom, lam, water, dist, hf, mu_level, mu_gain_log2, shift)
        click, err = event_probs(ps, p_noise, epol)
        e_given = np.divide(err, click, out=np.zeros_like(err), where=click > 0)
        tag = "c%d_" % k
        out[tag + "attempts"] = n_k
        tz = rng.binomial(n_k, 0.25)
        tx = rng.binomial(n_k - tz, 1.0 / 3.0)
        cz = rng.binomial(tz, click)
        cx = rng.binomial(tx, click)
        out[tag + "trials_z"], out[tag + "trials_x"] = tz, tx
        out[tag + "clicks_z"], out[tag + "clicks_x"] = cz, cx
        out[tag + "errors_z"] = rng.binomial(cz, e_given)
        out[tag + "errors_x"] = rng.binomial(cx, e_given)
        # Mismatched-basis attempts click independently of the sifted ones.
        out[tag + "clicks_all"] = cz + cx + rng.binomial(n_k - tz - tx, click)
    return pd.DataFrame(out)


if __name__ == "__main__":
    rom = load_rom()
    rows = []
    for turb in range(1, 8):
        r = point_defects(rom, turb=turb)
        r["turb"] = turb
        rows.append(r)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(pd.DataFrame(rows).set_index("turb").to_string(float_format=lambda v: "%.4g" % v))
