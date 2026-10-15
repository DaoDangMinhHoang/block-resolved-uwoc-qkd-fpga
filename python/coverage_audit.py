#!/usr/bin/env python3
"""
coverage_audit.py -- does the phase-2 estimator deliver what it claims?

Synthetic captures are drawn from the ROM-exact population (rom_exact_model),
so the estimand is known exactly: the population defect at the same number of
equal-population bins. Every synthetic capture goes through the estimator used
in the paper:

    8 equal-population bins of h -> exact defect
    -> parametric null (counts redrawn at the pooled operating point)
    -> within-bin block bootstrap, each replicate paired with a null replicate
    -> 95% percentile interval; "resolved" when its lower end is above zero

For every regime the audit reports the mean bias of the null-subtracted
estimate, the empirical coverage of the 95% interval, how often the level is
declared resolved, and, for Y0 = 0 regimes, the size of the G-test at 5%.
Null and bootstrap use --reps-inner replicates (the paper uses 600).

    python coverage_audit.py                    # all regimes, ~20 min
    python coverage_audit.py --quick            # fewer replicates
"""

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rom_exact_model as rm                         # noqa: E402

OUT_DIR = Path("phase2_results")
SEED = 20260916
LN2 = math.log(2.0)
BINS = 8

# (name, turbulence row, Y0 code, blocks, outer replicates)
REGIMES = [
    ("L3 null (Y0=0)", 3, 0, 20000, 150),
    ("L5 null (Y0=0)", 5, 0, 20000, 150),
    ("L2", 2, 218, 20000, 150),
    ("L3", 3, 218, 20000, 150),
    ("L4", 4, 218, 20000, 150),
    ("L5", 5, 218, 20000, 150),
    ("L1 null (Y0=0)", 1, 0, 286000, 30),
    ("L1", 1, 218, 286000, 30),
]


def h2(e):
    e = np.asarray(e, dtype=float)
    out = np.zeros_like(e)
    m = (e > 0.0) & (e < 1.0)
    x = e[m]
    out[m] = -(x * np.log2(x) + (1.0 - x) * np.log2(1.0 - x))
    return out


def defect_rows(Nb, Cb, Eb):
    """Exact defect and pooled V for arrays shaped (reps, bins)."""
    G = np.divide(Cb, Nb, out=np.zeros_like(Cb), where=Nb > 0)
    T = np.divide(Eb, Nb, out=np.zeros_like(Eb), where=Nb > 0)
    w = Nb / Nb.sum(axis=1, keepdims=True)
    e = np.divide(T, G, out=np.zeros_like(T), where=G > 0)
    Vb = G * (1.0 - h2(e))
    G0, T0 = (w * G).sum(axis=1), (w * T).sum(axis=1)
    V0 = G0 * (1.0 - h2(T0 / G0))
    return (w * Vb).sum(axis=1) - V0, V0


def bin_sums(lab, nb, N, C, E):
    return (np.bincount(lab, weights=N, minlength=nb)[None, :],
            np.bincount(lab, weights=C, minlength=nb)[None, :],
            np.bincount(lab, weights=E, minlength=nb)[None, :])


def null_replicates(lab, nb, N, A, C, E, reps, rng, chunk=20):
    pA, g, e0 = A.sum() / N.sum(), C.sum() / A.sum(), E.sum() / C.sum()
    out = []
    for start in range(0, reps, chunk):
        r = min(chunk, reps - start)
        Ai = rng.binomial(np.broadcast_to(N.astype(np.int64), (r, len(N))), pA)
        Ci = rng.binomial(Ai, g)
        Ei = rng.binomial(Ci, e0)
        flat = (np.arange(r)[:, None] * nb + lab[None, :]).ravel()
        Nb = np.bincount(flat, weights=np.broadcast_to(N, (r, len(N))).ravel(),
                         minlength=r * nb).reshape(r, nb)
        Cb = np.bincount(flat, weights=Ci.ravel(), minlength=r * nb).reshape(r, nb)
        Eb = np.bincount(flat, weights=Ei.ravel(), minlength=r * nb).reshape(r, nb)
        out.append(defect_rows(Nb, Cb.astype(float), Eb.astype(float))[0])
    return np.concatenate(out)


def bootstrap_replicates(lab, nb, N, C, E, reps, rng):
    Nb = np.zeros((reps, nb))
    Cb = np.zeros((reps, nb))
    Eb = np.zeros((reps, nb))
    for b in range(nb):
        members = np.flatnonzero(lab == b)
        if len(members) == 0:
            continue
        idx = members[rng.integers(0, len(members), (reps, len(members)))]
        Nb[:, b], Cb[:, b], Eb[:, b] = N[idx].sum(1), C[idx].sum(1), E[idx].sum(1)
    return defect_rows(Nb, Cb, Eb)[0]


def one_capture(rom, turb, p_noise, blocks, reps_inner, rng):
    sim = rm.simulate_blocks(rom, blocks, rng, turb=turb, p_noise=p_noise)
    N = sim["c0_attempts"].to_numpy(float)
    A = (sim["c0_trials_z"] + sim["c0_trials_x"]).to_numpy(float)
    C = (sim["c0_clicks_z"] + sim["c0_clicks_x"]).to_numpy(float)
    E = (sim["c0_errors_z"] + sim["c0_errors_x"]).to_numpy(float)
    lab = rm.quantile_labels(sim["h_mean"].to_numpy(), BINS)
    nb = int(lab.max()) + 1
    d_obs, V_obs = defect_rows(*bin_sums(lab, nb, N, C, E))
    d_obs, V_obs = float(d_obs[0]), float(V_obs[0])
    dnull = null_replicates(lab, nb, N, A, C, E, reps_inner, rng)
    dboot = bootstrap_replicates(lab, nb, N, C, E, reps_inner, rng)
    corr = dboot - dnull[rng.integers(0, len(dnull), len(dboot))]
    lo, hi = np.percentile(corr, [2.5, 97.5])
    # G-test on the same table
    Cb = np.bincount(lab, weights=C, minlength=nb)
    Eb = np.bincount(lab, weights=E, minlength=nb)
    m = Cb > 0
    Ct = Cb[m].sum()
    e0 = Eb[m].sum() / Ct
    I = float(h2(np.array([e0]))[0] - np.sum(Cb[m] / Ct * h2(Eb[m] / Cb[m])))
    G = 2.0 * Ct * LN2 * I
    mm = (m.sum() - 1) / (2.0 * Ct * LN2)
    phi0 = 1.0 - float(h2(np.array([e0]))[0])
    return {"excess_V": (d_obs - dnull.mean()) / V_obs, "lo_V": lo / V_obs,
            "hi_V": hi / V_obs, "mm_excess_V": (I - mm) / phi0,
            "p_gtest": float(chi2.sf(G, m.sum() - 1))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps-inner", type=int, default=200)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    OUT_DIR.mkdir(exist_ok=True)
    rom = rm.load_rom()
    rows = []
    for gi, (name, turb, p_noise, blocks, reps) in enumerate(REGIMES):
        if args.quick:
            reps = max(10, reps // 5)
        truth = 0.0 if p_noise == 0 else rm.point_defects(rom, turb=turb, bins=BINS)["d_bins_over_V"]
        rng = np.random.default_rng(SEED + gi)
        t0 = time.time()
        res = pd.DataFrame([one_capture(rom, turb, p_noise, blocks, args.reps_inner, rng)
                            for _ in range(reps)])
        cover = float(np.mean((res.lo_V <= truth) & (truth <= res.hi_V)))
        row = {"regime": name, "blocks": blocks, "replicates": reps,
               "inner_reps": args.reps_inner, "true_over_V": truth,
               "mean_estimate_over_V": float(res.excess_V.mean()),
               "bias_over_V": float(res.excess_V.mean() - truth),
               "sd_estimate_over_V": float(res.excess_V.std()),
               "mean_ci_width_over_V": float((res.hi_V - res.lo_V).mean()),
               "coverage95": cover,
               "resolved_rate": float(np.mean(res.lo_V > 0)),
               "mm_bias_over_V": float(res.mm_excess_V.mean() - truth),
               "gtest_reject_5pct": float(np.mean(res.p_gtest < 0.05)),
               "seconds": round(time.time() - t0, 1)}
        rows.append(row)
        print("%-16s coverage=%.3f resolved=%.3f bias=%+.2e (truth %.3e) G-test 5%%: %.3f  [%ss]"
              % (name, cover, row["resolved_rate"], row["bias_over_V"], truth,
                 row["gtest_reject_5pct"], row["seconds"]), flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "coverage_audit.csv", index=False)
    with pd.option_context("display.width", 250):
        print(out.to_string(index=False, float_format=lambda v: "%.4g" % v))


if __name__ == "__main__":
    main()
