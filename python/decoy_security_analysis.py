#!/usr/bin/env python3
"""
decoy_security_analysis.py -- block-resolved versus pooled decoy-state key rates.

INPUT
  v6 captures (capture_security_blocks.py rev E, rev B bitstream) with decoy
  intensities enabled. Every block carries basis-resolved trials / clicks /
  errors for the signal (class 0 = total - decoy - vacuum), decoy (class 1) and
  vacuum (class 2) attempts, all interleaved attempt by attempt.

BOUNDS
  Vacuum + weak decoy, Ma, Qi, Zhao and Lo, Phys. Rev. A 72, 012326 (2005):
      Y1L = mu/(mu nu - nu^2) [Q_nu e^nu - Q_mu e^mu nu^2/mu^2 - (mu^2-nu^2)/mu^2 Y0]
      T1U = (E_nu Q_nu e^nu - e0 Y0) / nu,          e1U = T1U / Y1L
  Key from basis b, single-photon phase error from the other basis b':
      R_b = mu e^-mu Y1L^b (1 - H2(e1U^b')) - f Q_mu^b H2(E_mu^b)
  per basis-matched signal attempt, and per attempt of the capture
      R = sum_b (signal trials in b / attempts) R_b .
  Y1L and T1U are affine in the measured gains, and R is the sum of a convex
  perspective and the negative of a concave one, so the coarse-graining
  inequality of Proposition 1 carries over:  sum_j p_j R(m_j) >= R(sum_j p_j m_j).

MODES
  asymptotic  point estimates of every gain.
  finite      every gain replaced by its worst-case Clopper-Pearson bound; the
              failure probability is split over all bounds of all bins (union
              bound). Counting within a bin is treated as i.i.d. given the
              channel -- this is a statistical-fluctuation analysis, NOT a
              composable finite-key proof.

LABELS
  h     the latent on-chip fading word (upper bound on what any label recovers)
  obs   clicks of basis-mismatched attempts of the block, all classes: announced
        by Bob, public, and independent of the sifted counts given the fading.

MODEL
  --model simulates captures from rom_exact_model for a planned configuration,
  to choose mu, the decoy ratio, the class preset and the block budget before
  touching the board. Its output is a prediction and must be labelled as such.

USAGE
  python decoy_security_analysis.py --csv data/q2b/decoy_L1.csv data/q2b/decoy_L3.csv ...
  python decoy_security_analysis.py --model --turb 1 3 5 --blocks 40000 \
      --mu-gain 2 --shift 3 --preset 0
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rom_exact_model as rm                         # noqa: E402

OUT_DIR = Path("phase2_results")
SEED = 20260915
F_EC = 1.16
E0 = 0.5
EPS_TOTAL = 1e-9
F_REP = 1e7
BASES = ("z", "x")
CLASS_FIELDS = ("attempts", "clicks_all", "trials_z", "clicks_z", "errors_z",
                "trials_x", "clicks_x", "errors_x")
PRESETS = {0: (0.50, 0.25, 0.25), 1: (0.70, 0.20, 0.10),
           2: (0.80, 0.15, 0.05), 3: (1 / 3, 1 / 3, 1 / 3)}
GROUP_COLS = ["water", "dist_idx", "turb", "lambda_sel", "mu_level",
              "mu_gain_log2", "decoy_shift", "decoy_preset", "p_noise"]
COUNT_KEYS = ["attempts_total"] + ["%s_%s" % (c, f) for c in "sdv"
                                   for f in ("trials_z", "clicks_z", "errors_z",
                                             "trials_x", "clicks_x", "errors_x")]


def h2s(x):
    x = float(x)
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -(x * math.log2(x) + (1.0 - x) * math.log2(1.0 - x))


def cp_lo(k, n, eps):
    return 0.0 if k <= 0 or n <= 0 else float(beta.ppf(eps, k, n - k + 1))


def cp_hi(k, n, eps):
    return 1.0 if n <= 0 or k >= n else float(beta.ppf(1.0 - eps, k + 1, n - k))


# ===========================================================================
# Key rate of one set of summed counts
# ===========================================================================

def key_rate(cnt, mu, nu, finite=False, eps=1e-10, f_ec=F_EC):
    """Per-attempt asymptotic or fluctuation-aware decoy-state key rate.

    cnt holds summed counts: attempts_total and, for class c in s/d/v and basis
    b in z/x, c_trials_b, c_clicks_b, c_errors_b. Returns (R, R_clipped, parts)
    where R_clipped discards a basis contribution that is negative."""
    em, en = math.exp(mu), math.exp(nu)
    per = {}
    for b in BASES:
        ts, cs, es = cnt["s_trials_" + b], cnt["s_clicks_" + b], cnt["s_errors_" + b]
        td, cd, ed = cnt["d_trials_" + b], cnt["d_clicks_" + b], cnt["d_errors_" + b]
        tv, cv = cnt["v_trials_" + b], cnt["v_clicks_" + b]
        if min(ts, td, tv) <= 0 or cs <= 0:
            return float("nan"), float("nan"), {}
        Qs, EQs = cs / ts, es / ts
        if finite:
            Qs_hi = cp_hi(cs, ts, eps)
            Qd_lo = cp_lo(cd, td, eps)
            EQd_hi = cp_hi(ed, td, eps)
            Y0_hi, Y0_lo = cp_hi(cv, tv, eps), cp_lo(cv, tv, eps)
        else:
            Qs_hi, Qd_lo, EQd_hi = Qs, cd / td, ed / td
            Y0_hi = Y0_lo = cv / tv
        Y1L = mu / (mu * nu - nu * nu) * (Qd_lo * en - Qs_hi * em * nu * nu / (mu * mu)
                                          - (mu * mu - nu * nu) / (mu * mu) * Y0_hi)
        T1U = (EQd_hi * en - E0 * Y0_lo) / nu
        per[b] = {"Qs": Qs, "EQs": EQs, "Y1L": max(Y1L, 0.0), "T1U": max(T1U, 0.0),
                  "trials": ts}
    R = Rc = 0.0
    parts = {}
    for b, bp in (("z", "x"), ("x", "z")):
        kb, pb = per[b], per[bp]
        e1U = min(pb["T1U"] / pb["Y1L"], 0.5) if pb["Y1L"] > 0 else 0.5
        priv = mu * math.exp(-mu) * kb["Y1L"] * (1.0 - h2s(e1U))
        ec = f_ec * kb["Qs"] * h2s(kb["EQs"] / kb["Qs"])
        rb = (priv - ec) * kb["trials"] / cnt["attempts_total"]
        R += rb
        Rc += max(rb, 0.0)
        parts["Y1L_" + b], parts["e1U_" + b], parts["R_" + b] = kb["Y1L"], e1U, rb
    return R, Rc, parts


def class_counts(df):
    """Per-block counts per class in the key_rate naming."""
    out = pd.DataFrame({"attempts_total": df["attempts"].astype(float)})
    for f in ("trials_z", "clicks_z", "errors_z", "trials_x", "clicks_x", "errors_x"):
        out["d_" + f] = df["dec_" + f].astype(float)
        out["v_" + f] = df["vac_" + f].astype(float)
        out["s_" + f] = (df[f] - df["dec_" + f] - df["vac_" + f]).astype(float)
    out["h"] = df["h_mean"].to_numpy()
    out["obs"] = (df["clicks_all"] - df["clicks_z"] - df["clicks_x"]).to_numpy()
    return out


def _labels(values, q):
    if q == 1:
        return np.zeros(len(values), dtype=np.int64)
    return pd.qcut(pd.Series(values), q=q, labels=False,
                   duplicates="drop").to_numpy().astype(np.int64)


def binned_rates(cc, label_col, q, mu, nu, finite, idx=None):
    lab = _labels(cc[label_col].to_numpy(), q)
    if idx is not None:
        lab = lab[idx]
    nb = int(lab.max()) + 1
    sums = {k: np.bincount(lab, weights=(cc[k].to_numpy() if idx is None
                                         else cc[k].to_numpy()[idx]), minlength=nb)
            for k in COUNT_KEYS}
    n_all = sums["attempts_total"].sum()
    eps = EPS_TOTAL / (12.0 * nb)
    R = Rc = 0.0
    for b in range(nb):
        cnt = {k: sums[k][b] for k in COUNT_KEYS}
        if cnt["attempts_total"] <= 0:
            continue
        r, rc, _ = key_rate(cnt, mu, nu, finite, eps)
        w = cnt["attempts_total"] / n_all
        if np.isnan(r):
            continue
        R += w * r
        Rc += w * rc
    return R, Rc, nb, EPS_TOTAL


def analyse_capture(df, boot_reps, rng, bins=(1, 2, 4, 8, 16, 32)):
    mu = rm.MU_NOM * float(df["mu_level"].iloc[0]) / 8.0 * (1 << int(df["mu_gain_log2"].iloc[0]))
    nu = mu / (1 << int(df["decoy_shift"].iloc[0]))
    cc = class_counts(df)
    rows = []
    for label_col in ("h", "obs"):
        for q in bins:
            if label_col == "obs" and q == 1:
                continue
            for finite in (False, True):
                R, Rc, nb, eps_tot = binned_rates(cc, label_col, q, mu, nu, finite)
                lo = hi = np.nan
                if boot_reps and not finite:
                    n = len(cc)
                    vals = []
                    for _ in range(boot_reps):
                        idx = rng.integers(0, n, n)
                        vals.append(binned_rates(cc, label_col, q, mu, nu, False, idx)[1])
                    lo, hi = np.percentile(vals, [2.5, 97.5])
                rows.append({"label": label_col, "bins": q, "n_bins": nb,
                             "mode": "finite" if finite else "asymptotic",
                             "mu": mu, "nu": nu, "R_per_attempt": R,
                             "R_clipped_per_attempt": Rc,
                             "R_clipped_bit_per_s": Rc * F_REP,
                             "ci_lo_clipped": lo, "ci_hi_clipped": hi,
                             "eps_total": eps_tot if finite else np.nan})
    out = pd.DataFrame(rows)
    pooled = out[(out.bins == 1)].set_index("mode")["R_clipped_per_attempt"]
    out["gain_over_pooled"] = [r / pooled.get(m, np.nan) if pooled.get(m, 0) > 0 else np.nan
                               for r, m in zip(out["R_clipped_per_attempt"], out["mode"])]
    return out


# ===========================================================================
# Measured captures
# ===========================================================================

def run_measured(args):
    OUT_DIR.mkdir(exist_ok=True)
    frames = [pd.read_csv(p, low_memory=False) for p in args.csv]
    df = pd.concat(frames, ignore_index=True)
    df = df[(df.get("capture_valid", 1) == 1) & (df["decoy_enable"] == 1)].copy()
    if df.empty:
        raise SystemExit("no valid decoy-enabled v6 blocks in the input")
    rows = []
    for gi, (key, g) in enumerate(df.groupby(GROUP_COLS, sort=True)):
        rng = np.random.default_rng(SEED + gi)
        res = analyse_capture(g, args.bootstrap, rng)
        for c, v in zip(GROUP_COLS, key):
            res[c] = v
        res["n_blocks"] = len(g)
        rows.append(res)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT_DIR / "decoy_rates.csv", index=False)
    _print(out)


# ===========================================================================
# Model prediction
# ===========================================================================

def synthetic_capture(rom, turb, blocks, args, rng):
    ps, pd_, pv = PRESETS[args.preset]
    sim = rm.simulate_blocks(rom, blocks, rng, turb=turb, mu_level=args.mu_level,
                             mu_gain_log2=args.mu_gain_log2, p_noise=args.p_noise,
                             classes=((ps, None), (pd_, args.shift), (pv, "vac")))
    df = pd.DataFrame({"h_mean": sim["h_mean"], "attempts": sim["attempts"]})
    for f in CLASS_FIELDS:
        if f == "attempts":
            df["dec_attempts"], df["vac_attempts"] = sim["c1_attempts"], sim["c2_attempts"]
            continue
        df[f] = sim["c0_" + f] + sim["c1_" + f] + sim["c2_" + f]
        df["dec_" + f], df["vac_" + f] = sim["c1_" + f], sim["c2_" + f]
    df["mu_level"], df["mu_gain_log2"], df["decoy_shift"] = args.mu_level, args.mu_gain_log2, args.shift
    return df


def run_model(args):
    OUT_DIR.mkdir(exist_ok=True)
    rom = rm.load_rom()
    rows = []
    for turb in args.turb:
        for rep in range(args.reps):
            rng = np.random.default_rng(SEED + 100 * turb + rep)
            res = analyse_capture(synthetic_capture(rom, turb, args.blocks, args, rng), 0, rng)
            res["turb"], res["rep"] = turb, rep
            rows.append(res)
    out = pd.concat(rows, ignore_index=True)
    agg = (out.groupby(["turb", "label", "bins", "mode"])
              .agg(mu=("mu", "first"), nu=("nu", "first"),
                   R_clipped=("R_clipped_per_attempt", "mean"),
                   R_clipped_sd=("R_clipped_per_attempt", "std"),
                   R_unclipped=("R_per_attempt", "mean"),
                   gain=("gain_over_pooled", "mean"))
              .reset_index())
    agg["bit_per_s"] = agg["R_clipped"] * F_REP
    tag = "decoy_model_m%d_g%d_s%d_p%d" % (args.mu_level, args.mu_gain_log2, args.shift, args.preset)
    agg.to_csv(OUT_DIR / (tag + ".csv"), index=False)
    print("MODEL PREDICTION (ROM-exact, %d blocks x %d reps per level), mu=%.4g nu=%.4g preset %d"
          % (args.blocks, args.reps, agg["mu"].iloc[0], agg["nu"].iloc[0], args.preset))
    _print(agg.rename(columns={"R_clipped": "R_clipped_per_attempt"}))


def _print(tab):
    cols = [c for c in ("turb", "label", "bins", "mode", "R_clipped_per_attempt",
                        "R_clipped_sd", "ci_lo_clipped", "ci_hi_clipped",
                        "bit_per_s", "R_clipped_bit_per_s", "gain", "gain_over_pooled")
            if c in tab.columns]
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(tab[cols].to_string(index=False, float_format=lambda v: "%.4g" % v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs="*", default=[])
    ap.add_argument("--bootstrap", type=int, default=200)
    ap.add_argument("--model", action="store_true")
    ap.add_argument("--turb", type=int, nargs="*", default=[1, 3, 5])
    ap.add_argument("--blocks", type=int, default=40000)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--mu-level", type=int, default=8)
    ap.add_argument("--mu-gain-log2", "--mu-gain", type=int, default=2, dest="mu_gain_log2")
    ap.add_argument("--shift", type=int, default=3)
    ap.add_argument("--preset", type=int, default=0, choices=range(4))
    ap.add_argument("--p-noise", type=int, default=rm.Y0_CODE)
    args = ap.parse_args()
    if args.model:
        run_model(args)
    elif args.csv:
        run_measured(args)
    else:
        ap.error("give --csv captures or --model")


if __name__ == "__main__":
    main()
