#!/usr/bin/env python3
"""
defect_from_clicklog.py (rev A)

Technical Assessment section 10 item 2, "Strongly recommended": compute the
actual coarse-graining defect for turbulence levels L1-L5 from the raw block
data that ALREADY EXISTS, without repeating the 34-hour campaign. Item 3,
"Recommended if data permit", is the exact-versus-BKM comparison, produced by
the same run.

WHY THIS SCRIPT EXISTS SEPARATELY FROM phase2_security_analysis.py

  phase2_security_analysis.py consumes the 96-byte binary security blocks
  emitted by block_reporter_bin.v. No such capture exists yet; it needs the DE1
  board. What does exist is data/clicks_fixed_B_turb_L*.csv, one line per
  detected click, collected with --coh 11 so that one coherence block equals one
  monitor window of exactly 65536 attempts and carries exactly one frozen
  fading sample. That is already block-resolved data in the sense of the
  assessment, section 7. This script reconstructs the blockwise observation
  vectors from it and then calls the SAME estimator, so there is one
  implementation of the statistics, not two.

WHAT IS RECONSTRUCTED PER BLOCK j

  block index      attempt >> (coh + 5)        (coh = 11 -> 65536 attempts)
  attempts N_j     2^(coh+5), by construction of the campaign
  clicks   C_j     rows with bmatch == 1       (sifted detections)
  errors   E_j     rows with bmatch == 1 and err == 1
  h_j              irrad / 128                 (the log scale is 128 = 1.0)

  G_j = C_j / N_j and T_j = E_j / N_j are per PHYSICAL ATTEMPT, and the mixture
  weight is N_j, exactly as Proposition 1 requires.

FOUR LIMITATIONS, PRINTED AT THE END OF EVERY RUN

  1. 160 blocks per level. The planning table says L4 needs ~2.9e3 and L3
     ~2.0e4, so L1-L3 are expected NOT to resolve. That is a sample-size
     statement, not a null result.
  2. irrad saturates at 255 (h = 2.0), so the top bin is compressed at L5. The
     security path uses the 12-bit h_mean (256 = 1.0, ceiling 16.0) and does
     not have this problem.
  3. A block with no sifted click has no observable h. Those blocks are the
     deep fades that carry most of the defect. --drop-empty measures how much
     that choice matters; it is a real selection bias either way.
  4. trials (basis-matched attempts) are not in the per-click log, since lost
     attempts are never logged. They are set to the nominal N/2. This does not
     touch the per-attempt estimator; it only makes the parametric null
     slightly less accurate.

USAGE

  python python/defect_from_clicklog.py
  python python/defect_from_clicklog.py --bins 8 --robust
  python python/defect_from_clicklog.py --drop-empty --no-figures
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from phase2_security_analysis import (      # noqa: E402
    SEED, SIGMA2_HO_20M, analyse_group, moments_from_h, operating_point,
    sample_h, V_priv,
)

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "phase2_results"
FIG_DIR = ROOT / "Images"
DEFAULT_GLOB = "data/clicks_fixed_B_turb_L*.csv"
LEVEL_RE = re.compile(r"_L(\d+)\.csv$", re.IGNORECASE)


# ===========================================================================
# 1. Rebuild blocks from the per-click log
# ===========================================================================

def blocks_from_clicklog(path, coh, drop_empty):
    """One row per coherence block, in the column names analyse_group expects."""
    d = pd.read_csv(path)
    need = {"attempt", "bmatch", "err", "irrad"}
    missing = need - set(d.columns)
    if missing:
        raise SystemExit("%s is missing column(s): %s"
                         % (path, ", ".join(sorted(missing))))

    block_len = 1 << (coh + 5)
    d["blk"] = d["attempt"].to_numpy(dtype=np.int64) >> (coh + 5)

    sift = d["bmatch"].to_numpy(dtype=np.int64) == 1
    d["is_sift"] = sift.astype(np.int64)
    d["is_err"] = (sift & (d["err"].to_numpy(dtype=np.int64) == 1)).astype(np.int64)

    g = d.groupby("blk", sort=True)
    rec = pd.DataFrame({
        "clicks": g["is_sift"].sum(),
        "errors": g["is_err"].sum(),
        "clicks_all": g.size(),
        "irrad_n": g["irrad"].nunique(),
        "irrad": g["irrad"].median(),
    }).reset_index()

    # The campaign ran a fixed 160 x 65536 attempts per level, so the blocks
    # that produced no click at all still happened and still carry weight.
    # Reindexing over the full span puts them back with C = E = 0.
    n_blocks = int(np.ceil((d["attempt"].max() + 1) / block_len))
    full = pd.DataFrame({"blk": np.arange(n_blocks, dtype=np.int64)})
    rec = full.merge(rec, on="blk", how="left").fillna(
        {"clicks": 0, "errors": 0, "clicks_all": 0, "irrad_n": 0, "irrad": 0.0})

    empty = int((rec["clicks_all"] == 0).sum())
    split = int((rec["irrad_n"] > 1).sum())
    if drop_empty:
        rec = rec[rec["clicks_all"] > 0].copy()

    rec["attempts"] = block_len
    rec["trials"] = block_len // 2          # nominal, see limitation 4
    rec["h_q8"] = rec["irrad"] / 128.0      # the log scale is 128 = 1.0
    rec = rec[["attempts", "trials", "clicks", "errors", "h_q8", "clicks_all"]]
    rec = rec.astype({"attempts": float, "trials": float,
                      "clicks": float, "errors": float})
    return rec, {"n_blocks": n_blocks, "empty": empty, "split": split,
                 "kept": len(rec), "block_len": block_len}


def level_of(path):
    m = LEVEL_RE.search(Path(path).name)
    return int(m.group(1)) if m else -1


# ===========================================================================
# 2. Model prediction for the same levels
# ===========================================================================

def model_curve(distance, n_draw=400000):
    """true_over_V per level, the same quantity phase2 --model-only prints."""
    rng = np.random.default_rng(SEED)
    op = operating_point(distance, 450)
    out = {}
    for lvl, s2 in SIGMA2_HO_20M.items():
        h = sample_h(s2, n_draw, rng)
        G, T = moments_from_h(op, h)
        V0 = float(V_priv(np.array([G.mean()]), np.array([T.mean()]))[0])
        out[lvl] = (float(V_priv(G, T).mean()) - V0) / V0
    return out


# ===========================================================================
# 3. Figures
# ===========================================================================

def fig_defect_vs_turbulence(tab, model, path):
    fig, ax = plt.subplots(figsize=(6.8, 4.8))

    res = tab[tab["resolved"] == 1]
    unres = tab[tab["resolved"] == 0]

    if len(res):
        yerr = np.abs(np.vstack([res["excess_over_V"] - res["ci_lo_over_V"],
                                 res["ci_hi_over_V"] - res["excess_over_V"]]))
        ax.errorbar(res["level"], res["excess_over_V"], yerr=yerr,
                    marker="o", capsize=3, ls="-", color="tab:blue",
                    label=r"measured $\Delta_{cg}/V$, resolved")
    if len(unres):
        # Proposition 1 gives Delta >= 0. An unresolved level is a 95% upper
        # limit, never a negative measurement.
        ax.plot(unres["level"], unres["upper_limit95_over_V"],
                marker="v", ls="none", mfc="none", color="tab:blue",
                label="95% upper limit (sample-size limited)")

    lv = sorted(model)
    ax.plot(lv, [model[k] for k in lv], "k--", lw=1.2, marker="s", ms=4,
            label="analytical channel model")

    ax.set_yscale("log")
    ax.set_xticks(lv)
    ax.set_xlabel("ocean turbulence level")
    ax.set_ylabel(r"$\Delta_{cg}\,/\,V$, counting-noise null subtracted")
    ax.set_title("Coarse-graining defect versus turbulence, clear ocean 25 m")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def fig_exact_vs_bkm(tab, path):
    sel = tab[(tab["excess_over_V"] > 0) & (tab["bkm_excess_over_V"] > 0)]
    if len(sel) < 2:
        return False
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.errorbar(sel["level"], sel["excess_over_V"],
                yerr=np.abs(np.vstack([
                    sel["excess_over_V"] - sel["ci_lo_over_V"],
                    sel["ci_hi_over_V"] - sel["excess_over_V"]])),
                marker="o", capsize=3, label="exact defect")
    ax.errorbar(sel["level"], sel["bkm_excess_over_V"],
                yerr=np.abs(np.vstack([
                    sel["bkm_excess_over_V"] - sel["bkm_ci_lo_over_V"],
                    sel["bkm_ci_hi_over_V"] - sel["bkm_excess_over_V"]])),
                marker="s", ls="--", capsize=3,
                label=r"BKM $\frac{1}{2}\mathrm{Tr}(\Gamma_0\mathrm{Cov})$")
    ax.set_yscale("log")
    ax.set_xticks(sorted(sel["level"]))
    ax.set_xlabel("ocean turbulence level")
    ax.set_ylabel(r"$\Delta\,/\,V$, both null-subtracted")
    ax.set_title("Exact defect versus its local BKM quadratic approximation")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return True


# ===========================================================================
# 4. Driver
# ===========================================================================

def analyse(files, args):
    rows, meta = [], {}
    for i, f in enumerate(sorted(files, key=level_of)):
        lvl = level_of(f)
        rec, info = blocks_from_clicklog(f, args.coh, args.drop_empty)
        meta[lvl] = info
        if len(rec) < max(8, args.bins * 2):
            print("  L%d: only %d blocks, skipped" % (lvl, len(rec)))
            continue
        rng = np.random.default_rng(SEED + i)
        r = analyse_group(rec, args.bins, args.null_reps, args.bootstrap, rng)
        if r is None:
            print("  L%d: fewer than two populated bins, skipped" % lvl)
            continue
        V0 = r["stat"]["V0"]
        rows.append({
            "level": lvl,
            "source": Path(f).name,
            "n_blocks": info["kept"],
            "n_empty_blocks": info["empty"],
            "n_bins": r["n_bins"],
            "G_pool": r["stat"]["G0"],
            "QBER_pool": r["stat"]["QBER"],
            "V_pool": V0,
            "raw_over_V": r["stat"]["delta"] / V0,
            "null_over_V": r["delta_null"] / V0,
            "excess_over_V": r["excess"] / V0,
            "ci_lo_over_V": r["d_lo"] / V0,
            "ci_hi_over_V": r["d_hi"] / V0,
            "upper_limit95_over_V": r["d_ul95"] / V0,
            "resolved": int(r["d_lo"] > 0),
            "bkm_excess_over_V": r["excess_bkm"] / V0,
            "bkm_ci_lo_over_V": r["b_lo"] / V0,
            "bkm_ci_hi_over_V": r["b_hi"] / V0,
        })
    return pd.DataFrame(rows), meta


def robustness(files, args):
    scan = []
    for nb in (5, 8, 12):
        for drop in (False, True):
            a = argparse.Namespace(**vars(args))
            a.bins = nb
            a.drop_empty = drop
            # A quarter of the replicates: this table is about the direction of
            # the answer, not about its error bar.
            a.null_reps = max(args.null_reps // 4, 100)
            a.bootstrap = max(args.bootstrap // 4, 100)
            tab, _ = analyse(files, a)
            for _, r in tab.iterrows():
                scan.append({"bins": nb, "empty_blocks": "dropped" if drop
                             else "kept", "level": int(r["level"]),
                             "excess_over_V": r["excess_over_V"],
                             "resolved": int(r["resolved"])})
    return pd.DataFrame(scan)


def main():
    ap = argparse.ArgumentParser(
        description="Coarse-graining defect from the existing per-click logs.")
    ap.add_argument("--csv", nargs="*", default=[],
                    help="click logs; default %s" % DEFAULT_GLOB)
    ap.add_argument("--coh", type=int, default=11,
                    help="coherence selector used at collection time; the "
                         "block is 2^(coh+5) attempts (default 11 = 65536)")
    ap.add_argument("--bins", type=int, default=8)
    ap.add_argument("--null-reps", type=int, default=400)
    ap.add_argument("--bootstrap", type=int, default=400)
    ap.add_argument("--distance", type=float, default=25.0,
                    help="range of the section-B sweep, for the model overlay")
    ap.add_argument("--drop-empty", action="store_true",
                    help="discard blocks with no sifted click instead of "
                         "entering them at h = 0")
    ap.add_argument("--robust", action="store_true",
                    help="also scan bins x empty-block handling")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    files = args.csv or sorted(str(p) for p in ROOT.glob(DEFAULT_GLOB))
    if not files:
        raise SystemExit("no click logs matched %s" % DEFAULT_GLOB)

    OUT_DIR.mkdir(exist_ok=True)
    print("coarse-graining defect from %d click log(s), block = %d attempts, "
          "%d bins\n" % (len(files), 1 << (args.coh + 5), args.bins))

    tab, meta = analyse(files, args)
    if tab.empty:
        raise SystemExit("nothing analysable in the supplied logs")

    for lvl in sorted(meta):
        m = meta[lvl]
        print("  L%d  blocks=%d  kept=%d  no-click=%d  multi-h blocks=%d"
              % (lvl, m["n_blocks"], m["kept"], m["empty"], m["split"]))
    print()

    show = tab[["level", "n_blocks", "G_pool", "QBER_pool", "raw_over_V",
                "null_over_V", "excess_over_V", "ci_lo_over_V", "ci_hi_over_V",
                "resolved", "bkm_excess_over_V"]]
    with pd.option_context("display.width", 200,
                           "display.float_format", lambda v: "%.4g" % v):
        print(show.to_string(index=False))

    out_csv = OUT_DIR / "defect_from_clicklog.csv"
    tab.to_csv(out_csv, index=False)
    print("\nwrote %s" % out_csv.relative_to(ROOT))

    model = model_curve(args.distance)
    print("\nanalytical model at %.0f m, same levels:" % args.distance)
    for lvl in sorted(model):
        meas = tab.loc[tab["level"] == lvl, "excess_over_V"]
        inside = ""
        if len(meas):
            r = tab[tab["level"] == lvl].iloc[0]
            if r["ci_lo_over_V"] <= model[lvl] <= r["ci_hi_over_V"]:
                inside = "   <- model inside the measured 95% interval"
        print("  L%d  model %.3e%s" % (lvl, model[lvl], inside))

    if not args.no_figures:
        FIG_DIR.mkdir(exist_ok=True)
        p1 = FIG_DIR / "fig_uwoc_defect_vs_turbulence.png"
        fig_defect_vs_turbulence(tab, model, p1)
        print("\nwrote %s" % p1.relative_to(ROOT))
        p2 = FIG_DIR / "fig_uwoc_exact_vs_bkm.png"
        if fig_exact_vs_bkm(tab, p2):
            print("wrote %s" % p2.relative_to(ROOT))
        else:
            print("exact-vs-BKM figure skipped: fewer than two levels have "
                  "both estimates positive after null subtraction")

    if args.robust:
        print("\nrobustness scan (excess_over_V, * = interval excludes zero):")
        sc = robustness(files, args)
        sc.to_csv(OUT_DIR / "defect_from_clicklog_robustness.csv", index=False)
        piv = sc.pivot_table(index=["bins", "empty_blocks"], columns="level",
                             values="excess_over_V")
        with pd.option_context("display.float_format", lambda v: "%.3e" % v):
            print(piv.to_string())
        print("wrote %s"
              % (OUT_DIR / "defect_from_clicklog_robustness.csv")
              .relative_to(ROOT))

    print("""
LIMITATIONS OF THIS RECONSTRUCTION -- quote them with the numbers.
  1. 160 blocks per level. L1-L3 are expected not to resolve on this sample
     size; that is the planning table, not a null result.
  2. irrad saturates at 255 (h = 2.0), compressing the top bin at L5. The
     12-bit h_mean on the security path does not saturate until h = 16.
  3. Blocks with no sifted click have no observable h and enter at h = 0.
     They are the deep fades that carry the defect. Compare --drop-empty.
  4. trials are not in the per-click log and are set to the nominal N/2. The
     per-attempt estimator is unaffected; the parametric null is slightly
     less accurate.
  There is no Y0 = 0 control here. That run needs the DE1 board and is the
  experimental proof that the defect is the background floor, not counting
  noise.""")


if __name__ == "__main__":
    main()
