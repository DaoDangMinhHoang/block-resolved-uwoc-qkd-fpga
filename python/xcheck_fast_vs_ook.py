#!/usr/bin/env python3
"""Phase 2 — compare channel statistics between the OOK path (--fast off) and the fast path.

WHY THIS EXISTS. The security campaign runs on the fast path, which bypasses the
OOK serialisation stage that the original bitstream validated. Nothing in the
data proves the shortcut leaves the channel statistics alone; this script is that
proof, and it is the first thing a referee will ask for.

WHY A BLOCK BOOTSTRAP AND NOT A BINOMIAL TEST. Both legs sample the SAME fading
distribution, and every coherence block is one frozen draw of h. Attempts inside
a block are therefore not independent samples of the channel — they are 65 536
repeats of one draw. Treating them as independent (the naive Wald/Clopper-Pearson
route) understates the error by CV(h)·sqrt(attempts/blocks) and manufactures
disagreements: on the 2026-08-08 vs 2026-09-08 data it reported |z| = 7.7 at L4
and 5.3 at L5 for gains that in fact agree. The resampling unit must be the
BLOCK, in both legs.

The OOK leg is the one that limits the test. It runs at ~3 700 qubit/s against
the fast leg's 10 Mattempt/s, so it buys blocks roughly 2 700x more slowly, and
the ratio's confidence interval is set almost entirely by its block count:

    relative half-width  ~  1.96 * CV(h) / sqrt(n_blocks_ook)

CV(h) climbs from ~0.2 at L1 to ~1.7 at L5, so a --scale 0.05 collection
(8 blocks) can only resolve a factor-2 discrepancy at L5, while the full 160-block
budget reaches ~+/-26%. The script prints that resolving power next to every
verdict: a PASS from an interval that cannot fail is not evidence, and saying so
is more useful than a green tick.

Usage:

    python python\\xcheck_fast_vs_ook.py

    # if the OOK log has not yet been renamed to data/q2/
    python python\\xcheck_fast_vs_ook.py --ook data\\clicks_fixed_B_turb_L%d.csv
"""

import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(ROOT, "data")
POINTS_CSV = os.path.join(DATA_DIR, "fpga_points.csv")

# ATTEMPT_LOG2 in security_block_monitor.v, and coh_sel = 11 gives 2^(11+5) in
# uwoc_channel. The two must agree — that is what makes h_min == h_max hold — so
# one constant covers both legs.
ATT_LOG2 = 16
N_ATT = 1 << ATT_LOG2

LEVELS = (1, 2, 3, 4, 5)

# Below this many resampling units a block bootstrap percentile interval is not
# trustworthy: the resampled sums take only a handful of distinct values, so the
# 2.5/97.5 percentiles land on lattice points rather than on the distribution.
# At --scale 0.05 a level holds eight blocks, which is how L4 produced a "FAIL"
# on a QBER difference the binomial column puts at z = 1.0.
MIN_BLOCKS = 20


# ═══════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════
def ook_budget(level):
    """attempts the host actually ran for section-B level `level`, from the
    checkpoint file. The click log holds only the clicks, so the denominator has
    to come from here — deriving it from the highest attempt index would silently
    shrink whenever a point ends on a dark stretch."""
    tag = "fixed_B_turb_L%d" % level
    if not os.path.exists(POINTS_CSV):
        raise SystemExit("cannot find %s" % POINTS_CSV)
    with open(POINTS_CSV, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["tag"] == tag:
                return int(r["n_qubit"])
    raise SystemExit(
        "no row '%s' in fpga_points.csv — has the OOK leg not finished running?" % tag)


def load_ook(path, n_qubit):
    """→ (clicks_all, sifted, errors) per block, each an array of n_blocks."""
    nb = n_qubit // N_ATT
    if nb < 1:
        raise SystemExit("n_qubit = %d < one block" % n_qubit)
    ca = np.zeros(nb, dtype=np.int64)
    sf = np.zeros(nb, dtype=np.int64)
    er = np.zeros(nb, dtype=np.int64)
    with open(path, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            j = int(r["attempt"]) >> ATT_LOG2
            if j >= nb:                 # tail past the budget, if any
                continue
            ca[j] += 1
            if r["bmatch"] == "1":
                sf[j] += 1
                if r["err"] == "1":
                    er[j] += 1
    return ca, sf, er


def load_fast(path):
    """→ (clicks_all, sifted, errors) per VALID block.

    Invalid blocks are dropped, not repaired: capture_valid = 0 means one of the
    twelve gates failed, and a block that spans two fading realisations or lost
    UART ownership mid-packet is not a measurement of anything.
    """
    ca, sf, er = [], [], []
    nbad = 0
    with open(path, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["capture_valid"] != "1":
                nbad += 1
                continue
            if int(r["attempts"]) != N_ATT:
                nbad += 1
                continue
            ca.append(int(r["clicks_all"]))
            sf.append(int(r["clicks_z"]) + int(r["clicks_x"]))
            er.append(int(r["errors_z"]) + int(r["errors_x"]))
    return (np.array(ca, dtype=np.int64), np.array(sf, dtype=np.int64),
            np.array(er, dtype=np.int64), nbad)


# ═══════════════════════════════════════════════════════════════════════════
# Statistics
# ═══════════════════════════════════════════════════════════════════════════
def pooled(ca, sf, er):
    n = len(ca) * N_ATT
    return (ca.sum() / n,                              # Q_all
            sf.sum() / n,                              # G
            er.sum() / max(sf.sum(), 1))               # QBER


def boot(ca, sf, er, reps, rng):
    """Block bootstrap → (Q_all, G, QBER) each of shape (reps,)."""
    nb = len(ca)
    idx = rng.integers(0, nb, size=(reps, nb))
    n = nb * N_ATT
    s = sf[idx].sum(1)
    return (ca[idx].sum(1) / n, s / n, er[idx].sum(1) / np.maximum(s, 1))


def mark(ok, lo, hi, row, tol):
    """Verdict string. A PASS is annotated, never bare: an interval too wide to
    fail, or one built on too few blocks, is not evidence and must not read like
    it is."""
    if not ok:
        return "FAIL" + ("  (too few blocks)" if row["nb_o"] < MIN_BLOCKS else "")
    if row["nb_o"] < MIN_BLOCKS:
        return "PASS  (too few blocks)"
    if (hi - lo) / 2 > tol * 4:                  # +-20% by default
        return "PASS  (weak)"
    return "PASS"


def cv_h(sf):
    """Coefficient of variation of the per-block sifted count — the empirical
    stand-in for CV(h), and the number that sets how sharp this test can be."""
    m = sf.mean()
    return float(sf.std() / m) if m else float("nan")


# ═══════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser(
        description="Phase 2: compare the OOK path and the fast path")
    ap.add_argument("--ook", default=os.path.join(DATA_DIR, "q2",
                                                  "xcheck_ook_L%d.csv"),
                    help="OOK-leg click-log path pattern, %%d is the turbulence level")
    ap.add_argument("--fast", default=os.path.join(DATA_DIR, "q2",
                                                   "xcheck_fast_L%d.csv"),
                    help="fast-leg block-CSV path pattern")
    ap.add_argument("--reps", type=int, default=4000,
                    help="number of bootstrap repetitions (default 4000)")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--tol", type=float, default=0.05,
                    help="warning threshold: if the confidence interval is wider than +-tol, the test is considered weak")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    rows = []

    for L in LEVELS:
        po, pf = args.ook % L, args.fast % L
        for p in (po, pf):
            if not os.path.exists(p):
                raise SystemExit("missing file: %s" % p)

        oca, osf, oer = load_ook(po, ook_budget(L))
        fca, fsf, fer, nbad = load_fast(pf)
        if len(fsf) == 0:
            raise SystemExit("%s: no valid blocks" % pf)

        bo = boot(oca, osf, oer, args.reps, rng)
        bf = boot(fca, fsf, fer, args.reps, rng)
        rows.append(dict(L=L, tag="L%d" % L,
                         nb_o=len(osf), nb_f=len(fsf), nbad=nbad,
                         cv=cv_h(osf),
                         po=pooled(oca, osf, oer), pf=pooled(fca, fsf, fer),
                         bo=bo, bf=bf,
                         raw_o=(oca, osf, oer), raw_f=(fca, fsf, fer)))

    # ---- pool the 5 levels --------------------------------------------------
    # THIS IS THE ROW THE CLAIM ACTUALLY NEEDS. What phase 2 defends is a DATAPATH
    # property: the fast engine bypasses the OOK serialiser. The datapath does not
    # know which turbulence level it is running, so splitting the evidence five
    # ways buys no information and costs a factor sqrt(5) in precision — worse,
    # eight blocks is below where a block bootstrap is trustworthy at all, which
    # is how a level can come out "FAIL" on a discrepancy a binomial test puts
    # at z ~ 1. Pooling restores the resampling unit count and states the claim
    # at the altitude it is actually made.
    cat = lambda key, i: np.concatenate([r[key][i] for r in rows])
    gca, gsf, ger = cat("raw_o", 0), cat("raw_o", 1), cat("raw_o", 2)
    hca, hsf, her = cat("raw_f", 0), cat("raw_f", 1), cat("raw_f", 2)
    rows.append(dict(L=0, tag="POOLED",
                     nb_o=len(gsf), nb_f=len(hsf),
                     nbad=sum(r["nbad"] for r in rows), cv=cv_h(gsf),
                     po=pooled(gca, gsf, ger), pf=pooled(hca, hsf, her),
                     bo=boot(gca, gsf, ger, args.reps, rng),
                     bf=boot(hca, hsf, her, args.reps, rng),
                     raw_o=(gca, gsf, ger), raw_f=(hca, hsf, her)))

    # ---- table --------------------------------------------------------------
    print()
    print("=" * 78)
    print("  PHASE 2 — OOK leg (--fast off) compared with fast leg (--fast on)")
    print("  BLOCK BOOTSTRAP on both legs, %d repetitions, seed %d"
          % (args.reps, args.seed))
    print("=" * 78)
    print("\n  Sample size: OOK = number of fading blocks; this determines the test's resolving power")
    print("  %-6s %10s %11s %14s %10s %s"
          % ("level", "OOK blocks", "fast blocks", "discarded fast", "CV(h)", ""))
    for r in rows:
        print("  %-6s %10d %11d %14d %10.2f %s"
              % (r["tag"], r["nb_o"], r["nb_f"], r["nbad"], r["cv"],
                 "<-- fewer than 20 blocks: bootstrap NOT yet reliable"
                 if r["nb_o"] < MIN_BLOCKS else ""))

    names = [("Q_all = click/attempt", 0),
             ("G     = sifted/attempt", 1)]
    verdicts = []

    for title, k in names:
        print("\n  %s" % title)
        print("  %-6s %12s %12s %8s %22s %8s"
          % ("level", "OOK", "fast", "ratio", "95% CI of ratio", "verdict"))
        for r in rows:
            a, b = r["po"][k], r["pf"][k]
            ratio = a / b if b else float("nan")
            d = r["bo"][k] / np.maximum(r["bf"][k], 1e-300)
            lo, hi = np.percentile(d, [2.5, 97.5])
            ok = lo <= 1.0 <= hi
            verdicts.append((r["tag"], title.split()[0], ok, r["nb_o"]))
            print("  %-6s %12.4e %12.4e %8.3f   [%6.3f, %6.3f]   %s"
                  % (r["tag"], a, b, ratio, lo, hi, mark(ok, lo, hi, r, args.tol)))

    # QBER carries an independent binomial column ON PURPOSE. The block bootstrap
    # is the right tool for a GAIN, which fading drives directly; QBER is nearly
    # h-independent — that is the paper's central observation — so a binomial
    # comparison is close to correct for it and does not need 20+ blocks to be
    # trusted. When the two columns disagree at low block count, the binomial is
    # the one to believe, and printing both is what makes that visible.
    print("\n  QBER (difference, in percentage points)")
    print("  %-6s %12s %12s %10s %22s %10s %8s"
          % ("level", "OOK", "fast", "difference", "95% bootstrap CI", "binomial z", "verdict"))
    for r in rows:
        a, b = r["po"][2], r["pf"][2]
        d = (r["bo"][2] - r["bf"][2]) * 100
        lo, hi = np.percentile(d, [2.5, 97.5])
        ok = lo <= 0.0 <= hi
        ns = int(r["raw_o"][1].sum())
        nf = int(r["raw_f"][1].sum())
        se = ((a * (1 - a) / max(ns, 1)) + (b * (1 - b) / max(nf, 1))) ** 0.5
        z = (a - b) / se if se else float("nan")
        verdicts.append((r["tag"], "QBER", ok, r["nb_o"]))
        print("  %-6s %11.3f%% %11.3f%% %9.3f   [%6.3f, %6.3f] %10.2f   %s"
              % (r["tag"], a * 100, b * 100, (a - b) * 100, lo, hi, z,
                 mark(ok, lo, hi, r, args.tol)))

    # ---- conclusion ---------------------------------------------------------
    g = [v for v in verdicts if v[0] == "POOLED"]
    npass = sum(1 for v in verdicts if v[2])
    print("\n" + "-" * 78)
    print("  POOLED ROW — the 3/3 tests on the POOLED row are the conclusion of phase 2:")
    print("      %s" % ("All 3 PASS => the fast path DOES NOT change channel statistics."
                        if all(v[2] for v in g) else
                        "A TEST FAILS on the POOLED row — this is the real concern."))
    print("  Individual-level details: %d/%d PASS." % (npass, len(verdicts)))
    bad = [v for v in verdicts if not v[2] and v[0] != "POOLED"]
    if bad:
        print("  Failed at individual levels: %s"
              % ", ".join("%s/%s" % (v[0], v[1]) for v in bad))
        if any(v[3] < MIN_BLOCKS for v in bad):
            print("  But those levels have < %d blocks, so the bootstrap confidence interval is not yet"
                  % MIN_BLOCKS)
            print("  reliable. Compare the 'binomial z' column: |z| < 2 means there is no")
            print("  real discrepancy, only a bootstrap with too few resampling units.")

    print()
    print("  Resolving power of the measurement itself (ratio CI half-width, estimated")
    print("  as 1.96*CV/sqrt(n_blk_OOK)) — the test cannot detect smaller discrepancies:")
    for r in rows:
        hw = 1.96 * r["cv"] / max(r["nb_o"], 1) ** 0.5
        print("    %-6s +-%5.1f %%%s" % (r["tag"], 100 * hw,
              "   <-- too wide; a PASS here says nothing"
              if hw > 0.30 else ""))
    print("-" * 78)
    return 0 if all(v[2] for v in g) else 1


if __name__ == "__main__":
    sys.exit(main())
