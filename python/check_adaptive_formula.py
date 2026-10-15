#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_adaptive_formula.py — Score the ADAPTIVE phase against the right model
═══════════════════════════════════════════════════════════════════════════════

Reads  data/fpga_points.csv  and  data/clicks_<tag>.csv   (produced by fpga_collect.py)
Writes data/adaptive_formula_check.csv   and prints the two rows of README §6.5.

WHY THIS SCRIPT EXISTS

  `p_click_model` in fpga_points.csv is written by model_expect() (fpga_collect.py)
  at the λ the HOST REQUESTED and at the nominal μ-level 8. That is the correct
  reference for the fixed phase and the WRONG one for the adaptive phase, where
  the controller owns both variables and spends most of its time at μ-level 9–12.
  Scoring the adaptive points against that column charges the controller for
  bookkeeping the host never did:

      pooled measured / predicted, at λ_host and μ-level 8   →  1.284
      pooled measured / predicted, at each click's λ and μ   →  0.988

  The 28 % "excess" is entirely λ and μ accounting. Nothing is left over for the
  controller or the RTL to explain — which is the claim §6.5 makes.

HOW THE MODEL IS EVALUATED AT THE LOGGED KNOBS

  λ = LAMBDA_OPTIONS[lam_idx]. μ = LinkConfig.mu · mu_level / 8, mirroring the
  RTL, where p_sig is scaled by `mu_level >> 3` (uwoc_lut_gen.py:423, the
  generator for uwoc_channel.v). mu_level 8 is therefore nominal μ by definition.

THE ONE DEGREE OF FREEDOM — ATTEMPTS OF UNKNOWN STATE

  ⚠ THE PER-CLICK MEAN OF `mu` IS A BIASED ESTIMATOR OF THE APPLIED INTENSITY.
    mu is logged only when a click happened, and a higher μ makes a click likelier,
    so click-weighting over-counts the high-μ stretches. (The same warning is on
    fig_adaptive_controller() in paper_figs_uwoc.py.) Weighting must be by
    ATTEMPT, not by click.

  The controller holds one (mode, μ, λ) per coherence window of `coh_qubits`
  attempts, and this is verified here, not assumed: every window that contains
  clicks is checked for a single state, and a violation is reported as an
  integrity failure rather than averaged away. So a window with at least one
  click has an EXACTLY known state and contributes L_w · p_model(state).

  A window with no clicks does not. Two conventions are computed and both printed,
  because the choice is a judgement call and hiding it would be dishonest:

      covered   — score only the attempts whose state is known, and report the
                  coverage fraction. This is the convention §6.5 quotes.
      hold      — additionally attribute each unobserved window to the last known
                  state (the controller usually does hold), which is the
                  conservative reading: it never shrinks the denominator.

  They differ by ~1 % pooled (0.988 vs 0.979) and the conclusion — "≈ 1, not
  1.284" — does not depend on which one is taken. Where they DO diverge is on the
  far points (harbor 3.5 m: 65 % coverage, 1.54 vs 1.01), because there a window
  is empty for want of photons, not because the controller was probing a dead λ.
  That is why the per-point band is reported over points above --min-coverage
  only; quoting `covered` on a 65 %-coverage point is meaningless.

ERRORS NEED NO RECONSTRUCTION

  Every sifted click carries its own logged λ and μ, so the predicted error count
  is a direct sum of the model QBER over the sifted clicks — no attribution of
  unobserved attempts is involved, and the number is convention-free.

Run:
    python python/check_adaptive_formula.py
    python python/check_adaptive_formula.py --min-coverage 0.99 --verbose
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uwoc_channel_model import (                            # noqa: E402
    LAMBDA_OPTIONS, LinkConfig, UWOCChannel,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
POINTS_CSV = os.path.join(DATA, "fpga_points.csv")
OUT_CSV = os.path.join(DATA, "adaptive_formula_check.csv")

MU_NOMINAL = 8          # mu_level 8 ↔ LinkConfig.mu — the value the fixed phase holds
MIN_CLICKS = 50         # below this a point carries no usable ratio


# ═══════════════════════════════════════════════════════════════════════════
_MODEL_CACHE: dict[tuple, tuple[float, float]] = {}


def model_at(water: str, d: float, turb: int, lam_idx: int, mu_level: int):
    """(p_click, qber) from the channel model at the knobs the CONTROLLER chose.

    Same physics as model_expect() in fpga_collect.py; the difference is that λ
    and μ come from the click log instead of from the host's request.
    """
    key = (water, d, turb, lam_idx, mu_level)
    if key not in _MODEL_CACHE:
        cfg = LinkConfig(lam_nm=LAMBDA_OPTIONS[lam_idx],
                         mu=LinkConfig().mu * mu_level / MU_NOMINAL)
        m = UWOCChannel(water, d, turb, cfg).mean_metrics()
        _MODEL_CACHE[key] = (float(m["p_click"]), float(m["qber"]))
    return _MODEL_CACHE[key]


def load_points():
    if not os.path.exists(POINTS_CSV):
        sys.exit("Chưa có %s — chạy fpga_collect.py trước." % POINTS_CSV)
    rows = []
    with open(POINTS_CSV, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["phase"] != "adaptive":
                continue
            r["distance_m"] = float(r["distance_m"])
            r["p_click_model"] = float(r["p_click_model"])
            r["qber_model"] = float(r["qber_model"])
            for k in ("turb", "n_qubit", "n_click", "n_sift", "n_err", "coh_qubits"):
                r[k] = int(r[k])
            rows.append(r)
    return rows


def load_clicks(tag: str):
    """Per-click telemetry for one tag, or None when the log predates [v13].

    A pre-[v13] bitstream leaves mode/mu/lam_idx as empty strings. That is a
    normal outcome for the fixed_A points and must be skipped, never read as a
    zero — int("") would abort the whole run over a log this script cannot score.
    """
    p = os.path.join(DATA, "clicks_%s.csv" % tag)
    if not os.path.exists(p):
        return None
    att, mu, lam, bmatch, err = [], [], [], [], []
    with open(p, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if not r.get("mode") or not r.get("mu") or not r.get("lam_idx"):
                return None
            att.append(int(r["attempt"]))
            mu.append(int(r["mu"]))
            lam.append(int(r["lam_idx"]))
            bmatch.append(int(r["bmatch"]))
            err.append(int(r["err"]))
    if not att:
        return None
    return dict(attempt=np.array(att), mu=np.array(mu), lam=np.array(lam),
                bmatch=np.array(bmatch), err=np.array(err))


def window_states(c, window: int):
    """{window index: (lam_idx, mu_level)} for every window holding a click.

    Also returns the number of windows carrying MORE THAN ONE state. That count
    must be zero: the controller latches its decision once per coherence window,
    so a mixed window means the log and the RTL disagree about the window grid,
    and every attempt-weighted number below would be built on sand.
    """
    win = c["attempt"] // window
    state, mixed = {}, set()
    for w, li, ml in zip(win.tolist(), c["lam"].tolist(), c["mu"].tolist()):
        if w in state and state[w] != (li, ml):
            mixed.add(w)
        state[w] = (li, ml)
    return state, len(mixed)


def score_point(r):
    """All the ratios for one adaptive point, or None if it cannot be scored."""
    c = load_clicks(r["tag"])
    if c is None or len(c["attempt"]) < MIN_CLICKS:
        return None

    water, d, turb = r["water"], r["distance_m"], r["turb"]
    n_qubit, window = r["n_qubit"], r["coh_qubits"]
    state, mixed = window_states(c, window)

    # ---- clicks: attempt-weighted over windows of known state ---------------
    n_win = (n_qubit + window - 1) // window
    pred_cov = pred_hold = 0.0
    covered = 0
    last = None
    for w in range(n_win):
        length = min(window, n_qubit - w * window)
        if w in state:
            last = state[w]
            p = model_at(water, d, turb, *last)[0]
            pred_cov += length * p
            pred_hold += length * p
            covered += length
        elif last is not None:
            pred_hold += length * model_at(water, d, turb, *last)[0]

    # ---- errors: one term per sifted click, at that click's own knobs -------
    # Gate the log against the point row first. A truncated clicks_<tag>.csv
    # still parses and still produces a plausible ratio — it just silently drops
    # the tail of the run, which is exactly where the controller has settled.
    sift = c["bmatch"] == 1
    log_ok = (len(c["attempt"]) == r["n_click"]
              and int(sift.sum()) == r["n_sift"]
              and int(c["err"].sum()) == r["n_err"])
    pred_err = float(sum(model_at(water, d, turb, int(li), int(ml))[1]
                         for li, ml in zip(c["lam"][sift], c["mu"][sift])))

    # ---- the reference the fixed phase is scored against --------------------
    pred_nom = n_qubit * r["p_click_model"]

    return dict(
        tag=r["tag"], water=water, distance_m=d, turb=turb,
        n_qubit=n_qubit, n_click=r["n_click"], n_sift=r["n_sift"],
        n_err=r["n_err"], coh_qubits=window, n_windows=n_win,
        n_windows_empty=n_win - len(state), n_windows_mixed=mixed,
        log_matches_point=int(log_ok), coverage=covered / n_qubit,
        mu_level_click_wtd=float(c["mu"].mean()),
        pred_clicks_nominal=pred_nom,
        pred_clicks_covered=pred_cov,
        pred_clicks_hold=pred_hold,
        ratio_nominal=r["n_click"] / pred_nom if pred_nom > 0 else float("nan"),
        ratio_covered=r["n_click"] / pred_cov if pred_cov > 0 else float("nan"),
        ratio_hold=r["n_click"] / pred_hold if pred_hold > 0 else float("nan"),
        pred_err=pred_err,
        ratio_err=r["n_err"] / pred_err if pred_err > 0 else float("nan"),
        qber_meas=r["n_err"] / max(r["n_sift"], 1),
        qber_pred=pred_err / max(r["n_sift"], 1),
    )


# ═══════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--min-coverage", type=float, default=0.95,
                    help="points below this attempt coverage are excluded from "
                         "the per-point band (default 0.95)")
    ap.add_argument("--verbose", action="store_true",
                    help="print the per-point table")
    ap.add_argument("--csv", default=OUT_CSV)
    args = ap.parse_args()

    rows = load_points()
    if not rows:
        sys.exit("Không có điểm adaptive nào trong %s." % POINTS_CSV)

    res = [s for s in (score_point(r) for r in rows) if s is not None]
    skipped = len(rows) - len(res)
    if not res:
        sys.exit("Không điểm adaptive nào có telemetry [v13] — không chấm được.")

    bad = [s for s in res if s["n_windows_mixed"]]
    if bad:
        print("⚠ %d điểm có window chứa nhiều hơn một trạng thái (lưới window "
              "của log lệch với RTL) — các số bên dưới KHÔNG dùng được:" % len(bad))
        for s in bad:
            print("    %s: %d window" % (s["tag"], s["n_windows_mixed"]))
    torn = [s for s in res if not s["log_matches_point"]]
    if torn:
        print("⚠ %d điểm có clicks_<tag>.csv không khớp hàng trong fpga_points.csv "
              "(log bị cắt?) — %s"
              % (len(torn), ", ".join(s["tag"] for s in torn)))

    if args.verbose:
        print("%-26s %7s %6s %6s %8s %8s %8s %8s"
              % ("tag", "clicks", "cov%", "mu_cw", "r_nom", "r_cov", "r_hold", "r_err"))
        for s in sorted(res, key=lambda s: (s["water"], s["distance_m"], s["turb"])):
            print("%-26s %7d %6.1f %6.2f %8.3f %8.3f %8.3f %8.3f"
                  % (s["tag"], s["n_click"], 100 * s["coverage"],
                     s["mu_level_click_wtd"], s["ratio_nominal"],
                     s["ratio_covered"], s["ratio_hold"], s["ratio_err"]))
        print()

    cl = sum(s["n_click"] for s in res)
    nom = sum(s["pred_clicks_nominal"] for s in res)
    cov = sum(s["pred_clicks_covered"] for s in res)
    hold = sum(s["pred_clicks_hold"] for s in res)
    er = sum(s["n_err"] for s in res)
    per = sum(s["pred_err"] for s in res)
    sift = sum(s["n_sift"] for s in res)

    print("Scoring of the %d adaptive points%s"
          % (len(res), " (%d skipped: no telemetry)" % skipped if skipped else ""))
    print("  pooled measured clicks                      %10d" % cl)
    print("  at lambda_host, mu-level %d  (p_click_model) %10.0f   ratio %.3f"
          % (MU_NOMINAL, nom, cl / nom))
    print("  at each click's logged lambda and mu")
    print("      covered attempts only                   %10.0f   ratio %.3f"
          % (cov, cl / cov))
    print("      + unobserved windows held at last state %10.0f   ratio %.3f"
          % (hold, cl / hold))

    keep = [s for s in res if s["coverage"] >= args.min_coverage]
    if keep:
        v = np.array([s["ratio_covered"] for s in keep])
        print("  per point (%d of %d with coverage >= %.0f%%): "
              "mean %.3f, sd %.3f, range [%.3f, %.3f]"
              % (len(keep), len(res), 100 * args.min_coverage,
                 v.mean(), v.std(ddof=1), v.min(), v.max()))
    low = [s for s in res if s["coverage"] < args.min_coverage]
    if low:
        print("  excluded from the band (attempt coverage too low): %s"
              % ", ".join("%s %.1f%%" % (s["tag"], 100 * s["coverage"])
                          for s in sorted(low, key=lambda s: s["coverage"])))

    print("  pooled errors  %d measured against %.0f predicted   ratio %.3f"
          % (er, per, er / per))
    print("  pooled QBER    measured %.4f, predicted %.4f  (%d sifted)"
          % (er / max(sift, 1), per / max(sift, 1), sift))

    total_att = sum(s["n_qubit"] for s in res)
    print("  attempt coverage overall: %.2f%% of %s attempts"
          % (100.0 * sum(s["coverage"] * s["n_qubit"] for s in res) / total_att,
             format(total_att, ",")))

    fields = list(res[0].keys())
    with open(args.csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(res, key=lambda s: (s["water"], s["distance_m"], s["turb"])))
    print("\nper-point table written to %s" % args.csv)


if __name__ == "__main__":
    main()
