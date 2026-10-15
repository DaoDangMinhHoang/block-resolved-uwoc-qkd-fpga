#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_vs_theory.py — Does the FPGA data match the analytical model?
═══════════════════════════════════════════════════════════════════════════════

Reads data/fpga_points.csv (produced by fpga_collect.py) and answers in NUMBERS,
not by feel. Two separate tests:

  [1] CLICK COUNT — validates the emulator's loss + fading stages.
      n_click ~ Poisson(N · P_click_model). Reports z = (measured − expected)/√expected.

  [2] ERROR COUNT — validates the polarisation-error stage e_pol(d) and the Y₀ floor.
      n_err ~ Binomial(n_sift, QBER_model). Compares the Clopper-Pearson interval of the
      measured QBER against the model value: model INSIDE the interval ⇒ match.

⚠ Why not compare "measured QBER = model QBER" directly? At n_sift of a few
  hundred, the Poisson fluctuation of the error count is ±30-50%. A single point
  off by 2× can still be perfectly normal. Only pooling has the resolving power.

Run:  python python/check_vs_theory.py
"""

from __future__ import annotations

import csv
import math
import os
import sys

from scipy.stats import beta as _beta, poisson

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONF = 0.95
MIN_SIFT = 16


def cp(k, n, conf=CONF):
    if n <= 0:
        return 0.0, 1.0
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else float(_beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k >= n else float(_beta.ppf(1 - a, k + 1, n - k))
    return lo, hi


def two_sided_poisson_p(k, mu):
    """Probability of observing a deviation 'at least as large as k' if the model holds."""
    if mu <= 0:
        return 1.0 if k == 0 else 0.0
    p = min(poisson.cdf(k, mu), poisson.sf(k - 1, mu)) * 2
    return min(1.0, p)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        ROOT, "data", "fpga_points.csv")
    if not os.path.exists(path):
        sys.exit("Chưa có %s" % path)

    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            for k in ("p_click", "qber", "p_click_model", "qber_model",
                      "distance_m", "seconds"):
                r[k] = float(r[k])
            for k in ("n_qubit", "n_click", "n_sift", "n_err", "turb", "lam_nm"):
                r[k] = int(r[k])
            rows.append(r)

    print("═" * 100)
    print("  KIỂM ĐỊNH 1 — SỐ CLICK  (suy hao Beer-Lambert + fading)")
    print("═" * 100)
    print("  %-22s %8s %8s %9s %7s %8s" %
          ("cấu hình", "click", "kỳ vọng", "tỉ số", "z", "p"))
    print("  " + "-" * 74)
    tot_c = tot_ce = 0
    for r in rows:
        exp = r["n_qubit"] * r["p_click_model"]
        z = (r["n_click"] - exp) / math.sqrt(max(exp, 1e-9))
        tot_c += r["n_click"]
        tot_ce += exp
        print("  %-22s %8d %8.1f %9.3f %+7.2f %8.3f" %
              (r["tag"].replace("fixed_", "").replace("adaptive_", ""),
               r["n_click"], exp, r["n_click"] / max(exp, 1e-9), z,
               two_sided_poisson_p(r["n_click"], exp)))
    z = (tot_c - tot_ce) / math.sqrt(max(tot_ce, 1e-9))
    print("  " + "-" * 74)
    print("  %-22s %8d %8.1f %9.3f %+7.2f %8.3f" %
          ("GỘP", tot_c, tot_ce, tot_c / max(tot_ce, 1e-9), z,
           two_sided_poisson_p(tot_c, tot_ce)))

    print()
    print("═" * 100)
    print("  KIỂM ĐỊNH 2 — SỐ LỖI  (lỗi phân cực e_pol(d) + nền Y₀)")
    print("═" * 100)
    print("  %-22s %6s %5s %8s %8s %-18s %8s %s" %
          ("cấu hình", "sift", "lỗi", "kỳ vọng", "QBER đo", "khoảng CP 95%",
           "mô hình", "kết luận"))
    print("  " + "-" * 96)
    tot_s = tot_e = tot_ee = 0
    n_out = 0
    for r in rows:
        exp = r["n_sift"] * r["qber_model"]
        lo, hi = cp(r["n_err"], r["n_sift"])
        inside = lo <= r["qber_model"] <= hi
        if r["n_sift"] >= MIN_SIFT:
            tot_s += r["n_sift"]
            tot_e += r["n_err"]
            tot_ee += exp
            if not inside:
                n_out += 1
            verdict = "khớp" if inside else "LỆCH"
        else:
            verdict = "thiếu mẫu"
        print("  %-22s %6d %5d %8.2f %7.2f%% [%6.2f%%,%6.2f%%] %7.2f%% %s" %
              (r["tag"].replace("fixed_", "").replace("adaptive_", ""),
               r["n_sift"], r["n_err"], exp, r["qber"] * 100, lo * 100,
               hi * 100, r["qber_model"] * 100, verdict))

    print("  " + "-" * 96)
    if tot_s:
        lo, hi = cp(tot_e, tot_s)
        q = tot_e / tot_s
        qm = tot_ee / tot_s
        inside = lo <= qm <= hi
        print("  %-22s %6d %5d %8.2f %7.2f%% [%6.2f%%,%6.2f%%] %7.2f%% %s" %
              ("GỘP", tot_s, tot_e, tot_ee, q * 100, lo * 100, hi * 100,
               qm * 100, "khớp" if inside else "LỆCH"))
        print("\n  p (hai phía) cho tổng số lỗi = %.3f"
              % two_sided_poisson_p(tot_e, tot_ee))

    print()
    print("═" * 100)
    print("  TỔNG KẾT")
    print("═" * 100)
    print("  %d cấu hình, %d qubit, %d click, %d bit sàng, %d lỗi, %.2f giờ"
          % (len(rows), sum(r["n_qubit"] for r in rows), tot_c, tot_s, tot_e,
             sum(r["seconds"] for r in rows) / 3600))
    print("  P_click đo / mô hình = %.3f" % (tot_c / max(tot_ce, 1e-9)))
    if tot_s:
        print("  QBER  đo / mô hình = %.3f" % (q / max(qm, 1e-12)))
        print("  %d/%d điểm đủ mẫu có mô hình NẰM NGOÀI khoảng CP 95%%"
              % (n_out, sum(1 for r in rows if r["n_sift"] >= MIN_SIFT)))
        print("  (kỳ vọng ~5% số điểm nằm ngoài ngay cả khi mô hình hoàn toàn đúng)")


if __name__ == "__main__":
    main()
