#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
paper_figs_uwoc.py — Build the paper figures from real FPGA data (UWOC channel)
═══════════════════════════════════════════════════════════════════════════════

Reads  data/fpga_points.csv  and  data/clicks_<tag>.csv   (produced by fpga_collect.py)
Writes Images/fig_uwoc_*.png  +  Images/table_uwoc_results.csv/.tex

PRESENTATION RULES — deliberately different from the atmospheric paper's script:
  [1] EVERY QBER point carries a two-sided Clopper-Pearson interval. Underwater
      P_click ~ 10⁻²…10⁻⁴, so each point holds only a few hundred sifted bits; drawing
      the bare point value would lie about how certain the measurement is.
  [2] Points with n_sift < MIN_SIFT are NOT drawn as measurements but as UPPER-BOUND arrows.
      A point of 0 errors / 3 sifted bits is not "QBER = 0%".
  [3] Any figure holding FPGA data also carries the model curve on the same axes — the
      paper has to show the hardware emulator reproduces the physics correctly, not
      merely that it runs.

Run:
    python python/paper_figs_uwoc.py
"""

from __future__ import annotations

import csv
import glob
import os
import sys

import numpy as np
from scipy.stats import beta as _beta

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                             # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uwoc_channel_model import (                            # noqa: E402
    WATER_TYPES, LinkConfig, UWOCChannel, h2,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "Images")

MIN_SIFT = 16
QBER_LIMIT = 0.11
CONF = 0.95

# ---- IEEE style, matching regen_paper_figs.py so both papers stay consistent ----
FIGW, FIGW2 = 3.5, 7.16
plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.labelsize": 11,
    "legend.fontsize": 8, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "lines.linewidth": 1.5, "lines.markersize": 5,
    "figure.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})

C_FPGA = "#d62728"
C_MODEL = "#1f77b4"
C_LIMIT = "#555555"
C_W = {"clear_ocean": "#1f77b4", "coastal": "#2ca02c", "harbor": "#d62728"}


# ═══════════════════════════════════════════════════════════════════════════
def cp_interval(k: int, n: int, conf: float = CONF):
    """Two-sided Clopper-Pearson interval for k errors out of n trials."""
    if n <= 0:
        return 0.0, 1.0
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else float(_beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k >= n else float(_beta.ppf(1 - a, k + 1, n - k))
    return lo, hi


def load_points():
    p = os.path.join(DATA, "fpga_points.csv")
    if not os.path.exists(p):
        sys.exit("Chưa có %s — chạy fpga_collect.py trước." % p)
    rows = []
    with open(p, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            for k in ("distance_m", "p_click", "sift_eff", "qber", "qber_hi",
                      "skr_per_pulse", "seconds", "rate_qps", "p_click_model",
                      "qber_model"):
                r[k] = float(r[k])
            for k in ("dist_idx", "turb", "lam_nm", "n_qubit", "n_click",
                      "n_sift", "n_err"):
                r[k] = int(r[k])
            r["secure"] = r["secure"] == "True"
            r["qber_lo"], r["qber_hi2"] = cp_interval(r["n_err"], r["n_sift"])
            r["valid"] = r["n_sift"] >= MIN_SIFT
            rows.append(r)
    return rows


def pick(rows, prefix, phase="fixed", turb=None):
    sel = [r for r in rows if r["phase"] == phase
           and r["tag"].startswith("%s_%s" % (phase, prefix))
           and (turb is None or r["turb"] == turb)]
    return sorted(sel, key=lambda r: (r["water"], r["dist_idx"], r["turb"],
                                      r["lam_nm"]))


def turb_levels(rows):
    """Turbulence levels present in the section-A sweep, ascending.

    A range sweep may be repeated at several levels (--turb 1, then --turb 2, …).
    Those are separate series: overlaying them as one would hide the very result
    section B measures — that the MEAN QBER barely moves with turbulence.
    """
    return sorted({r["turb"] for r in pick(rows, "A_dist")})


# Marker per turbulence level, so repeated sweeps stay distinguishable.
MK = {1: "o", 2: "s", 3: "^", 4: "D", 5: "v"}


def model_curve(water, turb, lam_nm, dmax, n=120):
    cfg = LinkConfig(lam_nm=lam_nm)
    dd = np.linspace(0.5, dmax, n)
    pc, qb, sk = [], [], []
    for d in dd:
        m = UWOCChannel(water, float(d), turb, cfg).mean_metrics()
        pc.append(m["p_click"])
        qb.append(m["qber"])
        sk.append(0.5 * m["p_click"] * max(0.0, 1 - 2 * h2(m["qber"])))
    return dd, np.array(pc), np.array(qb), np.array(sk)


def qber_errorbars(rows):
    """Return (x_ok, q_ok, yerr_ok), (x_lim, hi_lim) — measured points and upper-bound-only points."""
    ok = [r for r in rows if r["valid"]]
    lim = [r for r in rows if not r["valid"]]
    xo = np.array([r["distance_m"] for r in ok])
    qo = np.array([r["qber"] for r in ok]) * 100
    lo = qo - np.array([r["qber_lo"] for r in ok]) * 100
    hi = np.array([r["qber_hi2"] for r in ok]) * 100 - qo
    xl = np.array([r["distance_m"] for r in lim])
    hl = np.array([r["qber_hi2"] for r in lim]) * 100
    return (xo, qo, np.vstack([lo, hi])), (xl, hl)


def _save(fig, name):
    fig.tight_layout(pad=0.3)
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("  ✓", name)


# ═══════════════════════════════════════════════════════════════════════════
# FIG 1 — P_click vs range: FPGA against the model  (emulator validation)
# ═══════════════════════════════════════════════════════════════════════════
def fig_pclick(rows):
    A = pick(rows, "A_dist")
    if not A:
        return
    fig, ax = plt.subplots(figsize=(FIGW, 2.7))
    lv = turb_levels(rows)
    dd, pc, _, _ = model_curve("clear_ocean", lv[0], 450,
                               max(r["distance_m"] for r in A) * 1.1)
    ax.semilogy(dd, pc, "-", color=C_MODEL, label="Analytical model")

    for t in lv:
        S = pick(rows, "A_dist", turb=t)
        x = [r["distance_m"] for r in S]
        y = [r["p_click"] for r in S]
        # Poisson error bars on the click count itself
        ye = [r["p_click"] / np.sqrt(max(r["n_click"], 1)) for r in S]
        ax.errorbar(x, y, yerr=ye, fmt=MK.get(t, "o"), color=C_FPGA, capsize=2,
                    mfc="none" if t != lv[0] else None,
                    label="FPGA, L%d" % t if len(lv) > 1 else "FPGA measurement",
                    zorder=3)

    ax.set_xlabel("Link distance, $d$ (m)")
    ax.set_ylabel(r"Detection probability $P_{\rm click}$")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(framealpha=0.0)
    _save(fig, "fig_uwoc_pclick_vs_distance.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 2 — QBER vs range (the headline figure)
# ═══════════════════════════════════════════════════════════════════════════
def fig_qber_distance(rows):
    A = pick(rows, "A_dist")
    if not A:
        return
    fig, ax = plt.subplots(figsize=(FIGW, 2.9))
    dmax = max(r["distance_m"] for r in A) * 1.1
    lv = turb_levels(rows)
    dd, _, qb, _ = model_curve("clear_ocean", lv[0], 450, dmax)
    ax.plot(dd, qb * 100, "-", color=C_MODEL, label="Analytical model")

    top = 60.0
    for t in lv:
        S = pick(rows, "A_dist", turb=t)
        (xo, qo, ye), (xl, hl) = qber_errorbars(S)
        lab = ("FPGA L%d (95%% CP)" % t if len(lv) > 1
               else "FPGA (95% Clopper–Pearson)")
        if len(xo):
            ax.errorbar(xo, qo, yerr=ye, fmt=MK.get(t, "o"), color=C_FPGA,
                        capsize=2, mfc="none" if t != lv[0] else None,
                        label=lab, zorder=3)
            top = max(top, float((qo + ye[1]).max()) * 1.15)
        if len(xl):
            # uplims=True: the convention for "this is an UPPER BOUND"; matplotlib
            # draws a downward arrow. It must not be drawn as a measurement —
            # 0 errors / 3 sifted bits is not a QBER measurement. Only the first
            # level labels it, so repeated sweeps do not repeat the legend entry.
            ax.errorbar(xl, hl, yerr=hl * 0.30, uplims=True, fmt="_",
                        color=C_FPGA, capsize=2, lw=1.0,
                        label=(r"FPGA upper bound ($n_{\rm sift}<%d$)" % MIN_SIFT
                               if t == lv[0] else None),
                        zorder=3)
            top = max(top, float(hl.max()) * 1.15)

    ax.axhline(QBER_LIMIT * 100, color=C_LIMIT, ls="--", lw=1.2,
               label="Shor–Preskill bound (11%)")
    ax.set_xlabel("Link distance, $d$ (m)")
    ax.set_ylabel("QBER (%)")
    ax.set_ylim(0, min(top, 105))
    ax.grid(True, alpha=0.3)
    ax.legend(framealpha=0.0, loc="upper left")
    _save(fig, "fig_uwoc_qber_vs_distance.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 3 — SKR vs range
# ═══════════════════════════════════════════════════════════════════════════
def fig_skr(rows):
    A = pick(rows, "A_dist")
    if not A:
        return
    fig, ax = plt.subplots(figsize=(FIGW, 2.7))
    lv = turb_levels(rows)
    dd, _, _, sk = model_curve("clear_ocean", lv[0], 450,
                               max(r["distance_m"] for r in A) * 1.1)
    sk = np.where(sk > 0, sk, np.nan)
    ax.semilogy(dd, sk, "-", color=C_MODEL, label="Analytical model")

    for t in lv:
        S = pick(rows, "A_dist", turb=t)
        pos = [r for r in S if r["skr_per_pulse"] > 0]
        if pos:
            ax.semilogy([r["distance_m"] for r in pos],
                        [r["skr_per_pulse"] for r in pos], MK.get(t, "o"),
                        color=C_FPGA, mfc="none" if t != lv[0] else None,
                        label="FPGA, L%d" % t if len(lv) > 1
                        else "FPGA measurement")
    for r in [x for x in pick(rows, "A_dist", turb=lv[0])
              if x["skr_per_pulse"] <= 0]:
        ax.annotate("no key", (r["distance_m"], ax.get_ylim()[0] * 3),
                    ha="center", fontsize=7, color=C_LIMIT, rotation=90)
    ax.set_xlabel("Link distance, $d$ (m)")
    ax.set_ylabel("Secure key rate (bits/pulse)")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(framealpha=0.0)
    _save(fig, "fig_uwoc_skr_vs_distance.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 4 — QBER vs turbulence level (expectation: almost flat)
# ═══════════════════════════════════════════════════════════════════════════
def fig_turbulence(rows):
    B = sorted(pick(rows, "B_turb"), key=lambda r: r["turb"])
    if not B:
        return
    fig, ax = plt.subplots(figsize=(FIGW, 2.7))
    x = np.array([r["turb"] for r in B])

    qm = []
    for r in B:
        qm.append(UWOCChannel("clear_ocean", r["distance_m"], r["turb"],
                              LinkConfig(lam_nm=r["lam_nm"])
                              ).mean_metrics()["qber"] * 100)
    ax.plot(x, qm, "-s", color=C_MODEL, mfc="none", label="Analytical model")

    ok = [r for r in B if r["valid"]]
    if ok:
        xo = np.array([r["turb"] for r in ok])
        qo = np.array([r["qber"] for r in ok]) * 100
        lo = qo - np.array([r["qber_lo"] for r in ok]) * 100
        hi = np.array([r["qber_hi2"] for r in ok]) * 100 - qo
        ax.errorbar(xo, qo, yerr=np.vstack([lo, hi]), fmt="o", color=C_FPGA,
                    capsize=2, label="FPGA (95% CP)", zorder=3)

    ax.set_xlabel("Oceanic turbulence level")
    ax.set_ylabel("QBER (%)")
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xticklabels(["L1\nvery weak", "L2\nweak", "L3\nmod.",
                        "L4\nstrong", "L5\nsevere"], fontsize=7)
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(framealpha=0.0)
    _save(fig, "fig_uwoc_qber_vs_turbulence.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 5 — The optimal λ depends on the water type (lam ↔ red)
# ═══════════════════════════════════════════════════════════════════════════
def fig_wavelength(rows):
    C = pick(rows, "C_lam")
    if not C:
        return
    groups = {}
    for r in C:
        groups.setdefault(r["water"], []).append(r)
    if not groups:
        return

    fig, axes = plt.subplots(1, len(groups), figsize=(FIGW2, 2.6), sharey=False)
    axes = np.atleast_1d(axes)
    for ax, (water, rs) in zip(axes, sorted(groups.items())):
        rs = sorted(rs, key=lambda r: r["lam_nm"])
        lams = [r["lam_nm"] for r in rs]
        xi = np.arange(len(lams))
        meas = [r["p_click"] for r in rs]
        mod = [r["p_click_model"] for r in rs]
        err = [r["p_click"] / np.sqrt(max(r["n_click"], 1)) for r in rs]

        ax.bar(xi - 0.2, meas, 0.4, yerr=err, capsize=2, color=C_FPGA,
               edgecolor="k", lw=0.4, label="FPGA")
        ax.bar(xi + 0.2, mod, 0.4, color=C_MODEL, edgecolor="k", lw=0.4,
               label="Model")
        ax.set_xticks(xi)
        ax.set_xticklabels(["%d nm" % l for l in lams])
        ax.set_title("%s, $d$ = %g m" % (WATER_TYPES[water].name,
                                         rs[0]["distance_m"]), fontsize=9)
        ax.set_ylabel(r"$P_{\rm click}$")
        ax.grid(True, alpha=0.3, axis="y")
        ax.legend(framealpha=0.0, fontsize=7)
    _save(fig, "fig_uwoc_wavelength.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 6 — The three water types
# ═══════════════════════════════════════════════════════════════════════════
def fig_water_types(rows):
    pts = pick(rows, "A_dist") + pick(rows, "D_")
    if not pts:
        return
    groups = {}
    for r in pts:
        groups.setdefault(r["water"], []).append(r)

    fig, ax = plt.subplots(figsize=(FIGW, 2.7))
    for water, rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: r["distance_m"])
        c = C_W.get(water, "k")
        dmax = max(r["distance_m"] for r in rs) * 1.15
        dd, pc, _, _ = model_curve(water, rs[0]["turb"], rs[0]["lam_nm"], dmax)
        ax.semilogy(dd, pc, "-", color=c, lw=1.2, alpha=0.7)
        ax.semilogy([r["distance_m"] for r in rs], [r["p_click"] for r in rs],
                    "o", color=c, label=WATER_TYPES[water].name)
    ax.set_xlabel("Link distance, $d$ (m)")
    ax.set_ylabel(r"$P_{\rm click}$")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(framealpha=0.0, title="Markers: FPGA · Lines: model",
              title_fontsize=7)
    _save(fig, "fig_uwoc_water_types.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 7 — METHODOLOGY figure: why N = 5,000–10,000 qubits is not enough
# ═══════════════════════════════════════════════════════════════════════════
def fig_sample_size(rows):
    """Achievable QBER upper bound versus transmitted qubits, for several ranges.

    This figure answers the question 'why did the measured QBER come out 0%'
    directly. With 0 errors over n sifted bits the 95% upper bound is 1−0.05^(1/n)
    ≈ 3/n. Asserting QBER < 11% needs n_sift ≳ 27; MEASURING QBER ~1.6% needs n_sift ~10³.
    """
    fig, ax = plt.subplots(figsize=(FIGW, 2.7))
    N = np.logspace(3, 6.3, 200)
    for d, p in ((5.0, 1.039e-2), (15.0, 3.511e-3), (25.0, 1.18e-3)):
        nsift = N * p * 0.5
        hi = np.array([cp_interval(0, max(int(n), 1))[1] for n in nsift]) * 100
        ax.loglog(N, hi, label="$d$ = %g m" % d)

    ax.axhline(QBER_LIMIT * 100, color=C_LIMIT, ls="--", lw=1.2,
               label="11% bound")
    ax.axvspan(5e3, 1e4, color="#999", alpha=0.25)
    ax.annotate("$N$ of the\natmospheric study", xy=(7e3, 0.35),
                xytext=(2.2e4, 0.16), fontsize=6.5, color="#444",
                arrowprops=dict(arrowstyle="->", color="#444", lw=0.7))
    ax.set_xlabel("Transmitted pulses, $N$")
    ax.set_ylabel("QBER upper bound (%)")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(framealpha=0.0, fontsize=7)
    _save(fig, "fig_uwoc_sample_size.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 8 — QBER versus the SNR proxy (irrad), from the per-click data
# ═══════════════════════════════════════════════════════════════════════════
def fig_qber_vs_irradiance():
    recs = []
    for p in glob.glob(os.path.join(DATA, "clicks_*.csv")):
        with open(p, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if int(r["bmatch"]):
                    recs.append((int(r["irrad"]), int(r["err"])))
    if len(recs) < 200:
        print("  (bỏ qua fig_uwoc_qber_vs_snr: mới %d bit sàng)" % len(recs))
        return
    irr = np.array([a for a, _ in recs])
    err = np.array([b for _, b in recs])

    edges = np.quantile(irr, np.linspace(0, 1, 7))
    edges = np.unique(edges)
    xs, qs, los, his = [], [], [], []
    for i in range(len(edges) - 1):
        # ⚠ the parentheses are mandatory: Python reads `a & b if c else d` as
        #   `(a & b) if c else d`, which loses the lower bound in the else branch.
        upper = ((irr <= edges[i + 1]) if i == len(edges) - 2
                 else (irr < edges[i + 1]))
        m = (irr >= edges[i]) & upper
        n = int(m.sum())
        if n < 20:
            continue
        k = int(err[m].sum())
        lo, hi = cp_interval(k, n)
        xs.append(float(irr[m].mean()))
        qs.append(k / n * 100)
        los.append(qs[-1] - lo * 100)
        his.append(hi * 100 - qs[-1])
    if not xs:
        return

    fig, ax = plt.subplots(figsize=(FIGW, 2.7))
    ax.errorbar(xs, qs, yerr=np.vstack([los, his]), fmt="o-", color=C_FPGA,
                capsize=2, label="FPGA, all configurations pooled")
    ax.axhline(QBER_LIMIT * 100, color=C_LIMIT, ls="--", lw=1.2,
               label="11% bound")
    ax.set_xlabel("Received irradiance proxy (8-bit)")
    ax.set_ylabel("QBER (%)")
    ax.grid(True, alpha=0.3)
    ax.legend(framealpha=0.0, fontsize=7)
    _save(fig, "fig_uwoc_qber_vs_snr.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 8 — fixed vs adaptive, the only figure that uses the adaptive phase
# ═══════════════════════════════════════════════════════════════════════════
def fig_adaptive(rows):
    """QBER and SKR vs range for both phases on shared axes.

    In PC-input mode the PC supplies the bases, so adapt_basis_prob is bypassed
    (top_module.v:280-283). What the controller still changes is active_lambda,
    active_power, eff_gap and tx_permitted — so the honest claim this figure can
    support is "the controller picks λ/power per channel state", not "it biases
    the basis". Note that tx_permitted = 0 (controller PAUSE) produces a point
    with no clicks at all: that is a result, not a failed measurement.
    """
    F, A = pick(rows, "A_dist", "fixed"), pick(rows, "A_dist", "adaptive")
    if not F or not A:
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIGW2, 2.8))

    for rr, lab, col, mk in ((F, "Fixed", C_MODEL, "o"),
                             (A, "Adaptive", C_FPGA, "s")):
        (xo, qo, ye), (xl, hl) = qber_errorbars(rr)
        if len(xo):
            ax1.errorbar(xo, qo, yerr=ye, fmt=mk + "-", color=col, capsize=2,
                         mfc="none", label=lab, zorder=3)
        if len(xl):
            ax1.errorbar(xl, hl, yerr=hl * 0.30, uplims=True, fmt="_",
                         color=col, capsize=2, lw=1.0, zorder=3)
        pos = [r for r in rr if r["skr_per_pulse"] > 0]
        if pos:
            ax2.semilogy([r["distance_m"] for r in pos],
                         [r["skr_per_pulse"] for r in pos], mk + "-",
                         color=col, mfc="none", label=lab)

    ax1.axhline(QBER_LIMIT * 100, color=C_LIMIT, ls="--", lw=1.2,
                label="11% bound")
    ax1.set_xlabel("Link distance, $d$ (m)")
    ax1.set_ylabel("QBER (%)")
    ax2.set_xlabel("Link distance, $d$ (m)")
    ax2.set_ylabel("Secure key rate (bits/pulse)")
    for ax in (ax1, ax2):
        ax.grid(True, alpha=0.3, which="both")
        ax.legend(framealpha=0.0)
    _save(fig, "fig_uwoc_fixed_vs_adaptive.png")


# ═══════════════════════════════════════════════════════════════════════════
# CONTROLLER TELEMETRY — the [v13]/[v14] per-click mode / mu / lam_idx columns
# ═══════════════════════════════════════════════════════════════════════════
# These three columns were added to the report line in [v13], and res_lam was
# only actually assigned in [v14] (top_module.v:737). Logs taken before that
# carry EMPTY strings — every fixed_A point except d8 is in that state — so
# every reader below must treat "no telemetry" as a normal outcome and skip,
# never as a zero. A silent int("") would abort the whole figure run.

MODE_NAMES = ("AGGRESSIVE", "MODERATE", "CONSERVATIVE", "PAUSE")
# Sequential, dark = throttled. Deliberately NOT the categorical palette: the
# four modes are an ordered escalation, and a reader must see that at a glance.
MODE_COLORS = ("#4a7fb5", "#8bb2d9", "#f4a582", "#b2182b")
LAM_NAMES = ("450 nm", "532 nm", "650 nm")
LAM_COLORS = ("#2166ac", "#4daf4a", "#d62728")
MU_NOMINAL = 8          # mu_level 8 ↔ LinkConfig.mu, the value the fixed phase holds


def load_clicks(tag):
    """Per-click arrays for one tag, or None when the file has no telemetry.

    Returns a dict of numpy arrays. `mode`, `mu` and `lam` are present only when
    the log actually carries them; a pre-[v13] log yields None so the caller can
    skip the point rather than plot a fabricated zero.
    """
    p = os.path.join(DATA, "clicks_%s.csv" % tag)
    if not os.path.exists(p):
        return None
    att, irr, mode, mu, lam, bm, er = [], [], [], [], [], [], []
    with open(p, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            # A blank telemetry field means the bitstream predates the column.
            if not r.get("mode") or not r.get("mu"):
                return None
            att.append(int(r["attempt"]) if r.get("attempt") else 0)
            irr.append(int(r["irrad"]))
            mode.append(int(r["mode"]))
            mu.append(int(r["mu"]))
            lam.append(int(r["lam_idx"]) if r.get("lam_idx") else 0)
            bm.append(int(r["bmatch"]))
            er.append(int(r["err"]))
    if not mode:
        return None
    return dict(attempt=np.array(att), irrad=np.array(irr),
                mode=np.array(mode), mu=np.array(mu), lam=np.array(lam),
                bmatch=np.array(bm), err=np.array(er))


def _fractions(v, k):
    """Fraction of samples in each of k categories."""
    n = max(len(v), 1)
    return [float((v == i).sum()) / n for i in range(k)]


# ═══════════════════════════════════════════════════════════════════════════
# FIG 9 — what the adaptive controller actually DID, versus range
# ═══════════════════════════════════════════════════════════════════════════
def fig_adaptive_controller(rows):
    """Mode / intensity / wavelength chosen by the controller, versus range.

    This is the figure that shows the controller WORKS, as opposed to merely
    runs — and it is stronger evidence than the fixed-vs-adaptive QBER panel,
    because the two phases' QBER differ by little while these three curves move
    monotonically and for physically explainable reasons.

    ⚠ PAUSE IN PANEL (a) IS A DECISION, NOT AN ACTION. `tx_permitted` is read
    only in the non-PC branch of FSM_IDLE (top_module.v:674); under SW[9] = 1,
    which is the mode fpga_collect.py drives, the FSM goes straight to
    FSM_WAIT_CMD and transmits regardless. So the honest reading is "the
    controller identified these windows as unsafe", not "the link stopped".
    Anything stronger needs tx_permitted gating the PC branch too.

    ⚠ Panel (b) is a CLICK-WEIGHTED mean intensity, and that is a biased
    estimator of the intensity the controller applied: mu is only logged when a
    click happened, and a higher mu makes a click likelier. The true
    attempt-weighted mean is lower. The axis label says so; do not quote this
    number as "the intensity used" in the text.
    """
    A = pick(rows, "A_dist", "adaptive")
    if not A:
        print("  (bỏ qua fig_uwoc_adaptive_controller: chưa có điểm adaptive)")
        return
    d, mode_f, mu_m, lam_f = [], [], [], []
    for r in A:
        c = load_clicks(r["tag"])
        if c is None or len(c["mode"]) < 50:
            continue
        d.append(r["distance_m"])
        mode_f.append(_fractions(c["mode"], 4))
        lam_f.append(_fractions(c["lam"], 3))
        mu_m.append(float(c["mu"].mean()))
    if len(d) < 3:
        print("  (bỏ qua fig_uwoc_adaptive_controller: %d điểm có telemetry)" % len(d))
        return

    d = np.array(d)
    mode_f = np.array(mode_f).T * 100          # (4, n)
    lam_f = np.array(lam_f).T * 100            # (3, n)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(FIGW2, 2.5))

    ax1.stackplot(d, mode_f, labels=MODE_NAMES, colors=MODE_COLORS, alpha=0.9)
    ax1.set_ylim(0, 100)
    ax1.set_ylabel("Time in mode (%)")
    ax1.legend(framealpha=0.0, fontsize=6, loc="center left")

    ax2.plot(d, mu_m, "o-", color=C_FPGA, mfc="none", label="Adaptive")
    ax2.axhline(MU_NOMINAL, color=C_LIMIT, ls="--", lw=1.2,
                label="Fixed ($\\mu_{\\rm nom}$)")
    ax2.set_ylabel("$\\mu$ level (click-wtd)")
    ax2.set_ylim(bottom=min(MU_NOMINAL - 0.5, mu_m[0] - 0.5))
    ax2.legend(framealpha=0.0, fontsize=7)

    ax3.stackplot(d, lam_f, labels=LAM_NAMES, colors=LAM_COLORS, alpha=0.9)
    ax3.set_ylim(0, 100)
    ax3.set_ylabel("Time on $\\lambda$ (%)")
    ax3.legend(framealpha=0.0, fontsize=6, loc="center left")

    for ax in (ax1, ax2, ax3):
        ax.set_xlabel("Link distance, $d$ (m)")
        ax.grid(True, alpha=0.3)
    _save(fig, "fig_uwoc_adaptive_controller.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 10 — the controller reacting to fading WITHIN one measurement point
# ═══════════════════════════════════════════════════════════════════════════
def fig_adaptive_timeline(rows):
    """Irradiance, mode and intensity against the attempt index of one point.

    The `attempt` column ([v12], 24-bit, reset by the host's 0x01) is what makes
    a real time axis possible: total_cnt counts CLICKS, not attempts, because a
    lost photon never reaches FSM_PROCESS.

    Point selection is a fixed rule, not a hand-pick: take the adaptive point
    whose PAUSE fraction is closest to 0.3, so both the transmitting and the
    throttled regime occupy a visible share of the axis. A point at 94% PAUSE
    (45 m) would be a flat band and show nothing.

    ⚠ Clicks continue through the PAUSE bands, and that is not a bug in the
    plot: in PC-input mode tx_permitted never gates the FSM (see
    fig_adaptive_controller). The band marks the controller's verdict on the
    channel, and the panel above shows the fading dip that produced it.
    """
    A = pick(rows, "A_dist", "adaptive")
    best, best_score = None, None
    for r in A:
        c = load_clicks(r["tag"])
        if c is None or len(c["mode"]) < 300:
            continue
        pause = float((c["mode"] == 3).mean())
        score = abs(pause - 0.30)
        if best_score is None or score < best_score:
            best, best_score, best_c = r, score, c
    if best is None:
        print("  (bỏ qua fig_uwoc_adaptive_timeline: không có điểm đủ click)")
        return
    c = best_c

    # Attempts are the x axis, in millions of pulses.
    x = c["attempt"] / 1e6
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(FIGW2, 3.2), sharex=True)

    # Raw irrad is one sample per click and very noisy; the rolling mean is what
    # the controller effectively sees through channel_monitor's window.
    w = max(9, len(x) // 120)
    ker = np.ones(w) / w
    irr_s = np.convolve(c["irrad"].astype(float), ker, mode="same")
    ax1.plot(x, c["irrad"], ".", ms=1.2, color="#bbbbbb", label="per click")
    ax1.plot(x, irr_s, "-", color=C_MODEL, lw=1.2, label="rolling mean")
    ax1.set_ylabel("Irradiance proxy")
    ax1.legend(framealpha=0.0, fontsize=7, ncol=2)

    # Mode as a filled band: one vertical stripe per click, coloured by mode.
    for k in range(4):
        m = c["mode"] == k
        if m.any():
            ax2.plot(x[m], np.full(m.sum(), k), "|", ms=8,
                     color=MODE_COLORS[k], label=MODE_NAMES[k])
    ax2.set_yticks(range(4))
    ax2.set_yticklabels(MODE_NAMES, fontsize=7)
    ax2.set_ylim(-0.5, 3.5)
    ax2.set_xlabel("Attempt index ($10^6$ pulses)")
    ax2.invert_yaxis()

    for ax in (ax1, ax2):
        ax.grid(True, alpha=0.3)
    ax1.set_title("%s — $d$ = %g m, L%d"
                  % (best["tag"], best["distance_m"], best["turb"]),
                  fontsize=8)
    _save(fig, "fig_uwoc_adaptive_timeline.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 11 — fixed → adaptive gain, WITH the intensity confound made explicit
# ═══════════════════════════════════════════════════════════════════════════
def fig_adaptive_gain(rows):
    """P_click gain, security margin and time cost, fixed versus adaptive.

    ⚠ Panel (a) deliberately overlays the INTENSITY ratio on the P_click ratio.
    Since p_sig scales as mu_level/8 in uwoc_channel.v, most of the throughput
    gain is bought by raising mu, not by cleverness — and skr() is the
    decoy-free Shor-Preskill rate, which charges nothing for the extra
    multi-photon fraction that Eve can split. Plotting the two ratios together
    is what keeps the claim honest; a reviewer will compute it anyway.
    """
    F = {r["dist_idx"]: r for r in pick(rows, "A_dist", "fixed")}
    A = {r["dist_idx"]: r for r in pick(rows, "A_dist", "adaptive")}
    common = sorted(set(F) & set(A))
    if len(common) < 3:
        print("  (bỏ qua fig_uwoc_adaptive_gain: chỉ %d cự ly có cả hai pha)"
              % len(common))
        return

    d = np.array([F[i]["distance_m"] for i in common])
    pc_ratio = np.array([A[i]["p_click"] / max(F[i]["p_click"], 1e-30)
                         for i in common])
    mu_ratio, mu_d = [], []
    for i in common:
        ca, cf = load_clicks(A[i]["tag"]), load_clicks(F[i]["tag"])
        if ca is None:
            continue
        # The fixed phase holds mu_level at nominal; use its logged value when
        # the log has one, else the nominal constant.
        mf = float(cf["mu"].mean()) if cf is not None else float(MU_NOMINAL)
        mu_ratio.append(float(ca["mu"].mean()) / mf)
        mu_d.append(F[i]["distance_m"])

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(FIGW2, 2.5))

    ax1.plot(d, pc_ratio, "o-", color=C_FPGA, mfc="none",
             label="$P_{\\rm click}$ ratio")
    if mu_ratio:
        ax1.plot(mu_d, mu_ratio, "s--", color=C_MODEL, mfc="none",
                 label="$\\mu$ ratio (confound)")
    ax1.axhline(1.0, color=C_LIMIT, ls=":", lw=1.0)
    ax1.set_ylabel("Adaptive / fixed")
    ax1.legend(framealpha=0.0, fontsize=7)

    for src, lab, col, mk in ((F, "Fixed", C_MODEL, "o"),
                              (A, "Adaptive", C_FPGA, "s")):
        ax2.plot([src[i]["distance_m"] for i in common],
                 [src[i]["qber_hi2"] * 100 for i in common],
                 mk + "-", color=col, mfc="none", label=lab)
    ax2.axhline(QBER_LIMIT * 100, color=C_LIMIT, ls="--", lw=1.2,
                label="11% bound")
    ax2.set_ylabel("QBER upper bound (%)")
    ax2.legend(framealpha=0.0, fontsize=7)

    # Wall-clock per sifted bit — the cost side of the same trade.
    for src, lab, col, mk in ((F, "Fixed", C_MODEL, "o"),
                              (A, "Adaptive", C_FPGA, "s")):
        ax3.semilogy([src[i]["distance_m"] for i in common],
                     [src[i]["seconds"] / max(src[i]["n_sift"], 1)
                      for i in common],
                     mk + "-", color=col, mfc="none", label=lab)
    ax3.set_ylabel("Time per sifted bit (s)")
    ax3.legend(framealpha=0.0, fontsize=7)

    for ax in (ax1, ax2, ax3):
        ax.set_xlabel("Link distance, $d$ (m)")
        ax.grid(True, alpha=0.3, which="both")
    _save(fig, "fig_uwoc_adaptive_gain.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 12 — per-block QBER DISPERSION against the model (section B)
# ═══════════════════════════════════════════════════════════════════════════
def fig_block_dispersion():
    """Measured vs simulated spread of the per-block QBER, per turbulence level.

    Every other figure in this script validates a MEAN. This one validates a
    VARIANCE, which is the quantity the coherence-block fix of [v12] exists to
    make observable at all — with h resampled per qubit the measured spread
    collapses by √65536 and L1 becomes indistinguishable from L5.

    Reads data/block_table.csv (written by sim_table.py --blocks); that file only
    has rows for runs taken with --coh >= 1.
    """
    p = os.path.join(DATA, "block_table.csv")
    if not os.path.exists(p):
        print("  (bỏ qua fig_uwoc_block_dispersion: chưa có block_table.csv)")
        return
    # FIXED ONLY. "_B_turb_" alone would also match adaptive_B_turb_L*, and the
    # two phases would then be sorted together into one series with five
    # duplicated turbulence levels. This figure validates the emulator's variance
    # against the model, which is a fixed-phase question; the adaptive section-B
    # run has its own figure (fig_adaptive_turbulence).
    B = []
    with open(p, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["tag"].startswith("fixed_B_turb_"):
                B.append(r)
    if len(B) < 3:
        print("  (bỏ qua fig_uwoc_block_dispersion: %d điểm mục B)" % len(B))
        return
    B.sort(key=lambda r: int(r["turb"]))

    lv = np.array([int(r["turb"]) for r in B])
    std_m = np.array([float(r["qber_win_std_pct"]) for r in B])
    std_s = np.array([float(r["sim_qber_win_std_pct"]) for r in B])
    out_m = np.array([float(r["outage"]) for r in B]) * 100
    out_s = np.array([float(r["sim_outage"]) for r in B]) * 100
    mean_m = np.array([float(r["qber_win_mean_pct"]) for r in B])
    mean_s = np.array([float(r["sim_qber_win_mean_pct"]) for r in B])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIGW2, 2.6))
    w = 0.35

    # The per-block QBER is strongly right-skewed at high turbulence (Weibull h),
    # so a symmetric ±std bar would dip below zero — L5 has std 24.5% against a
    # mean of 14%. A negative QBER is not a thing; clip the lower whisker at 0 and
    # let the asymmetry itself carry the message that the spread is one-sided.
    def _clip_err(mean, std):
        return np.vstack([np.minimum(std, mean), std])

    ax1.bar(lv - w / 2, mean_m, w, color=C_FPGA, alpha=0.85, label="FPGA mean")
    ax1.bar(lv + w / 2, mean_s, w, color=C_MODEL, alpha=0.85, label="Model mean")
    ax1.errorbar(lv - w / 2, mean_m, yerr=_clip_err(mean_m, std_m), fmt="none",
                 ecolor="#333333", capsize=3, lw=1.0,
                 label="$\\pm$ per-block std")
    ax1.errorbar(lv + w / 2, mean_s, yerr=_clip_err(mean_s, std_s), fmt="none",
                 ecolor="#333333", capsize=3, lw=1.0)
    ax1.set_ylim(bottom=0)
    ax1.set_ylabel("Per-block QBER (%)")
    ax1.legend(framealpha=0.0, fontsize=7)

    ax2.plot(lv, std_m, "o-", color=C_FPGA, mfc="none", label="FPGA std")
    ax2.plot(lv, std_s, "s--", color=C_MODEL, mfc="none", label="Model std")
    ax2b = ax2.twinx()
    ax2b.plot(lv, out_m, "^:", color="#7b3294", mfc="none", label="FPGA outage")
    ax2b.plot(lv, out_s, "v:", color="#c2a5cf", mfc="none", label="Model outage")
    ax2b.set_ylabel("Outage (%)", fontsize=9)
    ax2.set_ylabel("Std of block QBER (%)")
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2b.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, framealpha=0.0, fontsize=6)

    for ax in (ax1, ax2):
        ax.set_xlabel("Ocean turbulence level")
        ax.set_xticks(lv)
        ax.set_xticklabels(["L%d" % v for v in lv])
        ax.grid(True, alpha=0.3)
    _save(fig, "fig_uwoc_block_dispersion.png")


# ═══════════════════════════════════════════════════════════════════════════
# FIG 13 — the controller against TURBULENCE at fixed range (section B)
# ═══════════════════════════════════════════════════════════════════════════
def fig_adaptive_turbulence(rows):
    """Controller response versus turbulence level, range held at 25 m.

    The counterpart to fig_adaptive_controller, and the more fundamental of the
    two. A range sweep degrades the link DETERMINISTICALLY, so a controller that
    merely tracked mean received power would produce the same curves. Section B
    holds the range fixed and moves only sigma^2_ho, i.e. only the VARIANCE —
    the quantity the controller exists to track. Mode and mu rising with L here,
    at constant mean P_click, is the result that a range sweep cannot give.

    Draws nothing until `--phase adaptive --sections B` has been collected.
    """
    A = pick(rows, "B_turb", "adaptive")
    if not A:
        print("  (bỏ qua fig_uwoc_adaptive_turbulence: chưa đo adaptive mục B)")
        return
    F = {r["turb"]: r for r in pick(rows, "B_turb", "fixed")}

    lv, mode_f, mu_m, pc_a, pc_f = [], [], [], [], []
    for r in sorted(A, key=lambda r: r["turb"]):
        c = load_clicks(r["tag"])
        if c is None or len(c["mode"]) < 50:
            continue
        lv.append(r["turb"])
        mode_f.append(_fractions(c["mode"], 4))
        mu_m.append(float(c["mu"].mean()))
        pc_a.append(r["p_click"])
        pc_f.append(F[r["turb"]]["p_click"] if r["turb"] in F else np.nan)
    if len(lv) < 2:
        print("  (bỏ qua fig_uwoc_adaptive_turbulence: %d mức có telemetry)" % len(lv))
        return

    lv = np.array(lv)
    mode_f = np.array(mode_f).T * 100

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(FIGW2, 2.5))

    ax1.stackplot(lv, mode_f, labels=MODE_NAMES, colors=MODE_COLORS, alpha=0.9)
    ax1.set_ylim(0, 100)
    ax1.set_ylabel("Time in mode (%)")
    ax1.legend(framealpha=0.0, fontsize=6, loc="center left")

    ax2.plot(lv, mu_m, "o-", color=C_FPGA, mfc="none", label="Adaptive")
    ax2.axhline(MU_NOMINAL, color=C_LIMIT, ls="--", lw=1.2,
                label="Fixed ($\\mu_{\\rm nom}$)")
    ax2.set_ylabel("$\\mu$ level (click-wtd)")
    ax2.legend(framealpha=0.0, fontsize=7)

    # The control panel: the MEAN P_click is flat in L (P_click is linear in h,
    # so the expectation cancels the fading). If this panel slopes, the point is
    # not measuring what section B claims to measure.
    ax3.plot(lv, pc_a, "s-", color=C_FPGA, mfc="none", label="Adaptive")
    ax3.plot(lv, pc_f, "o--", color=C_MODEL, mfc="none", label="Fixed")
    ax3.set_ylabel("$P_{\\rm click}$ (mean)")
    ax3.legend(framealpha=0.0, fontsize=7)

    for ax in (ax1, ax2, ax3):
        ax.set_xlabel("Ocean turbulence level")
        ax.set_xticks(lv)
        ax.set_xticklabels(["L%d" % v for v in lv])
        ax.grid(True, alpha=0.3)
    _save(fig, "fig_uwoc_adaptive_turbulence.png")


# ═══════════════════════════════════════════════════════════════════════════
def write_table(rows):
    cols = [("water", "Water"), ("distance_m", "d (m)"), ("turb", "L"),
            ("lam_nm", "lambda (nm)"), ("n_qubit", "N"), ("n_click", "clicks"),
            ("n_sift", "n_sift"), ("n_err", "errors"),
            ("p_click", "P_click"), ("p_click_model", "P_click (model)"),
            ("qber", "QBER"), ("qber_hi2", "QBER upper 95%"),
            ("qber_model", "QBER (model)"), ("skr_per_pulse", "SKR/pulse")]
    p = os.path.join(OUT, "table_uwoc_results.csv")
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([h for _, h in cols])
        for r in rows:
            w.writerow([r[k] for k, _ in cols])
    print("  ✓ table_uwoc_results.csv")

    p = os.path.join(OUT, "table_uwoc_results.tex")
    with open(p, "w", encoding="utf-8") as f:
        f.write("\\begin{tabular}{l r c c r r r r r r}\n\\hline\n")
        f.write("Water & $d$ (m) & $L$ & $\\lambda$ & $N$ & clicks & "
                "$n_{\\rm sift}$ & QBER & 95\\% u.b. & model \\\\\n\\hline\n")
        for r in rows:
            f.write("%s & %g & %d & %d & %d & %d & %d & %.2f\\%% & %.2f\\%% & "
                    "%.2f\\%% \\\\\n" % (
                        WATER_TYPES[r["water"]].name, r["distance_m"],
                        r["turb"], r["lam_nm"], r["n_qubit"], r["n_click"],
                        r["n_sift"], r["qber"] * 100, r["qber_hi2"] * 100,
                        r["qber_model"] * 100))
        f.write("\\hline\n\\end{tabular}\n")
    print("  ✓ table_uwoc_results.tex")


def summarise(rows):
    print("\n  ĐỐI CHIẾU FPGA ↔ MÔ HÌNH")
    tot_m = sum(r["p_click_model"] * r["n_qubit"] for r in rows)
    tot_c = sum(r["n_click"] for r in rows)
    print("    click đo được / click mô hình dự đoán = %d / %.0f = %.3f"
          % (tot_c, tot_m, tot_c / max(tot_m, 1e-9)))
    ok = [r for r in rows if r["valid"]]
    if ok:
        ns = sum(r["n_sift"] for r in ok)
        ne = sum(r["n_err"] for r in ok)
        lo, hi = cp_interval(ne, ns)
        print("    QBER gộp trên %d bit sàng = %.2f%%  [%.2f%%, %.2f%%]"
              % (ns, ne / ns * 100, lo * 100, hi * 100))
    print("    tổng thời gian đo = %.2f giờ, %d cấu hình"
          % (sum(r["seconds"] for r in rows) / 3600, len(rows)))


def main():
    global DATA, OUT
    import argparse
    ap = argparse.ArgumentParser(description="Dựng hình bài báo UWOC")
    ap.add_argument("--data", default=DATA, help="thư mục số liệu FPGA")
    ap.add_argument("--out", default=OUT, help="thư mục hình ra")
    a = ap.parse_args()

    DATA, OUT = os.path.abspath(a.data), os.path.abspath(a.out)
    os.makedirs(OUT, exist_ok=True)

    rows = load_points()
    print("Đã nạp %d điểm đo từ %s\n" % (len(rows), DATA))
    fig_pclick(rows)
    fig_qber_distance(rows)
    fig_skr(rows)
    fig_turbulence(rows)
    fig_wavelength(rows)
    fig_water_types(rows)
    fig_sample_size(rows)
    fig_adaptive(rows)
    fig_qber_vs_irradiance()
    fig_adaptive_controller(rows)
    fig_adaptive_timeline(rows)
    fig_adaptive_gain(rows)
    fig_block_dispersion()
    fig_adaptive_turbulence(rows)
    write_table(rows)
    summarise(rows)
    print("\nHình ở %s" % OUT)


if __name__ == "__main__":
    main()
