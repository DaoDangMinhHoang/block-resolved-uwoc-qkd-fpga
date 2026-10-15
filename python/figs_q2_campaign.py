#!/usr/bin/env python3
"""
figs_q2_campaign.py -- paper figures for the Q2 hardware campaign (rev E).

Consumes what phase2_security_analysis.py already wrote:

    phase2_results/security_summary.csv     one row per (turbulence, Y0) group
    phase2_results/bin_scan_nested.csv      nested resolution scan (--bin-scan)
    phase2_results/observable_label.csv     latent vs observable label
    data/q2/nom_L*.csv                      nominal blocks,  Y0 = 218/2^24
    data/q2/y0zero_L*.csv                   control blocks,  Y0 = 0

Writes to Images/:

    fig_uwoc_q2_defect_campaign.png   THE headline: pooled QBER is flat while
                                      Delta_cg spans two decades. Two stacked
                                      panels sharing the x axis -- never a
                                      second y axis on one plot. The model is
                                      the ROM-exact population defect at the
                                      SAME 8 bins as the measurement.
    fig_uwoc_q2_exact_vs_bkm.png      exact defect vs its BKM quadratic.
    fig_uwoc_q2_control_gt.png        (G,T) on log-log, L1, L3, L5 in one panel.
    fig_uwoc_q2_qber_vs_h.png         QBER against irradiance at L5.
    fig_uwoc_q2_qber_vs_h_levels.png  the same at L1, L3, L5 in one panel.
    fig_uwoc_q2_bin_scan.png          [rev E] nested resolution scan with
                                      bootstrap errors, up to every distinct h.
    fig_uwoc_q2_label.png             [rev E] latent fading word vs the
                                      protocol-observable label.
    fig_uwoc_q2_ablation.png          [rev E] matched-variance family pairs,
                                      only when rows 6/7 were measured.
    fig_uwoc_q2_decoy.png             [rev E] block-resolved vs pooled
                                      decoy-state rate, only when decoy
                                      captures were analysed.

Colours are categorical slots of the validated reference palette, in fixed
order: blue = Y0 nominal, orange = Y0 = 0 control -- except in the (G,T)
figure, where colour is the turbulence level and Y0 is carried by marker and
line style. Models are NOT series -- they are drawn in muted ink.

    python python/figs_q2_campaign.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.legend_handler import HandlerPatch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))

from phase2_security_analysis import (           # noqa: E402
    assign_bins, aggregate, SIGMA2_HO_20M, SEED,
    operating_point, moments_from_h, sample_h, V_priv,
)

FIG_DIR = ROOT / "Images"
RESULTS = ROOT / "phase2_results"
SUMMARY = RESULTS / "security_summary.csv"
Q2 = ROOT / "data" / "q2b"

# ---- validated reference palette, light surface ---------------------------
C_NOM = "#2a78d6"        # categorical slot 1
C_CTL = "#eb6834"        # categorical slot 2
INK = "#0b0b0b"          # primary ink
INK2 = "#52514e"         # secondary ink
MUTED = "#898781"        # axis / labels / non-series reference lines
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"

LEVELS = [1, 2, 3, 4, 5]

# ---- single-column authoring ----------------------------------------------
# main.tex places every figure at \columnwidth (~3.5 in in IEEEtran two-column),
# so the figures are authored at that width and dropped in at 1:1. Point sizes
# below are therefore the sizes that actually appear in the PDF: keep them at
# 6 pt or above.
COLW = 3.5
FS_TITLE = 8.0
FS_LABEL = 7.0
FS_TICK = 6.0
FS_LEGEND = 5.8
FS_NOTE = 5.8
FS_PANEL = 7.5

# (G,T) figure: one panel, colour = turbulence level. Categorical slots 1-3,
# the only three that validate all-pairs, so every level stays distinct.
GT_LEVELS = (1, 3, 5)
GT_COLOR = {1: "#2a78d6", 3: "#eb6834", 5: "#1baf7a", 4: "#8a5cc9"}
# QBER-vs-h figure: label offsets (points), ha, va for levels whose deep-fade
# end sits on the shared model curve; the rest take the default above-left.
QH_LABEL = {1: ((0, 26), "center", "bottom"), 3: ((-16, -20), "right", "top")}


def _swatch(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
    """Square colour key centred in the handle box. A Patch otherwise fills the
    whole (long) handle and reads as one more, thicker, solid line."""
    s = 0.8 * fontsize
    return Rectangle((width / 2 - s / 2 - xdescent, height / 2 - s / 2 - ydescent),
                     s, s)


# Pass as legend(handler_map=...) wherever Patch handles key the colours.
SWATCH = {Patch: HandlerPatch(patch_func=_swatch)}


def _panel(ax, tag):
    """Panel tag above the axes, left-aligned. Takes the place of a per-panel
    title: what the panel shows belongs in the LaTeX caption, not in the PNG."""
    ax.text(0.0, 1.015, tag, transform=ax.transAxes, color=INK,
            fontsize=FS_PANEL, fontweight="bold", ha="left", va="bottom")


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=FS_TICK)
    for lb in ax.get_xticklabels() + ax.get_yticklabels():
        lb.set_color(INK2)


def new_fig(*a, **k):
    fig, ax = plt.subplots(*a, **k)
    fig.patch.set_facecolor(SURFACE)
    for x in (ax.ravel() if isinstance(ax, np.ndarray) else [ax]):
        style(x)
    return fig, ax


def model_curve(distance=25.0, n=400000):
    """Continuous channel model at FULL resolution. Kept for reference only:
    comparing it with an 8-bin measurement mixes two resolutions."""
    rng = np.random.default_rng(SEED)
    op = operating_point(distance, 450)
    out = {}
    for lvl, s2 in SIGMA2_HO_20M.items():
        h = sample_h(s2, n, rng)
        G, T = moments_from_h(op, h)
        V0 = float(V_priv(np.array([G.mean()]), np.array([T.mean()]))[0])
        out[lvl] = (float(V_priv(G, T).mean()) - V0) / V0
    return out


def load_blocks(path):
    d = pd.read_csv(path, low_memory=False)
    d = d[d["capture_valid"] == 1].copy()
    d["trials"] = d["trials_z"] + d["trials_x"]
    d["clicks"] = d["clicks_z"] + d["clicks_x"]
    d["errors"] = d["errors_z"] + d["errors_x"]
    d["h_q8"] = d["h_mean"] / 256.0
    return d


def binned(path, n_bins=8):
    return aggregate(assign_bins(load_blocks(path), n_bins))


def errbars(sub, col, lo, hi):
    return np.abs(np.vstack([sub[col] - sub[lo], sub[hi] - sub[col]]))


def _main_levels(s):
    return s[s.turb.between(1, 5) & (s.get("decoy_enable", 0) == 0)]


# ===========================================================================
# 1. Headline: pooled QBER flat, defect over two decades
# ===========================================================================

def fig_campaign(s):
    s = _main_levels(s)
    nom = s[s.p_noise == 218].sort_values("turb")
    ctl = s[s.p_noise == 0].sort_values("turb")

    fig, (a0, a1) = plt.subplots(2, 1, figsize=(COLW, 3.7), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1.7]})
    fig.patch.set_facecolor(SURFACE)
    style(a0)
    style(a1)

    # ---- top: what a pooled campaign reports -----------------------------
    q = 100 * nom["QBER_pool"].to_numpy()
    n_sift = (nom["G_pool"] * nom["n_blocks"] * 65536).to_numpy()
    qerr = 100 * 1.96 * np.sqrt(nom["QBER_pool"] * (1 - nom["QBER_pool"]) / n_sift)
    a0.errorbar(nom["turb"], q, yerr=qerr, marker="o", ms=8, lw=2, capsize=4,
                color=C_NOM, zorder=3)
    a0.set_ylim(3.35, 4.05)
    a0.set_ylabel("pooled QBER (%)", color=INK2, fontsize=FS_LABEL)

    # ---- bottom: what the security functional reports --------------------
    r = nom[nom.resolved == 1]
    u = nom[nom.resolved == 0]
    a1.errorbar(r["turb"], r["excess_over_V"],
                yerr=errbars(r, "excess_over_V", "ci_lo_over_V", "ci_hi_over_V"),
                marker="o", ms=8, lw=2, capsize=4, color=C_NOM, zorder=4,
                label=r"$\Delta_{cg}/V$, $Y_0$ nominal")
    if len(u):
        a1.plot(u["turb"], u["upper_limit95_over_V"], "v", ms=9, mfc="none",
                mew=1.8, color=C_NOM, zorder=4,
                label="95% upper limit (unresolved)")

    rc = ctl[ctl.resolved == 1]
    uc = ctl[ctl.resolved == 0]
    if len(rc):
        a1.errorbar(rc["turb"], rc["excess_over_V"],
                    yerr=errbars(rc, "excess_over_V", "ci_lo_over_V", "ci_hi_over_V"),
                    marker="s", ms=8, lw=2, ls="--", capsize=4, color=C_CTL,
                    zorder=4, label=r"$Y_0 = 0$ radial-null control")
    if len(uc):
        a1.plot(uc["turb"], uc["upper_limit95_over_V"], "v", ms=9, mfc="none",
                mew=1.8, color=C_CTL, zorder=4,
                label="control, 95% upper limit")

    if "model_rom_bins_over_V" in nom.columns:
        a1.plot(nom["turb"], nom["model_rom_bins_over_V"], ls=":", lw=1.6,
                color=MUTED, marker="D", ms=5, mfc=SURFACE, zorder=2,
                label="ROM-exact model, same 8 bins")
    else:
        model = model_curve()
        a1.plot(LEVELS, [model[k] for k in LEVELS], ls=":", lw=1.6, color=MUTED,
                marker="D", ms=5, mfc=SURFACE, zorder=2,
                label="analytical channel model")

    a1.set_yscale("log")
    a1.set_xticks(LEVELS)
    a1.set_xlim(0.6, 5.6)
    a1.set_xlabel("ocean turbulence level", color=INK2, fontsize=FS_LABEL)
    a1.set_ylabel(r"$\Delta_{cg}\,/\,V$   (counting-noise null subtracted)",
                  color=INK2, fontsize=FS_LABEL)
    hs, ls_ = a1.get_legend_handles_labels()
    order = sorted(range(len(ls_)),
                   key=lambda i: ("upper limit" in ls_[i], "model" in ls_[i]))
    a1.legend([hs[i] for i in order], [ls_[i] for i in order],
              fontsize=FS_LEGEND, loc="upper left", framealpha=0.95)

    span = float(nom["excess_over_V"].iloc[-1] / nom["excess_over_V"].iloc[0])
    a1.annotate(r"$\times\,%.0f$" % span + "\nfrom L1 to L5",
                xy=(5.42, 6e-6), ha="right", fontsize=FS_NOTE, color=INK2)

    fig.tight_layout()
    p = FIG_DIR / "fig_uwoc_q2_defect_campaign.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


# ===========================================================================
# 2. Exact defect against its BKM quadratic
# ===========================================================================

def fig_bkm(s):
    s = _main_levels(s)
    nom = s[(s.p_noise == 218) & (s.resolved == 1)].sort_values("turb")

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(COLW, 2.0),
                                 gridspec_kw={"width_ratios": [1.35, 1]})
    fig.patch.set_facecolor(SURFACE)
    style(a0)
    style(a1)

    a0.errorbar(nom["turb"], nom["excess_over_V"],
                yerr=errbars(nom, "excess_over_V", "ci_lo_over_V", "ci_hi_over_V"),
                marker="o", ms=8, lw=2, capsize=4, color=C_NOM, zorder=4,
                label=r"exact  $\Delta_{cg}/V$")
    a0.errorbar(nom["turb"], nom["bkm_excess_over_V"],
                yerr=errbars(nom, "bkm_excess_over_V",
                             "bkm_ci_lo_over_V", "bkm_ci_hi_over_V"),
                marker="s", ms=8, lw=2, ls="--", capsize=4, color=C_CTL,
                zorder=3,
                label=r"BKM  $\frac{1}{2}\mathrm{Tr}(\Gamma_0\mathrm{Cov})/V$")
    a0.set_yscale("log")
    a0.set_xticks(sorted(nom["turb"]))
    a0.set_xlabel("ocean turbulence level", color=INK2, fontsize=FS_LABEL)
    a0.set_ylabel(r"$\Delta\,/\,V$", color=INK2, fontsize=FS_LABEL)
    # The curve runs bottom-left to top-right, so the lower-right corner is the
    # only one the data never enters at single-column size.
    a0.legend(fontsize=FS_LEGEND, loc="lower right", framealpha=0.95)
    _panel(a0, "(a)")

    ratio = (nom["bkm_excess_over_V"] / nom["excess_over_V"]).to_numpy()
    a1.axhline(1.0, color=MUTED, ls=":", lw=1.6, zorder=2)
    a1.plot(nom["turb"], ratio, marker="o", ms=8, lw=2, color=C_NOM, zorder=3)
    for x, y in zip(nom["turb"], ratio):
        a1.annotate("%.2f" % y, xy=(x, y), xytext=(0, 9),
                    textcoords="offset points", ha="center",
                    fontsize=FS_NOTE, color=INK2)
    # No "quadratic exact" tag: at one column its box spans levels 2-4 and lands
    # on the 0.95 / 0.93 value labels. The dotted line at 1.0 plus the axis
    # label already say it, and the caption repeats it.
    a1.set_ylim(0.3, 1.18)
    a1.set_xticks(sorted(nom["turb"]))
    a1.set_xlabel("ocean turbulence level", color=INK2, fontsize=FS_LABEL)
    a1.set_ylabel("BKM / exact", color=INK2, fontsize=FS_LABEL)
    _panel(a1, "(b)")

    fig.tight_layout()
    p = FIG_DIR / "fig_uwoc_q2_exact_vs_bkm.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


# ===========================================================================
# 3. The mechanism: (G,T) on log-log
# ===========================================================================

def _gt_files(L):
    """Nominal and control captures for one level. L1 was re-run at 286 000
    blocks under a _full suffix (the budget security_summary.csv reports), so
    prefer that pair whenever both halves of it exist."""
    full = (Q2 / ("nom_L%d_full.csv" % L), Q2 / ("y0zero_L%d_full.csv" % L))
    if all(p.exists() for p in full):
        return full
    return Q2 / ("nom_L%d.csv" % L), Q2 / ("y0zero_L%d.csv" % L)


def fig_control_gt(y0_nominal=218 / 2 ** 24, levels=GT_LEVELS):
    fig, ax = new_fig(figsize=(COLW, 2.9))

    runs = []
    trials = attempts = ctl_t = ctl_g = 0.0
    for L in levels:
        pn, pc = _gt_files(L)
        dn = load_blocks(pn)
        bn = aggregate(assign_bins(dn, 8))
        bc = binned(pc)
        bn = bn[(bn.G > 0) & (bn["T"] > 0)]
        bc = bc[(bc.G > 0) & (bc["T"] > 0)]
        runs.append((L, bn, bc))
        trials += float(dn["trials"].sum())
        attempts += float(dn["attempts"].sum())
        ctl_t += float(bc["T"].sum())
        ctl_g += float(bc["G"].sum())

    q_sift = trials / attempts
    g_inf, t_inf = q_sift * y0_nominal, q_sift * 0.5 * y0_nominal
    e0 = ctl_t / ctl_g

    g_lo = min(min(bn.G.min(), bc.G.min()) for _, bn, bc in runs)
    g_hi = max(max(bn.G.max(), bc.G.max()) for _, bn, bc in runs)
    gg = np.logspace(np.log10(g_lo * 0.5), np.log10(g_hi * 1.6), 50)
    ax.plot(gg, e0 * gg, ls=":", lw=1.6, color=MUTED, zorder=2)
    ax.axhline(t_inf, color=MUTED, lw=1, ls="--", zorder=1)
    ax.axvline(g_inf, color=MUTED, lw=1, ls="--", zorder=1)

    for k, (L, bn, bc) in enumerate(reversed(runs)):
        c = GT_COLOR[L]
        ax.plot(bc.G, bc["T"], marker="s", ms=7, lw=1.8, ls="--", color=c,
                mfc=SURFACE, mew=1.8, zorder=3 + 2 * k)
        ax.plot(bn.G, bn["T"], marker="o", ms=8, lw=2, color=c, mec=SURFACE,
                mew=1.2, zorder=4 + 2 * k)
        lo = bn.sort_values("G").iloc[0]
        ax.annotate("L%d" % L, xy=(float(lo.G), float(lo["T"])),
                    xytext=(-5, 7), textcoords="offset points", ha="right",
                    va="bottom", fontsize=FS_NOTE, color=INK2, fontweight="bold",
                    zorder=20)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("detection gain $G$ (per attempt)", color=INK2, fontsize=FS_LABEL)
    ax.set_ylabel("error gain $T$ (per attempt)", color=INK2, fontsize=FS_LABEL)
    ax.annotate(r"$q\,Y_0/2$", xy=(gg[0], t_inf), xytext=(2, 4),
                textcoords="offset points", fontsize=FS_LEGEND, color=MUTED)
    ax.annotate(r"$q\,Y_0$", xy=(g_inf, ax.get_ylim()[0]), xytext=(4, 8),
                textcoords="offset points", fontsize=FS_LEGEND, color=MUTED)
    lo5 = runs[-1][1].sort_values("G").iloc[0]
    ax.annotate("deep fades leave the ray\nonly when $Y_0>0$",
                xy=(float(lo5.G) * 1.15, float(lo5["T"])),
                xytext=(9e-5, 6e-7), fontsize=FS_NOTE, color=INK2,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.2))

    handles = [Line2D([], [], ls=":", lw=1.6, color=MUTED,
                      label=r"ray $T=e_{\rm pol}G$  (slope 1)")]
    handles += [Patch(facecolor=GT_COLOR[L], edgecolor="none",
                      label="turbulence level %d" % L) for L in levels]
    handles += [
        Line2D([], [], marker="o", ms=8, lw=2, color=INK2, mec=SURFACE,
               mew=1.2, label=r"$Y_0$ nominal"),
        Line2D([], [], marker="s", ms=7, lw=1.8, ls="--", color=INK2,
               mfc=SURFACE, mew=1.8, label=r"$Y_0 = 0$ control"),
    ]
    ax.legend(handles=handles, fontsize=FS_LEGEND, loc="upper left", framealpha=0.95,
              handlelength=3.4, handler_map=SWATCH)

    fig.tight_layout()
    p = FIG_DIR / "fig_uwoc_q2_control_gt.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


# ===========================================================================
# 4. QBER against irradiance
# ===========================================================================

def fig_qber_h():
    op_n = operating_point(25.0, 450)
    op_0 = operating_point(25.0, 450, y0=0.0)

    fig, ax = new_fig(figsize=(COLW, 2.6))
    hh = np.logspace(-2.4, 1.2, 200)
    Gn, Tn = moments_from_h(op_n, hh)
    G0, T0 = moments_from_h(op_0, hh)
    ax.plot(hh, 100 * Tn / Gn, lw=1.6, ls=":", color=MUTED, zorder=2,
            label="model, $Y_0$ nominal")
    ax.plot(hh, 100 * T0 / G0, lw=1.6, ls="--", color=MUTED, zorder=2,
            label="model, $Y_0=0$  ($=e_{\\rm pol}$)")

    for L, mk in ((5, "o"),):
        bn = binned(Q2 / ("nom_L%d.csv" % L))
        bc = binned(Q2 / ("y0zero_L%d.csv" % L))
        bn = bn[bn.G > 0]
        bc = bc[bc.G > 0]
        ax.plot(bn.h, 100 * bn["T"] / bn.G, marker=mk, ms=8, lw=2, color=C_NOM,
                zorder=4, label=r"FPGA L5, $Y_0$ nominal")
        ax.plot(bc.h, 100 * bc["T"] / bc.G, marker="s", ms=8, lw=2, ls="--",
                color=C_CTL, zorder=4, label=r"FPGA L5, $Y_0=0$ control")

    ax.set_xscale("log")
    ax.set_xlabel("block irradiance $h$   (12-bit, 256 $\\equiv$ 1.0)",
                  color=INK2, fontsize=FS_LABEL)
    ax.set_ylabel("block QBER (%)", color=INK2, fontsize=FS_LABEL)
    ax.set_title("With the background off, QBER is independent of fading",
                 color=INK, fontsize=FS_TITLE, pad=10)
    ax.legend(fontsize=FS_NOTE, framealpha=0.95)
    fig.tight_layout()
    p = FIG_DIR / "fig_uwoc_q2_qber_vs_h.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


def fig_qber_h_levels(levels=GT_LEVELS):
    op_n = operating_point(25.0, 450)
    op_0 = operating_point(25.0, 450, y0=0.0)

    fig, ax = new_fig(figsize=(COLW, 2.8))
    hh = np.logspace(-2.4, 1.2, 200)
    Gn, Tn = moments_from_h(op_n, hh)
    G0, T0 = moments_from_h(op_0, hh)
    ax.plot(hh, 100 * Tn / Gn, lw=1.6, ls=":", color=MUTED, zorder=2)
    ax.plot(hh, 100 * T0 / G0, lw=1.6, ls="--", color=MUTED, zorder=2)

    runs = []
    for L in levels:
        pn, pc = _gt_files(L)
        bn = binned(pn)
        bc = binned(pc)
        runs.append((L, bn[bn.C > 0], bc[bc.C > 0]))

    for k, (L, bn, bc) in enumerate(reversed(runs)):
        c = GT_COLOR[L]
        for b, nominal in ((bc, False), (bn, True)):
            q = (b.E / b.C).to_numpy()
            err = 100 * 1.96 * np.sqrt(q * (1 - q) / b.C.to_numpy())
            if nominal:
                kw = dict(marker="o", ms=7.5, lw=2, ls="-", mec=SURFACE,
                          mew=1.2, zorder=4 + 2 * k)
            else:
                kw = dict(marker="s", ms=6.5, lw=1.8, ls="--", mfc=SURFACE,
                          mew=1.8, zorder=3 + 2 * k)
            ax.errorbar(b.h, 100 * q, yerr=err, color=c, elinewidth=1.1,
                        capsize=0, **kw)
        lo = bn.sort_values("h").iloc[0]
        off, ha, va = QH_LABEL.get(L, ((-5, 7), "right", "bottom"))
        ax.annotate("L%d" % L, xy=(float(lo.h), 100 * float(lo.E / lo.C)),
                    xytext=off, textcoords="offset points", ha=ha, va=va,
                    fontsize=FS_NOTE, color=INK2, fontweight="bold", zorder=20,
                    arrowprops=None if L not in QH_LABEL else
                    dict(arrowstyle="-", color=MUTED, lw=0.9, shrinkB=5))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1.8, 48)
    ax.set_yticks([2, 3, 4, 5, 7, 10, 20, 30, 40])
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.tick_params(labelcolor=INK2)
    ax.set_xlabel("block irradiance $h$   (12-bit, 256 $\\equiv$ 1.0)",
                  color=INK2, fontsize=FS_LABEL)
    ax.set_ylabel("block QBER (%)", color=INK2, fontsize=FS_LABEL)

    handles = [
        Line2D([], [], ls=":", lw=1.6, color=MUTED, label="model, $Y_0$ nominal"),
        Line2D([], [], ls="--", lw=1.6, color=MUTED,
               label="model, $Y_0=0$  ($=e_{\\rm pol}$)"),
    ]
    handles += [Patch(facecolor=GT_COLOR[L], edgecolor="none",
                      label="turbulence level %d" % L) for L in levels]
    handles += [
        Line2D([], [], marker="o", ms=7.5, lw=2, color=INK2, mec=SURFACE,
               mew=1.2, label="FPGA, $Y_0$ nominal"),
        Line2D([], [], marker="s", ms=6.5, lw=1.8, ls="--", color=INK2,
               mfc=SURFACE, mew=1.8, label="FPGA, $Y_0=0$ control"),
    ]
    ax.legend(handles=handles, fontsize=FS_LEGEND, loc="upper right",
              framealpha=0.95, handlelength=3.4, handler_map=SWATCH)

    fig.tight_layout()
    p = FIG_DIR / "fig_uwoc_q2_qber_vs_h_levels.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


# ===========================================================================
# 5. [rev E] Resolution: nested partitions up to every distinct h
# ===========================================================================

SCAN_X = {2: 2, 4: 4, 8: 8, 16: 16, 32: 32, 64: 64, 128: 128, "all": 320}


def fig_bin_scan(levels=(1, 3, 4, 5)):
    """Null-corrected defect against the number of nested fading bins, with
    bootstrap standard errors. Refinement can only raise the population value
    along this chain; the last point groups blocks with identical fading word,
    i.e. full block resolution. One IEEE column wide."""
    path = RESULTS / "bin_scan_nested.csv"
    if not path.exists():
        return None
    d = pd.read_csv(path)
    d = d[(d.p_noise == 218) & (d.get("decoy_enable", 0) == 0)]
    fig, ax = new_fig(figsize=(3.5, 2.8))
    for L in levels:
        g = d[d.turb == L].copy()
        if g.empty:
            continue
        g["x"] = [SCAN_X["all" if str(b) == "all" else int(b)] for b in g.bins]
        g = g.sort_values("x")
        c = GT_COLOR[L]
        ax.errorbar(g.x, g.excess_over_V, yerr=g.se_over_V, marker="o", ms=4,
                    lw=1.4, capsize=2, color=c, mec=SURFACE, mew=0.6, zorder=4)
        ax.annotate("L%d" % L, xy=(g.x.iloc[-1], g.excess_over_V.iloc[-1]),
                    xytext=(5, 0), textcoords="offset points", va="center",
                    fontsize=FS_LEGEND, color=INK2, fontweight="bold")
    ax.axvline(8, color=MUTED, lw=1, ls="--", zorder=1)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([2, 8, 32, 128, 320])
    ax.set_xticklabels(["2", "8", "32", "128", "all h"])
    ax.set_xlim(1.6, 520)
    ax.tick_params(labelsize=FS_TICK, labelcolor=INK2)
    ax.annotate("reported", xy=(8, 0), xycoords=("data", "axes fraction"),
                xytext=(3, 3), textcoords="offset points", ha="left",
                va="bottom", fontsize=FS_LEGEND, color=MUTED)
    ax.set_xlabel("nested fading bins", color=INK2, fontsize=FS_LEGEND)
    ax.set_ylabel(r"$\Delta_{cg}\,/\,V$ (null corrected)", color=INK2, fontsize=FS_LEGEND)
    fig.tight_layout(pad=0.4)
    p = FIG_DIR / "fig_uwoc_q2_bin_scan.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


def fig_label(levels=(4, 5)):
    """Latent fading word against the protocol-observable label (clicks of
    basis-mismatched attempts) at the same number of bins."""
    path = RESULTS / "observable_label.csv"
    if not path.exists():
        return None
    d = pd.read_csv(path)
    # The dose-response captures share their turbulence row with the nominal
    # ones, so without this filter L3 and L5 would mix seven Y0 settings.
    if "p_noise" in d.columns:
        d = d[d.p_noise == 218]
    if "decoy_enable" in d.columns:
        d = d[d.decoy_enable == 0]
    fig, ax = new_fig(figsize=(3.5, 2.6))
    for L in levels:
        c = GT_COLOR[L]
        for lab, ls, mk in (("latent_h", "-", "o"), ("observable_mismatch_clicks", "--", "s")):
            g = d[(d.turb == L) & (d.block_label == lab)].sort_values("bins")
            if g.empty:
                continue
            ax.errorbar(g.bins, g.excess_over_V, yerr=g.se_over_V, marker=mk,
                        ms=4, lw=1.4, ls=ls, capsize=2, color=c,
                        mfc=SURFACE if lab != "latent_h" else c, zorder=4)
        gl = d[(d.turb == L) & (d.block_label == "latent_h")].sort_values("bins")
        if len(gl):
            ax.annotate("L%d" % L, xy=(gl.bins.iloc[-1], gl.excess_over_V.iloc[-1]),
                        xytext=(5, 0), textcoords="offset points", va="center",
                        fontsize=FS_LEGEND, color=INK2, fontweight="bold")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([2, 4, 8, 16])
    ax.set_xticklabels(["2", "4", "8", "16"])
    ax.tick_params(labelsize=FS_TICK, labelcolor=INK2)
    ax.set_xlabel("bins", color=INK2, fontsize=FS_LEGEND)
    ax.set_ylabel(r"$\Delta_{cg}\,/\,V$ (null corrected)", color=INK2, fontsize=FS_LEGEND)
    handles = [Line2D([], [], marker="o", ls="-", color=INK2, ms=4, label="latent fading word"),
               Line2D([], [], marker="s", ls="--", color=INK2, ms=4, mfc=SURFACE,
                      label="observable label")]
    ax.legend(handles=handles, fontsize=FS_LEGEND, loc="lower right", framealpha=0.95)
    fig.tight_layout(pad=0.4)
    p = FIG_DIR / "fig_uwoc_q2_label.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


def fig_ablation(s):
    """Matched-variance family pairs, only when rows 6/7 were measured."""
    nom = s[(s.p_noise == 218) & (s.turb.isin([3, 4, 6, 7]))]
    if not set(nom.turb) >= {6, 7}:
        return None
    pairs = [("$\\sigma^2_{h_o}=0.30$", 3, 6), ("$\\sigma^2_{h_o}=1.00$", 7, 4)]
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(COLW, 1.9))
    fig.patch.set_facecolor(SURFACE)
    style(a0)
    style(a1)
    for x, (name, ln, wb) in enumerate(pairs):
        for dx, row, c in ((-0.12, ln, C_NOM), (0.12, wb, C_CTL)):
            r = nom[nom.turb == row]
            if r.empty:
                continue
            r = r.iloc[0]
            a0.errorbar(x + dx, r.excess_over_V,
                        yerr=[[r.excess_over_V - r.ci_lo_over_V], [r.ci_hi_over_V - r.excess_over_V]],
                        marker="o", color=c, capsize=3)
            a1.plot(x + dx, r.bkm_excess_over_V / r.excess_over_V, "o", color=c)
    for ax in (a0, a1):
        ax.set_xticks([0, 1])
        ax.set_xticklabels([p[0] for p in pairs])
        # Without this the +-0.12 offsets land on the panel spines.
        ax.set_xlim(-0.5, 1.5)
    a0.set_yscale("log")
    a0.set_ylabel(r"$\Delta_{cg}/V$", color=INK2)
    a1.set_ylabel("BKM / exact", color=INK2)
    a1.set_ylim(0.3, 1.1)
    a0.legend(handles=[Patch(facecolor=C_NOM, label="log-normal"),
                       Patch(facecolor=C_CTL, label="Weibull")],
              fontsize=FS_LEGEND, handler_map=SWATCH)
    _panel(a0, "(a)")
    _panel(a1, "(b)")
    fig.tight_layout()
    p = FIG_DIR / "fig_uwoc_q2_ablation.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


def fig_decoy():
    path = RESULTS / "decoy_rates.csv"
    if not path.exists():
        return None
    d = pd.read_csv(path)
    # A capture is identified by (turb, p_noise): three of the six operating
    # points share turb = 5 and differ only in the background floor, so
    # grouping by turb alone silently merges them.
    d = d[d["label"] == "h"]
    pts = [((1, 218), "L1, $Y_0$", "#2a78d6"),
           ((3, 218), "L3, $Y_0$", "#eb6834"),
           ((5, 218), "L5, $Y_0$", "#1baf7a"),
           ((5, 872), r"L5, $4Y_0$", "#8a5cc9"),
           ((4, 3488), r"L4, $16Y_0$", "#c7a008"),
           ((5, 3488), r"L5, $16Y_0$", "#c4384f")]

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(COLW, 2.1))
    fig.patch.set_facecolor(SURFACE)
    style(a0)
    style(a1)
    for (turb, pn), name, col in pts:
        g0 = d[(d.turb == turb) & (d.p_noise == pn)]
        if g0.empty:
            continue
        for mode, ls, ax in (("asymptotic", "-", a0), ("finite", "--", a1)):
            g = g0[g0["mode"] == mode].sort_values("bins")
            ax.plot(g.bins, g.gain_over_pooled, marker="o", ms=4, lw=1.5,
                    ls=ls, color=col, mec=SURFACE, mew=0.6, zorder=4,
                    label=name)
    for ax, tag in ((a0, "(a)"), (a1, "(b)")):
        ax.axhline(1.0, color=MUTED, ls=":", lw=1.4, zorder=1)
        ax.set_xscale("log", base=2)
        ax.set_xticks([1, 2, 4, 8, 16, 32])
        ax.set_xticklabels(["1", "2", "4", "8", "16", "32"])
        ax.set_xlabel("fading bins (1 = pooled)", color=INK2, fontsize=FS_NOTE)
        ax.tick_params(labelsize=FS_TICK, labelcolor=INK2)
        _panel(ax, tag)
    a0.set_ylabel("key rate / pooled rate", color=INK2, fontsize=FS_NOTE)
    a1.set_ylabel("key rate / pooled rate", color=INK2, fontsize=FS_NOTE)
    # One shared key under both panels: inside panel (a) it would sit on top of
    # the steep 16*Y0 curve it is meant to identify.
    hs, ls_ = a0.get_legend_handles_labels()
    # Six entries in one row do not fit across a single column: the last one is
    # clipped. Two rows of three.
    fig.legend(hs, ls_, loc="lower center", ncol=3, fontsize=FS_LEGEND,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(pad=0.6, rect=(0, 0.17, 1, 1))
    p = FIG_DIR / "fig_uwoc_q2_decoy.png"
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


def main():
    FIG_DIR.mkdir(exist_ok=True)
    if not SUMMARY.exists():
        raise SystemExit("run phase2_security_analysis.py first: %s missing"
                         % SUMMARY)
    s = pd.read_csv(SUMMARY)
    for p in (fig_campaign(s), fig_bkm(s), fig_control_gt(),
              fig_qber_h(), fig_qber_h_levels(), fig_bin_scan(), fig_label(),
              fig_ablation(s), fig_decoy()):
        if p is not None:
            print("wrote %s" % p.relative_to(ROOT))


if __name__ == "__main__":
    main()
