#!/usr/bin/env python3
"""
fig1_block_diagram.py -- Fig. 1, the rev B testbed architecture.

Redraws QKD_UWOC_Block in the style of the original PowerPoint export (flat
pastel zones, solid blocks, short arrow labels) but with the rev B content:

  * the TRNG is labelled as serving the OOK path ONLY  (trng_advance is
    asserted only in FSM_NEXT of the per-qubit FSM, top_module.v:1184),
  * the real-time attempt engine and the two security_fast_mode multiplexers
    (top_module.v:945 and 668-673) are drawn explicitly,
  * the four block-acquisition modules get their own zone,
  * the emulator shows its 8-row h_o ROM, its five xorshift64 generators and
    the three interleaved intensity classes.

Every signal name on an arrow is a real wire in top_module.v.
Geometry is in "display units": 100 units = 1 inch of the drawn figure, so the
font sizes below land at roughly the same final point size as the original
export once the figure is placed at \\textwidth.

    python python/fig1_block_diagram.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
OUTS = [ROOT / "Images", ROOT / "Overleaf" / "figures"]

W, H = 1500.0, 720.0          # display units
DPI = 400

# ---- palette, matched to the original export -------------------------------
Z_GREEN = ("#e8f5e9", "#66bb6a")
Z_BLUE = ("#e3f2fd", "#64b5f6")
Z_PEACH = ("#fff3e0", "#ffb74d")
Z_AMBER = ("#fff8e1", "#ffa726")
Z_PINK = ("#f9e3f5", "#ce93d8")
Z_GREY = ("#f4f4f2", "#bdbdbd")

B_GREEN = "#a5d6a7"
B_TEAL = "#4db6ac"
B_LBLUE = "#90caf9"
B_MBLUE = "#42a5f5"
B_DBLUE = "#1565c0"
B_ORANGE = "#ef9a6a"
B_AMBER = "#ffcc80"
B_GREY = "#d7d7d7"
B_RED = "#ef9a9a"
B_INNER = "#cfe6fb"

A_RAND = "#43a047"
A_DATA = "#0d47a1"
A_FAST = "#00897b"
A_EVT = "#546e7a"
A_BLK = "#ef6c00"
A_CTRL = "#b25000"
A_CFG = "#8e24aa"

INK = "#111111"

FS_ZONE = 15
FS_TITLE = 14
FS_SUB = 10
FS_INNER = 9.5
FS_ARROW = 9


def zone(ax, x0, y0, x1, y1, title, colors, title_dx=0.0):
    fc, ec = colors
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=14",
                                fc=fc, ec=ec, lw=1.8, zorder=1))
    if title:
        ax.text(x0 + 6 + title_dx, y1 + 9, title, fontsize=FS_ZONE,
                fontweight="bold", color=INK, ha="left", va="bottom", zorder=6)


def box(ax, x0, y0, x1, y1, title, sub=None, fc=B_LBLUE, tc=INK,
        fs_title=FS_TITLE, fs_sub=FS_SUB, round_size=8, lw=0.0, ec="none"):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=%g" % round_size,
                                fc=fc, ec=ec, lw=lw, zorder=3))
    cx = 0.5 * (x0 + x1)
    if sub:
        ax.text(cx, 0.5 * (y0 + y1) + 0.16 * (y1 - y0), title, fontsize=fs_title,
                fontweight="bold", color=tc, ha="center", va="center", zorder=5)
        ax.text(cx, 0.5 * (y0 + y1) - 0.20 * (y1 - y0), sub, fontsize=fs_sub,
                color=tc, ha="center", va="center", zorder=5, linespacing=1.35)
    else:
        ax.text(cx, 0.5 * (y0 + y1), title, fontsize=fs_title, fontweight="bold",
                color=tc, ha="center", va="center", zorder=5, linespacing=1.3)


def inner(ax, x0, y0, x1, y1, text):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=B_INNER,
                           ec="#7fb8e8", lw=0.9, zorder=4))
    ax.text(0.5 * (x0 + x1), 0.5 * (y0 + y1), text, fontsize=FS_INNER,
            color=INK, ha="center", va="center", zorder=5, linespacing=1.3)


def arrow(ax, pts, color, label=None, lpos=None, lha="center", lva="bottom",
          dashed=False, lw=2.0, rot=0):
    """Orthogonal polyline with an arrowhead on the last segment."""
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.add_line(Line2D(xs[:-1] + [xs[-1]], ys[:-1] + [ys[-1]], color=color,
                       lw=lw, zorder=7, solid_capstyle="round",
                       linestyle=(0, (5, 3)) if dashed else "-"))
    ax.annotate("", xy=pts[-1], xytext=pts[-2], zorder=7,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                shrinkA=0, shrinkB=0,
                                mutation_scale=15))
    if label:
        lx, ly = lpos
        ax.text(lx, ly, label, fontsize=FS_ARROW, color=color, ha=lha, va=lva,
                zorder=8, rotation=rot, linespacing=1.25)


def main():
    fig = plt.figure(figsize=(W / 100.0, H / 100.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ================= zones =================
    zone(ax, 20, 105, 250, 665, "ATTEMPT SOURCES", Z_GREEN)
    zone(ax, 265, 105, 870, 665, "QKD DATAPATH", Z_BLUE)
    zone(ax, 950, 380, 1200, 665, "MONITOR & ADAPT", Z_PEACH)
    zone(ax, 950, 105, 1480, 355, "REAL-TIME BLOCK ACQUISITION", Z_AMBER)
    zone(ax, 1215, 380, 1480, 665, "", Z_GREY)
    zone(ax, 20, 18, 1480, 80, "", Z_PINK)

    # ================= host strip =================
    ax.text(750, 58, "BB84 CONTROLLER  &  HOST INTERFACE", fontsize=16,
            fontweight="bold", color=INK, ha="center", va="center", zorder=5)
    ax.text(750, 33, "UART RX (commands)   ·   UART TX arbiter, "
                     "block reporter has priority   ·   460.8 kbaud",
            fontsize=9.5, color="#4a3050", ha="center", va="center", zorder=5)

    # ================= zone 1 =================
    box(ax, 35, 545, 235, 650, "TRNG ×3",
        "ring oscillator\nOOK path only", fc=B_GREEN)
    box(ax, 35, 375, 235, 500, "REAL-TIME\nATTEMPT ENGINE",
        "10 MHz pacer · xorshift64 ×3\n607 LE / 322 reg", fc=B_TEAL, tc="white",
        fs_title=13)
    box(ax, 70, 230, 200, 310, "MUX", "real-time mode", fc=B_RED,
        fs_title=13, fs_sub=9)

    # ================= zone 2 =================
    box(ax, 280, 560, 390, 650, "ALICE", "BB84 encoder", fc=B_GREEN, fs_title=13)
    box(ax, 280, 440, 390, 520, "OOK TX", "4-slot serializer", fc=B_MBLUE,
        tc="white", fs_title=13, fs_sub=9)

    box(ax, 420, 195, 700, 650, "", fc=B_LBLUE, round_size=10)
    ax.text(560, 625, "UWOC Channel emulator", fontsize=14, fontweight="bold",
            color=INK, ha="center", va="center", zorder=5)
    inner(ax, 435, 545, 685, 597, "path-loss & polarization tables")
    inner(ax, 435, 480, 685, 535, "scattering table (8 classes)\n"
                                  "turbulence table (8 rows): L0–L5 + 2 ablation")
    inner(ax, 435, 415, 685, 470, "5 pseudo-random streams\n"
                                  "2 per fading draw · 3 per attempt")
    inner(ax, 435, 340, 685, 405, "click = signal ∨ background\n"
                                  "error = polarization / fair coin")
    inner(ax, 435, 265, 685, 330, "intensity class\nsignal · weak decoy · vacuum")

    box(ax, 745, 565, 870, 650, "OOK RX", "deserializer", fc=B_DBLUE,
        tc="white", fs_title=13, fs_sub=9)
    box(ax, 745, 450, 870, 530, "BOB", "BB84 decoder", fc=B_GREEN, fs_title=13)
    box(ax, 745, 325, 870, 415, "Sifting & Error", "basis match, error",
        fc=B_GREY, fs_title=12, fs_sub=9)
    box(ax, 745, 200, 870, 280, "MUX", "real-time mode", fc=B_RED,
        fs_title=13, fs_sub=9)

    # ================= zone 3 =================
    box(ax, 965, 545, 1185, 650, "Channel Monitor",
        "QBER, click rate, jitter", fc=B_ORANGE, fs_title=13, fs_sub=9)
    box(ax, 965, 405, 1185, 510, "Adaptive Controller",
        "mode · µ · λ · slot", fc=B_MBLUE, tc="white", fs_title=13, fs_sub=9)

    # ================= zone 4 =================
    box(ax, 965, 230, 1170, 340, "Security-block monitor",
        "per 2¹⁶ attempts: attempts,\ntrials, clicks, errors — by basis\n"
        "and intensity class · block fading\nno QBER formed · 1825 LE",
        fc=B_AMBER, fs_title=11.5, fs_sub=8.4)
    box(ax, 965, 120, 1170, 215, "Throughput meter",
        "measured attempt rate\n110 LE", fc=B_AMBER, fs_title=11.5,
        fs_sub=8.4)
    box(ax, 1205, 120, 1465, 340, "Binary block reporter",
        "packet v6 · 144 byte\nXOR checksum · watchdog\n3.13 ms = 48 % of a block\n1560 LE",
        fc=B_AMBER, fs_title=11.5, fs_sub=8.4)

    # ================= notes box =================
    ax.text(1231, 648, "Altera Cyclone II\nEP2C20F484C7 @ 50 MHz",
            fontsize=9.6, fontweight="bold", color=INK, ha="left", va="top",
            linespacing=1.4)
    ax.text(1231, 598, "11 385 LE (61 %) · 4 120 reg\n"
                       "setup slack +2.855 ns\n"
                       "block-acquisition path\n4 102 LE / 2 524 reg",
            fontsize=8.8, color="#333333", ha="left", va="top", linespacing=1.5)
    keys = [(A_DATA, "per-attempt datapath"),
            (A_FAST, "real-time path (10 MHz)"),
            (A_BLK, "block statistics & report"),
            (A_CTRL, "adaptive control"),
            (A_CFG, "configuration / run-time"),
            (A_RAND, "randomness")]
    y = 512
    for c, t in keys:
        ax.add_line(Line2D([1234, 1262], [y, y], color=c, lw=2.4, zorder=6))
        ax.text(1270, y, t, fontsize=8.8, color="#333333", ha="left",
                va="center", zorder=6)
        y -= 22

    # ================= arrows =================
    # Labels are readable phrases, as in the original export; the RTL wire
    # name behind each one is tabulated in docs/FIG1_REDRAW_GUIDE.md.
    # randomness
    arrow(ax, [(235, 600), (280, 600)], A_RAND, "Alice's data\n& basis",
          (257, 608), lva="bottom")
    arrow(ax, [(235, 565), (255, 565), (255, 165), (725, 165), (725, 370),
               (745, 370)], A_RAND, "Bob's measuring basis", (480, 172),
          lva="bottom")

    # OOK datapath
    arrow(ax, [(335, 560), (335, 520)], A_DATA, "Encoded qubit", (342, 540),
          lha="left", lva="center")
    arrow(ax, [(390, 480), (420, 480)], A_DATA, "Optical pulse", (405, 480),
          lha="center", lva="center", rot=90)
    arrow(ax, [(700, 605), (745, 605)], A_DATA, "Received pulse", (722, 605),
          lha="center", lva="center", rot=90)
    arrow(ax, [(807, 565), (807, 530)], A_DATA, "Detected qubit", (814, 548),
          lha="left", lva="center")
    # BOB -> Sifting: adjacent boxes, a label would only crowd the gap that
    # the adaptive-control return path also has to cross.
    arrow(ax, [(807, 450), (807, 415)], A_DATA)

    # real-time path
    arrow(ax, [(135, 375), (135, 310)], A_FAST, "Paced attempt", (128, 343),
          lha="right", lva="center")
    arrow(ax, [(160, 80), (160, 230)], A_FAST, "Host-paced\nattempt",
          (167, 150), lha="left", lva="center")
    arrow(ax, [(200, 270), (420, 270)], A_FAST, "Attempt strobe", (310, 277),
          lva="bottom")
    arrow(ax, [(235, 430), (420, 430)], A_FAST, "Intensity class", (327, 437),
          lva="bottom")
    arrow(ax, [(420, 355), (175, 355), (175, 375)], A_FAST,
          "Click / loss / error", (300, 347), lva="top", dashed=True)
    arrow(ax, [(50, 375), (50, 95), (807, 95), (807, 200)], A_FAST,
          "Real-time qubit events", (430, 102), lva="bottom")

    # OOK events into the event mux (adjacent boxes -- no label needed)
    arrow(ax, [(807, 325), (807, 280)], A_EVT)

    # block statistics
    arrow(ax, [(870, 265), (930, 265), (930, 620), (965, 620)], A_BLK,
          "Qubit events", (922, 455), lha="center", lva="center", rot=90)
    arrow(ax, [(870, 250), (965, 250)], A_BLK, "Qubit events", (917, 258),
          lva="bottom")
    arrow(ax, [(845, 200), (845, 168), (965, 168)], A_BLK, "Attempt count",
          (905, 176), lva="bottom")
    arrow(ax, [(700, 300), (965, 300)], A_BLK, "Fading word,\nattempt index",
          (878, 308), lva="bottom")
    arrow(ax, [(1170, 285), (1205, 285)], A_BLK, "Block\ncounters", (1187, 293),
          lva="bottom")
    arrow(ax, [(1170, 170), (1205, 170)], A_BLK, "Rate", (1187, 178),
          lva="bottom")
    arrow(ax, [(1335, 120), (1335, 80)], A_BLK, "Block packet\n(144 byte)",
          (1343, 100), lha="left", lva="center")

    # adaptive control
    arrow(ax, [(1075, 545), (1075, 510)], A_CTRL, "Channel quality\nreport",
          (1082, 527), lha="left", lva="center")
    arrow(ax, [(965, 435), (700, 435)], A_CTRL, "Transmitter settings",
          (830, 443), lva="bottom")
    arrow(ax, [(700, 547), (910, 547), (910, 585), (965, 585)], A_CTRL,
          "Expected click rate", (755, 555), lva="bottom")

    # configuration
    arrow(ax, [(520, 80), (520, 195)], A_CFG,
          "Water, range, λ,\nturbulence, background", (527, 138), lha="left",
          lva="center")
    arrow(ax, [(218, 80), (218, 375)], A_CFG,
          "Real-time rate,\ndecoy setup", (225, 210), lha="left", lva="center")

    for d in OUTS:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / "QKD_UWOC_Block.pdf")
        fig.savefig(d / "QKD_UWOC_Block.png", dpi=DPI)
        print("wrote %s" % (d / "QKD_UWOC_Block.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
