#!/usr/bin/env python3
"""
fig1_block_diagram_pptx.py -- Fig. 1 as an editable PowerPoint deck.

Ports python/fig1_block_diagram.py to native PowerPoint shapes so the figure can
be edited by hand: every zone, block and arrow is a real shape, every caption is
a real text box.  The geometry is the same "display unit" system as the
matplotlib source (100 units = 1 inch), so the slide is 15 x 7.2 in and the font
sizes below are the same point sizes as the PDF export.

    python python/fig1_block_diagram_pptx.py
"""

from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN, MSO_AUTO_SIZE
from pptx.oxml.ns import qn
from pptx.util import Emu, Pt

ROOT = Path(__file__).resolve().parent.parent
OUTS = [ROOT / "Images"]

W, H = 1500.0, 720.0          # display units, 100 units = 1 inch
EMU_PER_UNIT = 914400 / 100.0  # 1 unit = 0.01 in
PT_PER_UNIT = 72.0 / 100.0    # 1 unit = 0.72 pt
FONT = "Arial"

# ---- palette, matched to the original export -------------------------------
Z_GREEN = ("E8F5E9", "66BB6A")
Z_BLUE = ("E3F2FD", "64B5F6")
Z_PEACH = ("FFF3E0", "FFB74D")
Z_AMBER = ("FFF8E1", "FFA726")
Z_PINK = ("F9E3F5", "CE93D8")
Z_GREY = ("F4F4F2", "BDBDBD")

B_GREEN = "A5D6A7"
B_TEAL = "4DB6AC"
B_LBLUE = "90CAF9"
B_MBLUE = "42A5F5"
B_DBLUE = "1565C0"
B_ORANGE = "EF9A6A"
B_AMBER = "FFCC80"
B_GREY = "D7D7D7"
B_RED = "EF9A9A"
B_INNER = "CFE6FB"

A_RAND = "43A047"
A_DATA = "0D47A1"
A_FAST = "00897B"
A_EVT = "546E7A"
A_BLK = "EF6C00"
A_CTRL = "B25000"
A_CFG = "8E24AA"

INK = "111111"
WHITE = "FFFFFF"

FS_ZONE = 15
FS_TITLE = 14
FS_SUB = 10
FS_INNER = 9.5
FS_ARROW = 9


# ---- unit helpers ----------------------------------------------------------
def ex(x):
    """display x -> EMU"""
    return Emu(int(round(x * EMU_PER_UNIT)))


def ey(y):
    """display y (bottom-up) -> EMU (top-down)"""
    return Emu(int(round((H - y) * EMU_PER_UNIT)))


def rgb(h):
    return RGBColor.from_string(h)


def _ln(shape):
    """the <a:ln> element of a shape, created if missing"""
    spPr = shape._element.spPr
    ln = spPr.find(qn("a:ln"))
    if ln is None:
        ln = spPr.get_or_add_ln()
    return ln


def arrowhead(shape, w="med", length="med"):
    ln = _ln(shape)
    for tag in ("a:tailEnd",):
        el = ln.find(qn(tag))
        if el is not None:
            ln.remove(el)
    el = ln.makeelement(qn("a:tailEnd"), {"type": "triangle", "w": w, "len": length})
    ln.append(el)


# ---- primitives ------------------------------------------------------------
def _textframe(shape, anchor=MSO_ANCHOR.MIDDLE, wrap=False):
    tf = shape.text_frame
    tf.word_wrap = wrap
    tf.auto_size = MSO_AUTO_SIZE.NONE
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    return tf


def _fill_lines(tf, lines, size, bold, color, align=PP_ALIGN.CENTER,
                spacing=1.15, space_before=0):
    """lines: list of (text, size, bold, color) or plain strings"""
    first = True
    for item in lines:
        text, sz, bd, col, sb = item
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.alignment = align
        p.line_spacing = spacing
        if sb:
            p.space_before = Pt(sb)
        r = p.add_run()
        r.text = text
        r.font.size = Pt(sz)
        r.font.bold = bd
        r.font.name = FONT
        r.font.color.rgb = rgb(col)


def zone(slide, x0, y0, x1, y1, title, colors, title_dx=0.0):
    fc, ec = colors
    w, h = x1 - x0, y1 - y0
    sh = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                ex(x0), ey(y1), ex(w), ex(h))
    sh.adjustments[0] = 14.0 / min(w, h)
    sh.fill.solid()
    sh.fill.fore_color.rgb = rgb(fc)
    sh.line.color.rgb = rgb(ec)
    sh.line.width = Pt(1.8)
    sh.shadow.inherit = False
    sh.text_frame.text = ""
    if title:
        label(slide, x0 + 6 + title_dx, y1 + 9, title, FS_ZONE, INK,
              ha="left", va="bottom", bold=True)
    return sh


def box(slide, x0, y0, x1, y1, title, sub=None, fc=B_LBLUE, tc=INK,
        fs_title=FS_TITLE, fs_sub=FS_SUB, round_size=8):
    w, h = x1 - x0, y1 - y0
    sh = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                ex(x0), ey(y1), ex(w), ex(h))
    sh.adjustments[0] = float(round_size) / min(w, h)
    sh.fill.solid()
    sh.fill.fore_color.rgb = rgb(fc)
    sh.line.fill.background()
    sh.shadow.inherit = False
    tf = _textframe(sh)
    lines = [(t, fs_title, True, tc, 0) for t in title.split("\n")] if title else []
    if sub:
        for i, t in enumerate(sub.split("\n")):
            lines.append((t, fs_sub, False, tc, 3 if i == 0 else 0))
    if lines:
        _fill_lines(tf, lines, fs_title, True, tc, spacing=1.2)
    return sh


def inner(slide, x0, y0, x1, y1, text):
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE,
                                ex(x0), ey(y1), ex(x1 - x0), ex(y1 - y0))
    sh.fill.solid()
    sh.fill.fore_color.rgb = rgb(B_INNER)
    sh.line.color.rgb = rgb("7FB8E8")
    sh.line.width = Pt(0.9)
    sh.shadow.inherit = False
    tf = _textframe(sh)
    _fill_lines(tf, [(t, FS_INNER, False, INK, 0) for t in text.split("\n")],
                FS_INNER, False, INK, spacing=1.2)
    return sh


def label(slide, lx, ly, text, size, color, ha="center", va="bottom",
          bold=False, rot=0, spacing=1.2):
    """Free-standing text, positioned the way matplotlib's ha/va would."""
    lines = text.split("\n")
    n = len(lines)
    ncols = max(len(s) for s in lines)
    w = max(ncols * size * 0.56, 6.0) / PT_PER_UNIT   # display units
    h = (n * size * spacing + 1.0) / PT_PER_UNIT

    if ha == "left":
        x = lx
        align = PP_ALIGN.LEFT
    elif ha == "right":
        x = lx - w
        align = PP_ALIGN.RIGHT
    else:
        x = lx - w / 2.0
        align = PP_ALIGN.CENTER

    if va == "bottom":
        top = ly + h
    elif va == "top":
        top = ly
    else:
        top = ly + h / 2.0

    tb = slide.shapes.add_textbox(ex(x), ey(top), ex(w), ex(h))
    tf = _textframe(tb)
    _fill_lines(tf, [(t, size, bold, color, 0) for t in lines],
                size, bold, color, align=align, spacing=spacing)
    if rot:
        tb.rotation = (360 - rot) % 360
    return tb


def arrow(slide, pts, color, label_text=None, lpos=None, lha="center",
          lva="bottom", dashed=False, lw=2.0, rot=0):
    """Orthogonal polyline with an arrowhead on the last segment."""
    if len(pts) == 2:
        (x0, y0), (x1, y1) = pts
        sh = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT,
                                        ex(x0), ey(y0), ex(x1), ey(y1))
    else:
        b = slide.shapes.build_freeform(ex(pts[0][0]), ey(pts[0][1]), scale=1.0)
        b.add_line_segments([(ex(x), ey(y)) for x, y in pts[1:]], close=False)
        sh = b.convert_to_shape()
        sh.fill.background()
    sh.line.color.rgb = rgb(color)
    sh.line.width = Pt(lw)
    sh.shadow.inherit = False
    if dashed:
        ln = _ln(sh)
        for el in ln.findall(qn("a:prstDash")):
            ln.remove(el)
        ln.append(ln.makeelement(qn("a:prstDash"), {"val": "dash"}))
    arrowhead(sh)
    if label_text:
        label(slide, lpos[0], lpos[1], label_text, FS_ARROW, color,
              ha=lha, va=lva, rot=rot, spacing=1.25)
    return sh


def main():
    prs = Presentation()
    prs.slide_width = ex(W)
    prs.slide_height = ex(H)
    slide = prs.slides.add_slide(prs.slide_layouts[6])   # blank

    # ================= zones =================
    zone(slide, 20, 105, 250, 665, "ATTEMPT SOURCES", Z_GREEN)
    zone(slide, 265, 105, 870, 665, "QKD DATAPATH", Z_BLUE)
    zone(slide, 950, 380, 1200, 665, "MONITOR & ADAPT", Z_PEACH)
    zone(slide, 950, 105, 1480, 355, "REAL-TIME BLOCK ACQUISITION", Z_AMBER)
    zone(slide, 1215, 380, 1480, 665, "", Z_GREY)
    zone(slide, 20, 18, 1480, 80, "", Z_PINK)

    # ================= host strip =================
    label(slide, 750, 58, "BB84 CONTROLLER  &  HOST INTERFACE", 16, INK,
          ha="center", va="center", bold=True)
    label(slide, 750, 33, "UART RX (commands)   ·   UART TX arbiter, "
                          "block reporter has priority   ·   460.8 kbaud",
          9.5, "4A3050", ha="center", va="center")

    # ================= zone 1 =================
    box(slide, 35, 545, 235, 650, "TRNG ×3",
        "ring oscillator\nOOK path only", fc=B_GREEN)
    box(slide, 35, 375, 235, 500, "REAL-TIME\nATTEMPT ENGINE",
        "10 MHz pacer · xorshift64 ×3\n607 LE / 322 reg", fc=B_TEAL,
        tc=WHITE, fs_title=13)
    box(slide, 70, 230, 200, 310, "MUX", "real-time mode", fc=B_RED,
        fs_title=13, fs_sub=9)

    # ================= zone 2 =================
    box(slide, 280, 560, 390, 650, "ALICE", "BB84 encoder", fc=B_GREEN,
        fs_title=13)
    box(slide, 280, 440, 390, 520, "OOK TX", "4-slot serializer", fc=B_MBLUE,
        tc=WHITE, fs_title=13, fs_sub=9)

    box(slide, 420, 195, 700, 650, "", fc=B_LBLUE, round_size=10)
    label(slide, 560, 625, "UWOC Channel emulator", 14, INK,
          ha="center", va="center", bold=True)
    inner(slide, 435, 545, 685, 597, "path-loss & polarization tables")
    inner(slide, 435, 480, 685, 535, "scattering table (8 classes)\n"
                                     "turbulence table (8 rows): L0–L5 + 2 ablation")
    inner(slide, 435, 415, 685, 470, "5 pseudo-random streams\n"
                                     "2 per fading draw · 3 per attempt")
    inner(slide, 435, 340, 685, 405, "click = signal ∨ background\n"
                                     "error = polarization / fair coin")
    inner(slide, 435, 265, 685, 330, "intensity class\n"
                                     "signal · weak decoy · vacuum")

    box(slide, 745, 565, 870, 650, "OOK RX", "deserializer", fc=B_DBLUE,
        tc=WHITE, fs_title=13, fs_sub=9)
    box(slide, 745, 450, 870, 530, "BOB", "BB84 decoder", fc=B_GREEN,
        fs_title=13)
    box(slide, 745, 325, 870, 415, "Sifting & Error", "basis match, error",
        fc=B_GREY, fs_title=12, fs_sub=9)
    box(slide, 745, 200, 870, 280, "MUX", "real-time mode", fc=B_RED,
        fs_title=13, fs_sub=9)

    # ================= zone 3 =================
    box(slide, 965, 545, 1185, 650, "Channel Monitor",
        "QBER, click rate, jitter", fc=B_ORANGE, fs_title=13, fs_sub=9)
    box(slide, 965, 405, 1185, 510, "Adaptive Controller",
        "mode · µ · λ · slot", fc=B_MBLUE, tc=WHITE,
        fs_title=13, fs_sub=9)

    # ================= zone 4 =================
    box(slide, 965, 230, 1170, 340, "Security-block monitor",
        "per 2¹⁶ attempts: attempts,\ntrials, clicks, errors — by basis\n"
        "and intensity class · block fading\nno QBER formed · 1825 LE",
        fc=B_AMBER, fs_title=11.5, fs_sub=8.4)
    box(slide, 965, 120, 1170, 215, "Throughput meter",
        "measured attempt rate\n110 LE", fc=B_AMBER, fs_title=11.5, fs_sub=8.4)
    box(slide, 1205, 120, 1465, 340, "Binary block reporter",
        "packet v6 · 144 byte\nXOR checksum · watchdog\n"
        "3.13 ms = 48 % of a block\n1560 LE",
        fc=B_AMBER, fs_title=11.5, fs_sub=8.4)

    # ================= notes box =================
    label(slide, 1231, 648, "Altera Cyclone II\nEP2C20F484C7 @ 50 MHz", 9.6,
          INK, ha="left", va="top", bold=True, spacing=1.4)
    label(slide, 1231, 598, "11 385 LE (61 %) · 4 120 reg\n"
                            "setup slack +2.855 ns\n"
                            "block-acquisition path\n4 102 LE / 2 524 reg",
          8.8, "333333", ha="left", va="top", spacing=1.5)
    keys = [(A_DATA, "per-attempt datapath"),
            (A_FAST, "real-time path (10 MHz)"),
            (A_BLK, "block statistics & report"),
            (A_CTRL, "adaptive control"),
            (A_CFG, "configuration / run-time"),
            (A_RAND, "randomness")]
    y = 512
    for c, t in keys:
        ln = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT,
                                        ex(1234), ey(y), ex(1262), ey(y))
        ln.line.color.rgb = rgb(c)
        ln.line.width = Pt(2.4)
        ln.shadow.inherit = False
        label(slide, 1270, y, t, 8.8, "333333", ha="left", va="center")
        y -= 22

    # ================= arrows =================
    # randomness
    arrow(slide, [(235, 600), (280, 600)], A_RAND, "Alice's data\n& basis",
          (257, 608), lva="bottom")
    arrow(slide, [(235, 565), (255, 565), (255, 165), (725, 165), (725, 370),
                  (745, 370)], A_RAND, "Bob's measuring basis", (480, 172),
          lva="bottom")

    # OOK datapath
    arrow(slide, [(335, 560), (335, 520)], A_DATA, "Encoded qubit", (342, 540),
          lha="left", lva="center")
    arrow(slide, [(390, 480), (420, 480)], A_DATA, "Optical pulse", (405, 480),
          lha="center", lva="center", rot=90)
    arrow(slide, [(700, 605), (745, 605)], A_DATA, "Received pulse", (722, 605),
          lha="center", lva="center", rot=90)
    arrow(slide, [(807, 565), (807, 530)], A_DATA, "Detected qubit", (814, 548),
          lha="left", lva="center")
    arrow(slide, [(807, 450), (807, 415)], A_DATA)

    # real-time path
    arrow(slide, [(135, 375), (135, 310)], A_FAST, "Paced attempt", (128, 343),
          lha="right", lva="center")
    arrow(slide, [(160, 80), (160, 230)], A_FAST, "Host-paced\nattempt",
          (167, 150), lha="left", lva="center")
    arrow(slide, [(200, 270), (420, 270)], A_FAST, "Attempt strobe", (310, 277),
          lva="bottom")
    arrow(slide, [(235, 430), (420, 430)], A_FAST, "Intensity class", (327, 437),
          lva="bottom")
    arrow(slide, [(420, 355), (175, 355), (175, 375)], A_FAST,
          "Click / loss / error", (300, 347), lva="top", dashed=True)
    arrow(slide, [(50, 375), (50, 95), (807, 95), (807, 200)], A_FAST,
          "Real-time qubit events", (430, 102), lva="bottom")

    # OOK events into the event mux
    arrow(slide, [(807, 325), (807, 280)], A_EVT)

    # block statistics
    arrow(slide, [(870, 265), (930, 265), (930, 620), (965, 620)], A_BLK,
          "Qubit events", (922, 455), lha="center", lva="center", rot=90)
    arrow(slide, [(870, 250), (965, 250)], A_BLK, "Qubit events", (917, 258),
          lva="bottom")
    arrow(slide, [(845, 200), (845, 168), (965, 168)], A_BLK, "Attempt count",
          (905, 176), lva="bottom")
    arrow(slide, [(700, 300), (965, 300)], A_BLK, "Fading word,\nattempt index",
          (878, 308), lva="bottom")
    arrow(slide, [(1170, 285), (1205, 285)], A_BLK, "Block\ncounters", (1187, 293),
          lva="bottom")
    arrow(slide, [(1170, 170), (1205, 170)], A_BLK, "Rate", (1187, 178),
          lva="bottom")
    arrow(slide, [(1335, 120), (1335, 80)], A_BLK, "Block packet\n(144 byte)",
          (1343, 100), lha="left", lva="center")

    # adaptive control
    arrow(slide, [(1075, 545), (1075, 510)], A_CTRL, "Channel quality\nreport",
          (1082, 527), lha="left", lva="center")
    arrow(slide, [(965, 435), (700, 435)], A_CTRL, "Transmitter settings",
          (830, 443), lva="bottom")
    arrow(slide, [(700, 547), (910, 547), (910, 585), (965, 585)], A_CTRL,
          "Expected click rate", (755, 555), lva="bottom")

    # configuration
    arrow(slide, [(520, 80), (520, 195)], A_CFG,
          "Water, range, λ,\nturbulence, background", (527, 138),
          lha="left", lva="center")
    arrow(slide, [(218, 80), (218, 375)], A_CFG,
          "Real-time rate,\ndecoy setup", (225, 210), lha="left", lva="center")

    for d in OUTS:
        d.mkdir(parents=True, exist_ok=True)
        out = d / "QKD_UWOC_Block.pptx"
        prs.save(str(out))
        print("wrote %s" % out)


if __name__ == "__main__":
    main()
