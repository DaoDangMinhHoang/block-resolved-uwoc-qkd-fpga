#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fpga_collect.py — Collects BB84/UWOC data on the FPGA for the paper
═══════════════════════════════════════════════════════════════════════════════

DIFFERENT from bb84_uwoc_measure.py: that script sweeps fast to LOOK; this one
runs LONG to PRODUCE PUBLISHABLE DATA, so it has three things the other lacks:

  [1] A DUAL BUDGET per measurement point — stop on ENOUGH SAMPLES (n_sift ≥
      target) or on TIME OUT (cap), whichever comes first. Underwater P_click falls
      ~2 decades across the range of interest, so a fixed qubit count would either
      waste hours on the near points or never finish the far ones.
  [2] CHECKPOINTING — every finished point writes a CSV row immediately. A re-run
      skips points already present. A power cut at 3am does not lose the whole session.
  [3] PER-CLICK LOGGING — irrad (SNR proxy) + bmatch + err for every click, so the
      sliding-window / vs-SNR QBER figures (like Fig. 3 of the previous paper) can be built.

Run:
    python python/fpga_collect.py --port COM10 --phase fixed
    python python/fpga_collect.py --port COM10 --phase adaptive   # after flipping SW[1]=1

⚠ --chunk is always 1. The bitstream currently on the board has NO command FIFO;
  sending back-to-back loses commands and drives P_click down ~0.4× (see tb_cmd_fifo.v).
"""

from __future__ import annotations

import argparse
import csv
import glob
import io
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np                                          # noqa: E402

from bb84_uwoc_measure import (                             # noqa: E402
    CONF, MIN_SIFT, QBER_LIMIT, WATER_ORDER, _parse_full, qber_upper,
)
from uwoc_channel_model import (                            # noqa: E402
    LAMBDA_OPTIONS, WATER_TYPES, LinkConfig, UWOCChannel, h2,
)

import serial                                               # noqa: E402

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "data")
POINTS_CSV = os.path.join(DATA_DIR, "fpga_points.csv")

# Send commands this much slower than the FPGA's mean consumption rate. At a
# utilisation of exactly 1 the command FIFO is critically loaded and its
# occupancy random-walks into overflow; 15% of headroom keeps it bounded at a
# cost of 15% throughput.
PACE_MARGIN = 1.15

POINT_FIELDS = [
    "tag", "phase", "water", "dist_idx", "distance_m", "turb", "lam_nm",
    "n_qubit", "n_click", "n_sift", "n_err", "p_click", "sift_eff",
    "qber", "qber_hi", "secure", "skr_per_pulse", "seconds", "rate_qps",
    "p_click_model", "qber_model", "coh_sel", "coh_qubits", "dyn_walk",
    "timestamp",
]

# Per-click log. `attempt` is the [v12] 24-bit attempt index straight from the
# FPGA — the only thing that lets the analysis put a click in the right coherence
# block. Older logs have the column empty; the window analysis then skips them.
CLICK_FIELDS = ["attempt", "a_data", "a_basis", "b_basis", "bob", "bmatch",
                "err", "irrad", "mode", "mu", "lam_idx"]

# Adaptive controller modes, as reported in the [v13] per-qubit line.
MODE_NAMES = {0: "AGGRESSIVE", 1: "MODERATE", 2: "CONSERVATIVE", 3: "PAUSE"}

# coh_sel = k on the wire ⇒ the emulator freezes h for 2^(k+5) qubit events.
COH_BASE_LOG2 = 5


def coh_qubits(coh_sel: int) -> int:
    """Pulses per coherence block for a given selector (0 = legacy clock timer)."""
    return 0 if coh_sel == 0 else 1 << (coh_sel + COH_BASE_LOG2)


# ═══════════════════════════════════════════════════════════════════════════
# Link to the FPGA
# ═══════════════════════════════════════════════════════════════════════════
class Link:
    def __init__(self, port: str, baud: int = 460800, timeout: float = 0.004,
                 qubit_us: float = 220.0, report_us: float = 4000.0):
        # 4 ms: calibrated on COM10 — 86.6 qubits/s, and the P_click obtained is
        # 104% of the model value, so this timeout does NOT cut clicks off.
        self.timeout = timeout
        # Cost of a NO-CLICK qubit: TX 50 µs + rx_timeout 160 µs + FSM_GAP 10 µs
        # @ TURBO_SLOT = 500.
        self.qubit_s = qubit_us * 1e-6
        # EXTRA cost of a qubit that clicks: the FSM parks in FSM_REPORT for
        # 200_000 cycles @50 MHz (top_module.v:708) while the 34-byte report line
        # goes out — 34 B @115200 is 2.95 ms, so the 4 ms is real, not slack.
        # Mean cost per qubit is therefore qubit_s + P_click·report_s, and pacing
        # on qubit_s alone over-pumps the command FIFO on every clicking link.
        self.report_s = report_us * 1e-6
        self.ser = serial.Serial(port, baud, timeout=timeout)
        self._buf = b""
        time.sleep(0.5)
        self.ser.reset_input_buffer()

    def configure(self, water: int, dist: int, turb: int, lam: int,
                  coh: int = 0, dyn_walk: bool = False):
        """Load the channel configuration, then reset the statistics.

        The 0x01 reset MUST go last: on the [v12] bitstream it also restarts the
        coherence-block grid inside uwoc_channel, so that the `total` field of every
        subsequent report line and the fading block index share one origin. Sending
        it before the configuration bytes would leave the two off by however many
        qubits ran in between.
        """
        for b in (0x30 | (dist & 0x0F),
                  0x40 | ((water & 0x3) << 2) | (lam & 0x3),
                  0x50 | (0x08 if dyn_walk else 0x00) | (turb & 0x07),
                  0x60 | (coh & 0x0F),
                  0x01):
            self.ser.write(bytes([b]))
            time.sleep(0.02)
        time.sleep(0.05)
        self._buf = b""
        self.ser.reset_input_buffer()

    def _harvest(self):
        """Return every COMPLETE report line waiting in the OS buffer, non-blocking.

        Why not readline() per qubit: underwater >99% of qubits do not click, so a
        per-qubit readline() pays the full serial timeout on almost every qubit and
        the measurement rate collapses to 1/timeout regardless of how fast the FPGA
        runs. Reading the buffer in bulk decouples the two — the FPGA's own
        ~220 µs/qubit becomes the limit instead.
        """
        nb = self.ser.in_waiting
        if nb:
            self._buf += self.ser.read(nb)
        if b"\n" not in self._buf:
            return []
        parts = self._buf.split(b"\n")
        self._buf = parts[-1]                 # keep the partial tail for next time
        return [p.decode("ascii", "ignore") for p in parts[:-1]]

    def run_point(self, rng, target_sift: int, cap_s: float, p_basis_z=0.5,
                  progress_every=20.0, chunk=1, target_qubit: int = 0):
        """Run until target_sift sifted bits, target_qubit attempts, or cap_s seconds.

        ⚠ WHICH BUDGET TO USE. Stopping on `target_sift` is a rule CORRELATED WITH
        THE FADING: a point that happens to open on a bright stretch reaches its
        sifted quota sooner and ends there, so the sample over-represents high-h
        blocks. The pooled P_click is biased up by roughly σ²_h / n_blocks, and the
        dispersion statistics are biased DOWN because the run stops right after a
        run of good blocks. Measured: the L5 point at 25 m finished in 129 blocks
        where L1 needed 154, and its P_click came out 1.21× the model.

        For anything whose measurand is a DISPERSION (section B), pass
        `target_qubit` instead: a fixed number of attempts is independent of the
        fading, so the block sample is unbiased. The adaptive sifted-bit budget
        stays right for the range sweeps, where P_click spans two decades and a
        fixed attempt count would either waste hours or collect nothing.

        `chunk` sends that many qubit commands back-to-back, then paces the loop at
        the FPGA's own speed (chunk × qubit_s) while harvesting click lines in bulk.
        Only valid on a bitstream with the command FIFO (top_module.v [v11],
        CMD_DEPTH = 64) — on older RTL the FPGA silently drops commands and P_click
        comes out ~0.4× the model. Keep chunk ≤ 64 and verify with --chunk-check
        before committing to a long run.

        A click line can arrive after its own chunk's window closed; it is simply
        harvested by a later pass. Only aggregate counts are used, so the shifted
        arrival changes nothing — every line is still counted exactly once.
        """
        n = n_click = n_sift = n_err = 0
        clicks = []
        t0 = time.perf_counter()
        t_next = progress_every
        # Clamped to half the 64-deep command FIFO. Pacing at exactly the mean
        # consumption rate leaves the queue critically loaded, so its length
        # random-walks; a chunk of 64 written into a 64-deep FIFO then overflows
        # on the first excursion. Measured against a FIFO model: chunk 64 dropped
        # 10.5% of commands, chunk 32 with PACE_MARGIN drops none.
        chunk = max(1, min(int(chunk), 32))

        def account(lines):
            nonlocal n_click, n_sift, n_err
            for ln in lines:
                r = _parse_full(ln)
                if r is None:
                    continue
                n_click += 1
                clicks.append(r)
                if r["bmatch"]:
                    n_sift += 1
                    n_err += r["err"]

        while True:
            el = time.perf_counter() - t0
            if n_sift >= target_sift or el >= cap_s:
                break
            if target_qubit and n >= target_qubit:
                break

            cmds = bytearray()
            for _ in range(chunk):
                ad = int(rng.integers(0, 2))
                ab = 0 if rng.random() < p_basis_z else 1
                bb = 0 if rng.random() < p_basis_z else 1
                cmds.append(0x80 | (ad << 2) | (ab << 1) | bb)
            self.ser.write(bytes(cmds))
            n += chunk

            # Pace the next write to the FPGA's ACTUAL consumption rate, which
            # depends on how often it clicks. p_hat uses Laplace smoothing so the
            # very first chunks — before any click has been seen — are paced as if
            # P_click were 1/2 and then relax downwards. Starting optimistic
            # instead would overflow the FIFO before the estimate had any data.
            p_hat = (n_click + 1) / (n + 2)
            deadline = time.perf_counter() + PACE_MARGIN * chunk * (
                self.qubit_s + p_hat * self.report_s)
            while True:
                remain = deadline - time.perf_counter()
                if remain <= 0:
                    break
                # Windows sleep granularity is ~1 ms; spin out the last stretch
                # rather than overshooting the pace on small chunks.
                if remain > 0.002:
                    time.sleep(0.001)
                account(self._harvest())
            account(self._harvest())

            if el >= t_next:
                t_next += progress_every
                print("\r      %6.1f s  n=%d  click=%d  sift=%d  err=%d  %.0f q/s   "
                      % (el, n, n_click, n_sift, n_err, n / max(el, 1e-9)),
                      end="", flush=True)

        # Drain the lines still in flight: the last chunk's clicks, plus any line
        # whose 34 bytes (~3 ms @115200) had not finished arriving.
        t_drain = time.perf_counter() + 0.5
        while time.perf_counter() < t_drain:
            time.sleep(0.02)
            account(self._harvest())

        el = time.perf_counter() - t0
        print("\r" + " " * 62 + "\r", end="")
        return dict(n_qubit=n, n_click=n_click, n_sift=n_sift, n_err=n_err,
                    seconds=el, clicks=clicks)

    def close(self):
        self.ser.close()


# ═══════════════════════════════════════════════════════════════════════════
# Check the switch state before spending hours on a run
# ═══════════════════════════════════════════════════════════════════════════
def check_switches(link: Link, rng, n=8000) -> float:
    """Return the ratio P_click(λ=450)/P_click(λ=650) in clear ocean, d = 5 m.

    cfg_lambda ONLY takes effect when adaptive is off (top_module.v:133 —
    active_lambda = adaptive_enable ? adapt_lambda : cfg_lambda). Hence:
        ratio ≈ 3.0  ⇒ SW[1] = 0, fixed mode, cfg_lambda is in effect.
        ratio ≈ 1.0  ⇒ adaptive is overriding λ, so "wavelength sweep" data is meaningless.
    The value 3.05 comes from psig_rom: 1.038e-2 (450 nm) vs 3.40e-3 (650 nm).
    """
    out = []
    for lam in (0, 2):
        # coh = 0 ON PURPOSE. This check compares two click counts of 8000 qubits
        # each; with a qubit-counted block of 65 536 the whole check would sit
        # inside ONE fading sample, so the ratio would carry the full spread of h
        # (CV ≈ 0.46 at L5) and could abort a perfectly good session. The clock
        # timer averages h away here, which is exactly what a ratio test wants.
        # The measurement loop re-sends the real coh_sel for every point.
        link.configure(0, 0, 5, lam, coh=0)
        c = 0
        for i in range(n):
            ad = int(rng.integers(0, 2))
            link.ser.write(bytes([0x80 | (ad << 2) | (i & 0x03)]))
            if _parse_full(link.ser.readline().decode("ascii", "ignore")):
                c += 1
        out.append(c)
        print("      λ = %d nm : %d click / %d = %.3e"
              % (LAMBDA_OPTIONS[lam], c, n, c / n))
    return (out[0] / max(out[1], 1))


# ═══════════════════════════════════════════════════════════════════════════
# Measurement matrix
# ═══════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════
# HAND-TUNED PLAN for the headline range sweep (section A)
# ═══════════════════════════════════════════════════════════════════════════
# (dist_idx, target_sift, cap_minutes) — clear ocean d_grid = 5,10,15,…,80 m.
#
# WHY THE TARGET FALLS WITH RANGE. QBER is a binomial estimate: its standard
# error is √(q(1−q)/n_sift), so the figure needs n_sift, not qubits. But
# sift_eff = ½·P_click drops ~2 decades over 5→35 m, so a flat target would put
# 90% of the session on the last point. These numbers spend the budget where it
# buys resolution and let the Clopper-Pearson bound speak for the far points.
#
# Model QBER (L1, 450 nm) and the σ each target buys:
#     5 m 1.64% ±0.40   10 m 2.18% ±0.46   15 m 2.68% ±0.57   20 m 3.17% ±0.72
#    25 m 3.69% ±0.94   30 m 4.40% ±1.68   35 m 5.84% ±2.5
# The steps between adjacent points (~0.5%) are comparable to σ out past 20 m —
# the TREND is resolved, adjacent points are not. Say so in the caption.
#
# The caps are a SAFETY NET at 2× the expected time, not a schedule — every point
# normally stops on target_sift well before its cap. They assume the ~3300-3900
# qubit/s that bulk harvesting plus FIFO-paced writes sustain; run --chunk-check
# to confirm the rate on your board, and pass --cap-scale if it comes out lower.
#
# sigma = √(q(1−q)/n_sift) at the model QBER, i.e. the resolution each target buys:
#    5 m  1.64% ±0.18   10 m  2.18% ±0.21   15 m  2.68% ±0.23   20 m  3.17% ±0.28
#   25 m  3.69% ±0.34   30 m  4.40% ±0.46   35 m  5.84% ±0.68   40 m  8.45% ±1.14
#   45 m 12.97% ±1.94   50 m 19.73% ±3.25
# The ~0.5%/point steps out to 30 m are now 2-3 sigma apart, so adjacent points
# separate rather than merely tracing a trend. The far points are deliberately
# coarse: they exist to bracket where QBER crosses the 11% Shor-Preskill bound
# (between 40 and 45 m in the model), not to be resolved individually.
RANGE_PLAN = [
    (0, 5000, 10),      #  5 m — ~5 min
    (1, 5000, 15),      # 10 m — ~8 min
    (2, 5000, 26),      # 15 m — ~13 min
    (3, 4000, 34),      # 20 m — ~17 min
    (4, 3000, 43),      # 25 m — ~22 min
    (5, 2000, 54),      # 30 m — ~27 min
    (6, 1200, 73),      # 35 m — ~37 min
    (7,  600, 78),      # 40 m — ~39 min
    (8,  300, 76),      # 45 m — ~38 min
    (9,  150, 66),      # 50 m — ~33 min
]

# The same budget logic for the other two water types (section D). The plans are
# NOT the clear-ocean one reused: coastal and harbor have their own d_grid and
# their P_click falls off a cliff far sooner, so the index at which the target
# must be cut differs. Each plan is sized for ~4.5 h — the same as one clear-ocean
# sweep — and brackets where the model QBER crosses the 11% bound:
#   coastal  crosses between idx 5 (13 m, 8.6%) and idx 6 (16 m, 14.7%)
#   harbor   crosses between idx 3 ( 2 m, 9.4%) and idx 4 (2.5 m, 17.5%)
# Beyond those indices the link is already dead (P_click at the Y₀ floor), so
# spending hours there buys nothing but a flat line at QBER = 50%.
RANGE_PLAN_COASTAL = [
    (0, 5000, 12),      #  2 m — ~6 min
    (1, 5000, 24),      #  4 m — ~12 min
    (2, 4000, 40),      #  6 m — ~20 min
    (3, 3000, 62),      #  8 m — ~31 min
    (4, 2000, 85),      # 10 m — ~42 min
    (5,  800, 75),      # 13 m — ~49 min
    (6,  300, 70),      # 16 m — ~48 min
    (7,  150, 60),      # 19 m — ~50 min
]

RANGE_PLAN_HARBOR = [
    (0, 5000, 20),      # 0.5 m — ~10 min
    (1, 4000, 52),      # 1.0 m — ~26 min
    (2, 1500, 66),      # 1.5 m — ~33 min
    (3,  600, 70),      # 2.0 m — ~44 min
    (4,  250, 70),      # 2.5 m — ~51 min
    (5,  120, 65),      # 3.0 m — ~51 min
    (6,   80, 60),      # 3.5 m — ~48 min
]

RANGE_PLANS = {
    "clear_ocean": RANGE_PLAN,
    "coastal": RANGE_PLAN_COASTAL,
    "harbor": RANGE_PLAN_HARBOR,
}


def build_matrix(phase: str, turb: int = 5, water_idx: int = 0,
                 sections: str = "ABCD", cap_scale: float = 1.0):
    """(tag, water, dist_idx, turb, lam, target_sift, cap_seconds, target_qubit)

    target_qubit = 0 means "budget by sifted bits"; a non-zero value stops the
    point at a fixed attempt count instead — see Link.run_point on why a
    dispersion measurement must not use a fading-correlated stopping rule.

    Far points never reach the target — that is deliberate: the Clopper-Pearson upper
    bound then says outright that the sample is small, instead of faking a QBER.

    `turb` applies to sections A and D (they must share a level to be comparable).
    B sweeps L1…L5 by definition; C stays at L3. Which level A and D run at barely
    matters for the MEAN: P_click is linear in h, so the expectation cancels the
    fading — verified, L1 and L5 give identical mean_metrics to every printed
    digit. Turbulence shows up in the per-click irrad distribution and in the
    between-window variance, which is what section B is for.
    """
    M = []

    # ---- A. QBER/SKR vs range — the headline figure --------------------
    # The turbulence level MUST be in the tag. Without it, a --turb 2 run
    # produces byte-identical tags to a --turb 1 run, so the checkpoint skips
    # every point and the session silently collects nothing. The level goes at
    # the END so tags still start with "A_dist", which is what paper_figs
    # matches on.
    if "A" in sections:
        for di, target, cap_min in RANGE_PLAN:
            M.append(("A_dist_d%d_L%d" % (di, turb), water_idx, di, turb, 0,
                      target, cap_min * 60, 0))

    # ---- B. QBER vs turbulence level, fixed range -----------------------
    # Expectation: the MEAN QBER stays nearly constant (P_click is linear in h,
    # so the expectation cancels the fading). Measuring that is itself a result —
    # it is what makes repeating the whole range sweep at L1…L5 unnecessary.
    #
    # WHERE the point sits matters, and it is NOT the near field. The turbulence
    # signature the hardware can actually resolve is the spread of the per-click
    # irrad field, i.e. CV(h) = √[(1+σ²_s)(1+σ²_ho) − 1]. Measured on the board,
    # CV bottoms out at ≈0.14 (12-bit h + 256-point inverse-CDF ROM), and at
    # d = 5 m the model puts L1 at 0.031 and L2 at 0.061 — both under that floor,
    # which is exactly why the collected d0 L1 and d0 L2 came out at 0.140 and
    # 0.141, indistinguishable. At dist_idx 4 (25 m clear ocean) the model gives
    # CV = 0.18 (L1) … 2.1 (L5), all above the floor and monotone in the level.
    #
    # [v12] With the coherence block counted in QUBITS (--coh 11 → 65 536), this
    # section no longer settles for the per-click irrad spread: each block is one
    # frozen fading sample, so the PER-BLOCK QBER becomes a real measurement and
    # its std / outage can be compared straight against
    #     UWOCChannel.window_statistics(window = 65536).
    #
    # Sizing is driven by the number of BLOCKS, and the budget is in ATTEMPTS, not
    # sifted bits. The first run of this section used a sifted-bit target and the
    # bias showed immediately: L1 needed 154 blocks to reach 6000 sifted while L5,
    # opening on bright fading, got there in 129 — so the strongest-turbulence point
    # had the SMALLEST and most favourably-selected sample, exactly backwards.
    #   dist_idx 4 (25 m), P_click ≈ 1.19e-3 → 78 clicks and 39 sifted per block,
    #   comfortably above channel_monitor's MIN_SIFT = 16.
    #   relative s.e. of a std estimate ≈ 1/√(2·N_blocks) → 160 blocks gives ±5.6%,
    #   against a model separation of std(L1) = 3.1% vs std(L3) = 4.1% vs
    #   std(L5) = 22.1%. 160 blocks ⇒ 1.05e7 attempts ⇒ ~47 min at 3700 qubit/s.
    if "B" in sections:
        for lv in range(1, 6):
            M.append(("B_turb_L%d" % lv, 0, 4, lv, 0,
                      10 ** 9, 70 * 60, 160 * 65536))

    # ---- C. Optimal wavelength depends on the water type ----------------
    # clear ocean favours blue; turbid harbour shifts to red (450:5.1e-3 < 650:7.8e-3).
    if "C" in sections:
        for li in range(3):
            M.append(("C_lam_clear_%d" % LAMBDA_OPTIONS[li], 0, 0, 3, li,
                      400, 720, 0))
        for li in range(3):
            M.append(("C_lam_harbor_%d" % LAMBDA_OPTIONS[li], 2, 0, 3, li,
                      400, 720, 0))

    # ---- D. Comparison of the three water types ------------------------
    # A FULL range sweep per water, not three sample points: the paper's water
    # comparison is a set of curves, and three points cannot show where each one
    # crosses the 11% bound. Uses the same turbulence level as section A so the
    # three curves are directly comparable — and the level goes in the tag for the
    # same reason it does in A (without it the checkpoint silently skips a re-run
    # at a different level).
    if "D" in sections:
        for wi, wname in ((1, "coastal"), (2, "harbor")):
            for di, target, cap_min in RANGE_PLANS[wname]:
                M.append(("D_%s_d%d_L%d" % (wname, di, turb), wi, di, turb, 0,
                          target, cap_min * 60, 0))

    return [(("%s_%s" % (phase, m[0]),) + tuple(m[1:6])
             + (m[6] * cap_scale, m[7]))
            for m in M]


# ═══════════════════════════════════════════════════════════════════════════
def model_expect(water: str, d: float, turb: int, lam_nm: int):
    cfg = LinkConfig(lam_nm=lam_nm)
    m = UWOCChannel(water, d, turb, cfg).mean_metrics()
    return m["p_click"], m["qber"]


def lambda_diagnosis(water: str, d: float, turb: int, lam_nm: int, p_meas: float):
    """Explain a low P_click as a WAVELENGTH problem, or return None. → list of lines

    In the adaptive phase λ is NOT the one the host asked for: the controller's
    hill-climber owns it, and while it probes a candidate the link runs on that
    candidate. In clear ocean 650 nm is 27× worse than 450 nm at 15 m and 92× at
    25 m, so a stuck probe drops P_click straight onto the Y₀ floor — the same
    symptom as a command FIFO overflow, which is why the 2026-08-09 session was
    spent chasing --chunk instead of the controller. The two are separable
    because dropped commands scale P_click by a factor INDEPENDENT of range,
    while a wrong λ scales it by exp(−Δc·d) and therefore lands exactly on the
    model of another wavelength.

    Two signatures, because a probe can be stuck for a whole point or only for
    part of it:
      STUCK — P_meas matches the model at some other λ (within 2×).
      FLAP  — P_meas sits between the best and the worst λ, i.e. the point is a
              time-average of a good λ and a dead probe. The mixing fraction is
              reported: on 2026-08-09 both 15 m and 20 m came out at 14%, the
              duty cycle of one 12-window hill-climb cycle.
    """
    pm = {L: model_expect(water, d, turb, L)[0] for L in LAMBDA_OPTIONS}
    p_ref = pm[lam_nm]
    fit = min(LAMBDA_OPTIONS,
              key=lambda L: abs(np.log(max(p_meas, 1e-300) / max(pm[L], 1e-300))))
    if fit != lam_nm and 0.5 <= p_meas / max(pm[fit], 1e-300) <= 2.0:
        return ["→ Số đo khớp mô hình ở λ = %d nm (%.3e), KHÔNG phải %d nm."
                % (fit, pm[fit], lam_nm),
                "  Bộ leo đồi bước sóng đang KẸT ở λ sai."]

    p_bad = min(pm.values())
    if p_bad < p_meas < 0.5 * p_ref:
        # p_meas = f·p_ref + (1−f)·p_bad ⇒ f = fraction of time on the good λ
        f = (p_meas - p_bad) / max(p_ref - p_bad, 1e-300)
        return ["→ Số đo nằm GIỮA mô hình %d nm (%.3e) và %d nm (%.3e):"
                % (lam_nm, p_ref, min(pm, key=pm.get), p_bad),
                "  bộ leo đồi đang NHẢY λ, chỉ ~%.0f%% thời gian ở %d nm."
                % (100 * f, lam_nm)]
    return None


def load_done():
    if not os.path.exists(POINTS_CSV):
        return set()
    with open(POINTS_CSV, newline="", encoding="utf-8") as f:
        return {r["tag"] for r in csv.DictReader(f)}


def _csv_text(fieldnames, rows, header: bool) -> str:
    """Render rows to CSV text in memory, so the on-disk write is a single step."""
    buf = io.StringIO(newline="")
    # extrasaction="ignore": a click dict carries the running sifted/error totals
    # from the report line as well, and those are derivable — no reason to let an
    # extra key abort a write that is holding hours of measurement.
    w = csv.DictWriter(buf, fieldnames=fieldnames, extrasaction="ignore")
    if header:
        w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def _write_retry(path, mode, text, wait_s: float = 300.0):
    """Write `text` to `path`, surviving a file lock instead of throwing away hours
    of measurement.

    On Windows, opening a CSV in Excel takes a share-mode lock, so the next append
    raises PermissionError and — before this retry existed — killed the whole
    session, discarding the point that had just finished collecting. Waiting is
    the right response: the lock vanishes the moment the file is closed, and a
    measurement point is worth far more than the delay.

    The RETRY WRAPS THE WHOLE WRITE, not just the open. Byte-range locks (the
    `msvcrt.locking` kind) let open() succeed and fail at write time instead, so
    guarding open() alone still crashes. Rendering the text up front makes the
    retry safe to repeat: a failed attempt writes nothing.

    Falls back to a sidecar file if the lock outlives `wait_s`, so data is never
    lost. Merge sidecars back with:  python python/fpga_collect.py --merge
    """
    t0 = time.perf_counter()
    warned = False
    while True:
        try:
            with open(path, mode, newline="", encoding="utf-8") as f:
                f.write(text)
            return path
        except PermissionError:
            el = time.perf_counter() - t0
            if el >= wait_s:
                break
            if not warned:
                warned = True
                print("\n  ⚠ %s đang bị KHÓA (thường là đang mở trong Excel)."
                      % os.path.basename(path))
                print("    Đóng file đó đi, phép đo sẽ tự ghi tiếp. Chờ tối đa "
                      "%.0f phút..." % (wait_s / 60))
            print("\r    ...chờ %.0f s" % el, end="", flush=True)
            time.sleep(2.0)

    # Lock outlived the patience budget — never drop the data on the floor.
    stem, ext = os.path.splitext(path)
    alt = "%s.part%d%s" % (stem, int(time.time()), ext)
    print("\n  ⚠ vẫn bị khóa sau %.0f phút — ghi tạm sang %s"
          % (wait_s / 60, os.path.basename(alt)))
    with open(alt, "w", newline="", encoding="utf-8") as f:
        f.write(text)
    return alt


def append_point(row):
    new = not os.path.exists(POINTS_CSV) or os.path.getsize(POINTS_CSV) == 0
    text = _csv_text(POINT_FIELDS, [row], header=new)
    out = _write_retry(POINTS_CSV, "a", text)
    if out != POINTS_CSV and not new:
        # The sidecar starts empty, so it needs a header the main file already has.
        with open(out, "w", newline="", encoding="utf-8") as f:
            f.write(_csv_text(POINT_FIELDS, [row], header=True))


def save_clicks(tag, clicks):
    p = os.path.join(DATA_DIR, "clicks_%s.csv" % tag)
    _write_retry(p, "w", _csv_text(CLICK_FIELDS, clicks, header=True))


def merge_sidecars():
    """Fold any fpga_points.part*.csv back into fpga_points.csv, skipping dup tags."""
    parts = sorted(glob.glob(os.path.join(
        DATA_DIR, "fpga_points.part*.csv")))
    if not parts:
        print("  Không có file .part nào cần gộp.")
        return
    done = load_done()
    added = 0
    for p in parts:
        with open(p, newline="", encoding="utf-8") as f:
            rows = [r for r in csv.DictReader(f) if r["tag"] not in done]
        for r in rows:
            append_point(r)
            done.add(r["tag"])
            added += 1
        print("  %s: gộp %d dòng" % (os.path.basename(p), len(rows)))
        os.rename(p, p + ".merged")
    print("  Đã thêm %d điểm vào %s" % (added, os.path.basename(POINTS_CSV)))


def main():
    ap = argparse.ArgumentParser(description="Thu số liệu BB84/UWOC trên FPGA")
    ap.add_argument("--port", help="cổng COM của FPGA; không cần khi --merge")
    ap.add_argument("--merge", action="store_true",
                    help="gộp các file fpga_points.part*.csv (sinh ra khi CSV "
                         "bị khóa lúc đang đo) trở lại fpga_points.csv rồi thoát")
    # [Q2] 460800, matching uart_rx/uart_tx in top_module.v. The old 115200
    # default silently decodes nothing on the current bitstream.
    ap.add_argument("--baud", type=int, default=460800)
    ap.add_argument("--timeout", type=float, default=0.004)
    ap.add_argument("--phase", default="fixed", choices=("fixed", "adaptive"))
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--skip-check", action="store_true",
                    help="bỏ qua bước kiểm tra SW[1] (~3 phút)")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="nhân trần thời gian và chia mục tiêu số mẫu. "
                         "Dùng --scale 0.01 để chạy thử toàn bộ ma trận trong "
                         "vài phút trước khi giao máy cho phiên qua đêm.")
    ap.add_argument("--turb", type=int, default=5, choices=range(1, 6),
                    help="mức nhiễu loạn cho MỤC A (quét cự ly). Mục B luôn "
                         "quét L1…L5; C/D cố định L3.")
    ap.add_argument("--sections", default="ABCD",
                    help="chọn mục cần đo, vd 'A' hoặc 'AB'")
    ap.add_argument("--cap-scale", type=float, default=1.0,
                    help="chỉ nhân TRẦN THỜI GIAN, giữ nguyên mục tiêu n_sift. "
                         "Cổng còn Latency Timer 16 ms (~63 qubit/s) thì dùng 4.")
    ap.add_argument("--chunk", type=int, default=1,
                    help="số lệnh qubit gửi liên tiếp trước khi đọc trả lời. "
                         ">1 CHỈ an toàn với bitstream có FIFO lệnh "
                         "(top_module.v [v11], CMD_DEPTH=64) và phải ≤ 64. "
                         "Kiểm tra bằng --chunk-check trước khi chạy phiên dài.")
    ap.add_argument("--qubit-us", type=float, default=220.0,
                    help="thời gian FPGA cần cho qubit KHÔNG click, µs (TX 50 + "
                         "rx_timeout 160 + gap 10 @ TURBO_SLOT=500). Tăng lên "
                         "nếu --chunk-check báo P_click thấp hơn mô hình.")
    ap.add_argument("--report-us", type=float, default=4400.0,
                    help="thời gian FPGA tốn THÊM cho mỗi qubit CÓ click "
                         "(FSM_REPORT chờ 220000 chu kỳ @50MHz kể từ [v13], khi "
                         "dòng báo cáo dài lên 40 byte). "
                         "Nhịp gửi = qubit_us + P_click·report_us.")
    ap.add_argument("--chunk-check", action="store_true",
                    help="đo P_click ở chunk=1 và chunk=--chunk rồi thoát. "
                         "Hai giá trị phải khớp nhau trong ~10%%; lệch nhiều "
                         "⇒ bitstream trên board CHƯA có FIFO lệnh, hoặc "
                         "--qubit-us đặt quá thấp.")
    ap.add_argument("--coh", type=int, default=11, choices=range(0, 16),
                    help="[v12] độ dài KHỐI KẾT HỢP tính theo SỐ QUBIT: "
                         "2^(coh+5). Mặc định 11 → 65536 qubit/khối, đúng bằng "
                         "một cửa sổ của channel_monitor (NEXP_LOG2=16). "
                         "coh=0 quay về bộ đếm CHU KỲ CLOCK cũ — chỉ dùng khi "
                         "muốn tái lập số liệu đã thu trước [v12]. "
                         "CẦN bitstream [v12]; bitstream cũ bỏ qua lệnh này.")
    ap.add_argument("--dyn-walk", dest="dyn_walk", action="store_true",
                    help="cho mức nhiễu loạn tự đi ngẫu nhiên ±1 (hành vi cũ với "
                         "turb>=2). MẶC ĐỊNH TẮT: khi đang quét chính mức nhiễu "
                         "loạn thì nó là biến không kiểm soát được.")
    ap.add_argument("--out", default=None,
                    help="thư mục số liệu (mặc định <repo>/data)")
    args = ap.parse_args()

    global DATA_DIR, POINTS_CSV
    if args.out:
        DATA_DIR = os.path.abspath(args.out)
        POINTS_CSV = os.path.join(DATA_DIR, "fpga_points.csv")

    os.makedirs(DATA_DIR, exist_ok=True)

    if args.merge:
        merge_sidecars()
        return
    if not args.port:
        ap.error("cần --port <COM> (hoặc --merge)")

    rng = np.random.default_rng(args.seed)
    link = Link(args.port, args.baud, args.timeout, args.qubit_us,
                args.report_us)

    matrix = [(t, wi, di, tb, li,
               max(1, int(tg * args.scale)), max(2.0, cap * args.scale),
               int(tq * args.scale) if tq else 0)
              for (t, wi, di, tb, li, tg, cap, tq)
              in build_matrix(args.phase, args.turb, 0,
                              args.sections.upper(), args.cap_scale)]
    done = load_done()
    todo = [m for m in matrix if m[0] not in done]

    print("█" * 78)
    print("  THU SỐ LIỆU FPGA — pha '%s', cổng %s" % (args.phase, args.port))
    print("  mục %s, nhiễu loạn mục A = L%d, chunk = %d"
          % (args.sections.upper(), args.turb, args.chunk))
    nq = coh_qubits(args.coh)
    print("  khối kết hợp: %s | mức nhiễu loạn %s"
          % ("bộ đếm CHU KỲ CLOCK (chế độ cũ, ~18 qubit/khối)" if nq == 0
             else "%d qubit (coh=%d)" % (nq, args.coh),
             "TỰ ĐI NGẪU NHIÊN ±1" if args.dyn_walk else "cố định"))
    print("  %d điểm, %d đã có, %d cần chạy" % (len(matrix), len(matrix) - len(todo),
                                                len(todo)))
    print("  trần thời gian tổng: %.1f giờ" % (sum(m[6] for m in todo) / 3600))
    print("█" * 78)

    try:
        if args.chunk_check:
            # Same point, several chunk sizes. A FIFO-less bitstream drops commands
            # that arrive closer than ~220 µs apart, and every dropped command is a
            # qubit the PC counted but the FPGA never ran — so P_click, measured as
            # n_click/n_qubit, collapses. Comparing against the MODEL (not just
            # against chunk=1) is what makes the test conclusive.
            pc_m, _ = model_expect("clear_ocean", 5.0, args.turb, 450)
            print("\n  Kiểm tra FIFO lệnh — clear ocean, d = 5 m, L%d, 450 nm"
                  % args.turb)
            print("      mô hình: P_click = %.3e\n" % pc_m)
            print("      %-8s %10s %10s %12s %10s"
                  % ("chunk", "qubit/s", "dòng/s", "P_click", "so mô hình"))
            sizes = sorted({1, 8, max(2, min(args.chunk, 64))})
            best, lines_per_s, ratios = (0.0, 1), [], []
            for ch in sizes:
                link.configure(0, 0, args.turb, 0)
                r = link.run_point(np.random.default_rng(args.seed), 10 ** 9,
                                   90.0, chunk=ch)
                rate = r["n_qubit"] / max(r["seconds"], 1e-9)
                lps = r["n_click"] / max(r["seconds"], 1e-9)
                pc = r["n_click"] / max(r["n_qubit"], 1)
                ratio = pc / pc_m
                lines_per_s.append(lps)
                ratios.append(ratio)
                print("      %-8d %10.0f %10.1f %12.3e %9.2f×%s"
                      % (ch, rate, lps, pc, ratio,
                         "" if 0.85 <= ratio <= 1.2 else "  ✗"))
                if 0.85 <= ratio <= 1.2 and rate > best[0]:
                    best = (rate, ch)

            print()
            if best[0]:
                print("  ✓ Tỉ số nằm trong 0.85–1.2 ⇒ không mất lệnh.")
                print("  → Nhanh nhất mà vẫn đúng: --chunk %d (%.0f qubit/s)"
                      % (best[1], best[0]))
                return

            # Diagnose. Dropped commands can only push the ratio DOWN: the PC counts
            # a qubit the FPGA never ran, so it lands in the denominator with no
            # chance of a click. A ratio well ABOVE 1 therefore is not a FIFO
            # problem at all — the link itself is clicking far too often.
            spread = (max(lines_per_s) - min(lines_per_s)) / max(min(lines_per_s), 1e-9)
            if max(ratios) > 1.5:
                print("  ✗ P_click CAO hơn mô hình. Mất lệnh chỉ kéo tỉ số xuống,")
                print("    nên đây KHÔNG phải lỗi FIFO — kênh đang click quá nhiều.")
                if spread < 0.1:
                    print("    Số dòng/giây gần như không đổi giữa các chunk (~%.0f/s)."
                          % (sum(lines_per_s) / len(lines_per_s)))
                    print("    1/4 ms = 250/s là trần của FSM_REPORT ⇒ MỌI qubit đều click.")
                print("    → SW[4] (chan_enable) đang TẮT: kênh UWOC bị bypass,")
                print("      uwoc_channel.v:315 ép click = 1. Gạt SW[4] = 1 rồi reset.")
                print("      Kiểm tra trên board: LEDG[3] sáng = kênh đang bật,")
                print("      LEDG[7] sáng = FIFO đã từng tràn.")
            else:
                print("  ✗ P_click THẤP hơn mô hình ⇒ FPGA đang mất lệnh.")
                print("    Bitstream chưa có FIFO lệnh, hoặc --qubit-us/--report-us")
                print("    đặt thấp quá. Hạ --chunk, hoặc tăng --qubit-us lên 300.")
            return

        if not args.skip_check:
            print("\n  Kiểm tra SW[1] (adaptive) qua độ nhạy của cfg_lambda...")
            ratio = check_switches(link, rng)
            print("      tỉ số P_click(450)/P_click(650) = %.2f  (fixed ⇒ ~3.0)"
                  % ratio)
            if args.phase == "fixed" and ratio < 1.8:
                print("\n  ✗ DỪNG: λ không ảnh hưởng P_click ⇒ SW[1] đang BẬT.")
                print("    Gạt SW[1] = 0 rồi chạy lại, hoặc dùng --phase adaptive.")
                return
            if args.phase == "adaptive" and ratio > 1.8:
                print("\n  ✗ DỪNG: λ vẫn ảnh hưởng P_click ⇒ SW[1] đang TẮT.")
                print("    Gạt SW[1] = 1 rồi chạy lại.")
                return
            print("      ✓ trạng thái công tắc khớp với pha '%s'" % args.phase)

        t_start = time.perf_counter()
        for k, (tag, wi, di, turb, li, target, cap, tq) in enumerate(todo, 1):
            water = WATER_ORDER[wi]
            d = float(WATER_TYPES[water].d_grid[di])
            lam = LAMBDA_OPTIONS[li]
            budget = ("%d qubit = %.0f khối" % (tq, tq / max(coh_qubits(args.coh), 1))
                      if tq and coh_qubits(args.coh)
                      else ("%d qubit" % tq if tq else "%d sift" % target))
            print("\n  [%d/%d] %s — %s, d = %g m, L%d, λ = %d nm "
                  "(mục tiêu %s, trần %.0f phút)"
                  % (k, len(todo), tag, WATER_TYPES[water].name, d, turb, lam,
                     budget, cap / 60))

            link.configure(wi, di, turb, li, coh=args.coh,
                           dyn_walk=args.dyn_walk)
            r = link.run_point(rng, target, cap, chunk=args.chunk,
                               target_qubit=tq)

            q = r["n_err"] / r["n_sift"] if r["n_sift"] else 0.5
            q_hi = qber_upper(r["n_err"], r["n_sift"])
            secure = (r["n_sift"] >= MIN_SIFT) and (q_hi < QBER_LIMIT)
            pc_m, q_m = model_expect(water, d, turb, lam)

            row = dict(
                tag=tag, phase=args.phase, water=water, dist_idx=di,
                distance_m=d, turb=turb, lam_nm=lam,
                n_qubit=r["n_qubit"], n_click=r["n_click"],
                n_sift=r["n_sift"], n_err=r["n_err"],
                p_click=r["n_click"] / max(r["n_qubit"], 1),
                sift_eff=r["n_sift"] / max(r["n_qubit"], 1),
                qber=q, qber_hi=q_hi, secure=secure,
                skr_per_pulse=((r["n_sift"] / max(r["n_qubit"], 1))
                               * max(0.0, 1 - 2 * h2(q_hi)) if secure else 0.0),
                seconds=r["seconds"], rate_qps=r["n_qubit"] / max(r["seconds"], 1e-9),
                p_click_model=pc_m, qber_model=q_m,
                coh_sel=args.coh, coh_qubits=coh_qubits(args.coh),
                dyn_walk=args.dyn_walk,
                timestamp=time.strftime("%Y-%m-%d %H:%M:%S"),
            )
            append_point(row)
            save_clicks(tag, r["clicks"])

            print("        n=%d  click=%d (P=%.3e, mô hình %.3e)  sift=%d  err=%d"
                  % (r["n_qubit"], r["n_click"], row["p_click"], pc_m,
                     r["n_sift"], r["n_err"]))
            # How many fading blocks the point actually spanned, and how many
            # sifted bits each carries. Below ~32 sifted per block the per-block
            # QBER is shot noise and only the pooled number is usable.
            # [v13] Which mode the controller actually sat in. Without this the
            # adaptive-vs-fixed comparison cannot be explained: a change in
            # P_click is only interpretable once you know whether mu was 6, 9 or
            # 12 fifteenths. Modes are reported per CLICK, so the fractions are
            # click-weighted, not time-weighted — a paused stretch emits nothing
            # and shows up as missing attempts instead.
            modes = [c["mode"] for c in r["clicks"] if "mode" in c]
            if modes:
                tot = len(modes)
                frac = ", ".join(
                    "%s %.0f%%" % (MODE_NAMES[m], 100 * modes.count(m) / tot)
                    for m in sorted(set(modes)))
                mus = [c["mu"] for c in r["clicks"] if "mu" in c]
                print("        chế độ: %s | μ trung bình = %.2f/15"
                      % (frac, sum(mus) / len(mus) if mus else float("nan")))

            # [v14] λ THỰC SỰ ĐƯỢC DÙNG. Trong pha adaptive, bộ leo đồi giữ quyền
            # chọn bước sóng, nên λ mà host gửi qua 0x40 KHÔNG phải λ của kênh.
            # Cột này từng vô dụng: res_lam trong top_module.v được khai báo và nối
            # vào bộ báo cáo nhưng chưa bao giờ được gán, nên lam_idx luôn bằng 0 —
            # đúng cái cột lẽ ra đã chỉ thẳng ra chỗ hỏng. Cần bitstream [v14].
            lams = [c["lam_idx"] for c in r["clicks"] if "lam_idx" in c]
            if lams and args.phase == "adaptive":
                tot = len(lams)
                lfrac = ", ".join(
                    "%d nm %.0f%%" % (LAMBDA_OPTIONS[l], 100 * lams.count(l) / tot)
                    for l in sorted(set(lams)))
                print("        λ thực dùng: %s%s"
                      % (lfrac,
                         "" if set(lams) == {li} else "  ← BỘ ĐIỀU KHIỂN ĐÃ ĐỔI λ"))

            nqb = coh_qubits(args.coh)
            if nqb:
                nblk = r["n_qubit"] / nqb
                print("        %.0f khối kết hợp (%d qubit/khối), "
                      "%.1f bit sàng/khối %s"
                      % (nblk, nqb, r["n_sift"] / max(nblk, 1e-9),
                         "" if r["n_sift"] / max(nblk, 1e-9) >= 32
                         else "← quá ít, chỉ dùng được QBER gộp"))
            print("        QBER = %.2f%%  (cận trên %.0f%% = %.2f%%, mô hình %.2f%%)  %s"
                  % (q * 100, CONF * 100, q_hi * 100, q_m * 100,
                     "an toàn" if secure else "KHÔNG/thiếu mẫu"))

            # A point whose click rate is far below the model is almost never
            # physics — it is the command FIFO overflowing, and it LATCHES: the
            # pacing uses p_hat = n_click / n_SENT, so once commands start being
            # dropped that estimate falls, the host paces faster, and more get
            # dropped. An adaptive section-A run lost 8 points to this before
            # anyone noticed, because nothing on screen said anything was wrong.
            ratio = row["p_click"] / max(pc_m, 1e-300)
            if ratio < 0.5 and r["n_click"] > 20:
                print("        " + "!" * 62)
                print("        ⚠ P_click chỉ bằng %.2f lần mô hình — ĐIỂM NÀY ĐÁNG NGỜ."
                      % ratio)
                # HAI nguyên nhân cho cùng một triệu chứng, và chúng cần hai cách
                # xử lý ngược nhau. Mất lệnh FIFO kéo P_click xuống theo một hệ số
                # KHÔNG phụ thuộc cự ly; bộ leo đồi kẹt ở λ sai kéo nó xuống theo
                # exp(−Δc·d), tức là càng xa càng tệ — và khớp CHÍNH XÁC với mô
                # hình ở một bước sóng khác. Thử khớp trước khi đổ cho FIFO.
                lam_lines = lambda_diagnosis(water, d, turb, lam, row["p_click"])
                if lam_lines:
                    for ln in lam_lines:
                        print("          " + ln)
                    print("            Cần bitstream [v14] (adaptive_controller.v: hủy")
                    print("            phép đo ứng viên khi cửa sổ hết hợp lệ + cfg_rst).")
                    print("            Với bitstream cũ: nhấn KEY[3] reset board rồi")
                    print("            chạy lại — trạng thái sai KHÔNG tự thoát được.")
                else:
                    print("          Nhiều khả năng FIFO lệnh đang tràn: kiểm tra LEDG[7]")
                    print("          trên board, hạ --chunk (chunk 1 chỉ chậm hơn ~20%),")
                    print("          hoặc tăng --qubit-us. Dừng phiên và chạy --chunk-check")
                    print("          trước khi để nó chạy tiếp hàng giờ.")
                print("        " + "!" * 62)
            print("        %.1f phút, còn lại ~%.1f giờ"
                  % (r["seconds"] / 60,
                     sum(m[6] for m in todo[k:]) / 3600))

        print("\n  Xong. Tổng %.2f giờ. Số liệu ở %s"
              % ((time.perf_counter() - t_start) / 3600, DATA_DIR))
    finally:
        link.close()


if __name__ == "__main__":
    main()
