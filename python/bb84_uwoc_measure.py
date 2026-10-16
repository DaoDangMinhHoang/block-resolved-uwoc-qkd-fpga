#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bb84_uwoc_measure.py — Measure BB84 over the UWOC channel (real FPGA or simulation)
═══════════════════════════════════════════════════════════════════════════════

Replaces bb84_fpga_qber_snr_5_level.py (which used the atmospheric Gamma-Gamma model).

MEASUREMENT MATRIX (see §4 of the proposal):
  Axis 1  water type × range  : clear ocean / coastal / harbor
  Axis 2  turbulence level    : L1…L5 (ε, χ_T, w per Table II of [2])
  Axis 3  system scenario     : fixed vs adaptive, λ sweep, μ sweep

UART PROTOCOL (PC → FPGA), bytes with bit[7] = 0:
    0x01                          reset the statistics
    0x02                          request a status report
    0x30 | dist[3:0]              set the range index
    0x40 | {water[1:0], lam[1:0]} set water type + wavelength
    0x50 | turb[2:0]              set the turbulence level
Bytes with bit[7] = 1 are qubit commands: bit[2]=a_data, bit[1]=a_basis, bit[0]=b_basis

Run:
    # simulation (no hardware needed)
    python python/bb84_uwoc_measure.py --simulate --scan distance
    python python/bb84_uwoc_measure.py --simulate --scan turbulence
    python python/bb84_uwoc_measure.py --simulate --scan wavelength
    python python/bb84_uwoc_measure.py --simulate --scan mu

    # real FPGA
    python python/bb84_uwoc_measure.py --port COM28 --scan distance --water clear_ocean --batch 20000 --turb 1

THROUGHPUT (FPGA mode only)
    Every qubit that does NOT click leaves the FPGA silent, so the PC must wait out
    `--timeout`. That is the entire measurement time. The FPGA itself needs only
    ~220 µs/qubit (TX 50 µs + rx_timeout 160 µs + gap 10 µs @ TURBO_SLOT = 500), so
    the timeout is the only thing worth tuning. To go faster still: lower the
    USB-serial port's Latency Timer to 1 ms (Device Manager → Ports → Advanced).
    FTDI's 16 ms default is exactly why you often see ~63 qubits/s.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
from scipy.stats import beta as _beta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from uwoc_channel_model import (            # noqa: E402
    WATER_TYPES, OCEAN_TURB_LEVELS, LAMBDA_OPTIONS, LinkConfig, UWOCChannel,
    h2, max_secure_distance,
)

WATER_ORDER = ("clear_ocean", "coastal", "harbor")
QBER_LIMIT = 0.11          # BB84 security limit
CONF = 0.95                # ONE-SIDED confidence level for the QBER upper bound

# Minimum sifted qubits for a QBER estimate to mean anything.
# MUST match the MIN_SIFT parameter of channel_monitor.v. Without this gate, a
# long-range measurement would report "QBER = 0%, secure" just because exactly
# one photon arrived — a statistically meaningless conclusion.
MIN_SIFT = 16

try:
    import serial
    HAS_SERIAL = True
except ImportError:
    HAS_SERIAL = False


# ═══════════════════════════════════════════════════════════════════════════
# Data source: real FPGA
# ═══════════════════════════════════════════════════════════════════════════
def _parse_full(line: str):
    """Parse one per-qubit report line from top_module.

    Format (top_module.v, per_qubit_reporter):
        @a_data,a_basis,b_basis,bob_bit,bmatch,err,irrad,total,sifted,errors*

    `total` is the ATTEMPT INDEX of this click since the last 0x01 reset — 6 hex
    digits from [v12] on, 4 before that. It is what lets the host assign a click to
    a coherence block: block = total >> (coh_sel + 5). Parsing it is optional so
    that logs captured with the older 4-hex bitstream still load; `attempt` is then
    simply absent.

    Returns a dict, or None if the line is not a valid report line.
    """
    line = line.strip()
    if not line.startswith("@") or "*" not in line:
        return None
    p = line[1:line.index("*")].split(",")
    if len(p) < 7:
        return None
    try:
        r = dict(a_data=int(p[0]), a_basis=int(p[1]), b_basis=int(p[2]),
                 bob=int(p[3]), bmatch=int(p[4]), err=int(p[5]),
                 irrad=int(p[6]))
    except ValueError:
        return None
    if len(p) >= 10:
        try:
            r["attempt"] = int(p[7], 16)
            r["tot_sift"] = int(p[8], 16)
            r["tot_err"] = int(p[9], 16)
        except ValueError:
            pass
    # [v13] adaptive mode + applied intensity, appended after the [v12] fields.
    if len(p) >= 13:
        try:
            r["mode"] = int(p[10])
            r["mu"] = int(p[11], 16)
            r["lam_idx"] = int(p[12])
        except ValueError:
            pass
    return r


def _parse(line: str):
    """Short form of [_parse_full]: → (bmatch, err) or None."""
    r = _parse_full(line)
    return None if r is None else (bool(r["bmatch"]), bool(r["err"]))


class FPGASource:
    def __init__(self, port: str, baud: int = 460800, chunk: int = 1,
                 timeout: float = 0.005):
        if not HAS_SERIAL:
            raise RuntimeError("Required: pip install pyserial")
        # A click report line is 34 bytes ⇒ ~3 ms @115200. The timeout only needs
        # to exceed that; the surplus is time WASTED on every no-click qubit
        # (>99% of qubits underwater).
        self.timeout = max(0.004, float(timeout))
        self.ser = serial.Serial(port, baud, timeout=self.timeout)
        self.chunk = max(1, chunk)
        time.sleep(0.5)
        self.ser.reset_input_buffer()

    def configure(self, water: int, dist: int, turb: int, lam: int):
        """Load the channel configuration into top_module's registers (see [v10])."""
        for b in (0x30 | (dist & 0x0F),
                  0x40 | ((water & 0x3) << 2) | (lam & 0x3),
                  0x50 | (turb & 0x07)):
            self.ser.write(bytes([b]))
            time.sleep(0.01)
        self.ser.write(bytes([0x01]))          # reset the statistics
        time.sleep(0.05)
        self.ser.reset_input_buffer()

    def collect(self, n: int, p_basis_z: float, seed: int) -> dict:
        """Collect n qubits over the FPGA.

        ⚠ ABOUT --chunk. The previous version defaulted to chunk = 64, arguing "the
        surplus bytes sit in the UART RX FIFO". WRONG: uart_rx.v has no FIFO, and
        top_module.v held the command in a SINGLE register. A byte arriving before
        the FSM reaches FSM_WAIT_CMD overwrites the previous one ⇒ lost commands.
        Firing back-to-back, the PC spaces bytes 86.8 µs apart (10 bits @115200)
        while the FPGA needs ≥220 µs per qubit, so only ~40% of commands ran —
        measured P_click came out ~2.8× below the model, sifted bits falling likewise.

        top_module.v [v11] added a 64-deep command FIFO. tb_cmd_fifo.v measured:
        the old RTL ran 14/32 commands, the new RTL 32/32 and 64/64 with no loss; at
        160 the FIFO overflows (127/160) but cmd_drop_cnt counts every drop and
        LEDG[7] lights up — loss is no longer silent. The default is still 1: NO
        slower, because the loop is bounded by readline's timeout, not by the FPGA.
        """
        rng = np.random.default_rng(seed)
        click = np.zeros(n, dtype=bool)
        bmatch = np.zeros(n, dtype=bool)
        derr = np.zeros(n, dtype=bool)
        t0 = time.perf_counter()

        # ── Pre-build every qubit command as a byte array ───────────
        ad_arr = rng.integers(0, 2, size=n, dtype=np.uint8)
        rand_ab = rng.random(n)
        rand_bb = rng.random(n)
        ab_arr = np.where(rand_ab < p_basis_z, 0, 1).astype(np.uint8)
        bb_arr = np.where(rand_bb < p_basis_z, 0, 1).astype(np.uint8)
        cmds = (np.uint8(0x80)
                | (ad_arr << 2)
                | (ab_arr << 1)
                | bb_arr)

        # ── Send + collect, chunk by chunk ───────────────────────────
        chunk = self.chunk
        done = 0
        for start in range(0, n, chunk):
            end = min(start + chunk, n)
            csize = end - start

            # Send the whole chunk of commands at once — only valid when the bitstream
            # has the command FIFO (top_module.v [v11]) and chunk ≤ FIFO depth.
            self.ser.write(bytes(cmds[start:end]))

            # Collect responses: each clicking qubit → 1 line "@....*\r\n"
            #                    each no-click qubit → timeout (no line)
            for j in range(csize):
                i = start + j
                p = _parse(self.ser.readline().decode("ascii", errors="ignore"))
                if p is not None:
                    click[i], bmatch[i], derr[i] = True, p[0], p[1]
                # no response ⇒ timeout ⇒ photon lost (no-click)

            done = end
            if done % 500 < chunk or done == n:
                el = time.perf_counter() - t0
                print(f"\r      [{done}/{n}] {done/el:.0f} q/s  "
                      f"click={click[:done].sum()}", end="", flush=True)

        # ── Drain the buffer ─────────────────────────────────────────
        # The USB-serial converter batches bytes by its Latency Timer (FTDI default
        # 16 ms), so a click's line can arrive AFTER that qubit's own readline has
        # timed out. It is then read by the NEXT readline: the click is still counted
        # exactly once and bmatch/err still come from that very line, only the index
        # shifts — harmless, since the statistics merely accumulate. The only thing
        # genuinely lost is the line of the LAST qubits in a batch, hence this drain.
        self.ser.timeout = max(0.05, 8 * self.timeout)
        free = n - 1
        while True:
            p = _parse(self.ser.readline().decode("ascii", errors="ignore"))
            if p is None:
                break
            while free >= 0 and click[free]:
                free -= 1
            if free < 0:
                break
            click[free], bmatch[free], derr[free] = True, p[0], p[1]
            free -= 1
        self.ser.timeout = self.timeout

        print("\r" + " " * 60 + "\r", end="")
        return _metrics(click, bmatch, derr, p_basis_z)

    def close(self):
        self.ser.close()


# ═══════════════════════════════════════════════════════════════════════════
# Data source: simulation (UWOC model)
# ═══════════════════════════════════════════════════════════════════════════
class SimSource:
    def __init__(self, cfg: LinkConfig):
        self.cfg = cfg
        self.water = "clear_ocean"
        self.dist_idx = 0
        self.turb = 1

    def configure(self, water: int, dist: int, turb: int, lam: int):
        self.water = WATER_ORDER[water]
        self.dist_idx = dist
        self.turb = turb
        self.cfg = LinkConfig(**{**self.cfg.__dict__,
                                 "lam_nm": LAMBDA_OPTIONS[lam]})

    def collect(self, n: int, p_basis_z: float, seed: int) -> dict:
        d = float(WATER_TYPES[self.water].d_grid[self.dist_idx])
        ch = UWOCChannel(self.water, d, self.turb, self.cfg)
        r = ch.simulate(n, seed=seed, p_basis_z=p_basis_z)
        return _metrics(r["click"], r["basis_match"], r["error"], p_basis_z)

    def close(self):
        pass


def qber_upper(n_err: int, n_sift: int, conf: float = CONF) -> float:
    """ONE-SIDED Clopper-Pearson upper bound on QBER at confidence `conf`.

    ⚠ Do NOT use the Wald standard error √(q(1−q)/n) as the previous version did:
    it COLLAPSES TO 0 when n_err = 0, so the test `q + se < 0.11` automatically
    says "secure" even with only a few dozen sifted bits. A real example: 0 errors
    / 93 sifted bits gives Wald = 0.00% ± 0.00%, while the true upper bound is
    3.17% (≈ the rule of three: 3/93). At d = 5 m the model QBER is 1.64% ⇒ 1.5
    expected errors ⇒ a 21% chance of seeing 0 — "QBER = 0%" is sampling scatter, not a clean channel.
    """
    if n_sift <= 0 or n_err >= n_sift:
        return 1.0
    return float(_beta.ppf(conf, n_err + 1, n_sift - n_err))


def _metrics(click, bmatch, derr, p_basis_z) -> dict:
    n = len(click)
    sift = click & bmatch
    n_sift = int(sift.sum())
    n_err = int((sift & derr).sum())
    q = n_err / n_sift if n_sift else 0.5
    q_hi = qber_upper(n_err, n_sift)
    valid = n_sift >= MIN_SIFT
    # "Secure" is only declared when the measurement has ENOUGH SAMPLES and the
    # UPPER BOUND on QBER — not the point estimate — stays below the limit. SKR
    # also uses the upper bound (pessimistic): every bit Eve could have taken counts.
    secure = valid and (q_hi < QBER_LIMIT)
    return dict(
        n=n, n_click=int(click.sum()), n_sift=n_sift, n_err=n_err,
        p_click=float(click.mean()), sift_eff=n_sift / n, qber=q,
        qber_hi=q_hi, valid=valid, secure=secure,
        skr_per_pulse=(n_sift / n) * max(0.0, 1 - 2 * h2(q_hi)) if secure else 0.0,
    )


# ═══════════════════════════════════════════════════════════════════════════
# The sweeps
# ═══════════════════════════════════════════════════════════════════════════
def _hdr(t):
    print("\n" + "═" * 82)
    print("  " + t)
    print("═" * 82)


def _row_hdr():
    print(f"  {'configuration':<26}{'click':>8}{'sift':>7}{'P_click':>11}"
          f"{'QBER':>9}{'upper bound':>12}{'SKR/pulse':>11}{'verdict':>12}")
    print("  " + "-" * 84)


def _row(tag, m):
    if not m["valid"]:
        verdict = f"n/a (<{MIN_SIFT})"
    elif m["secure"]:
        verdict = "secure"
    else:
        verdict = "NOT SECURE"
    print(f"  {tag:<26}{m['n_click']:>8}{m['n_sift']:>7}{m['p_click']:>11.3e}"
          f"{m['qber']*100:>8.2f}%{m['qber_hi']*100:>11.2f}%"
          f"{m['skr_per_pulse']:>11.3e}{verdict:>12}")


def scan_distance(src, args, adaptive_pz):
    _hdr(f"DISTANCE SWEEP — {WATER_TYPES[args.water].name}, turbulence L{args.turb}")
    print(f"  c(λ) = {WATER_TYPES[args.water].c:.3f} m⁻¹, "
          f"λ = {LAMBDA_OPTIONS[args.lam]} nm, batch = {args.batch}")
    _row_hdr()
    wi = WATER_ORDER.index(args.water)
    grid = WATER_TYPES[args.water].d_grid
    last_secure = 0.0
    stop = False
    for di in range(0, len(grid), args.step):
        src.configure(wi, di, args.turb, args.lam)
        m = src.collect(args.batch, args.pz, args.seed + di)
        _row(f"d = {grid[di]:g} m", m)
        if m["secure"] and not stop:
            last_secure = grid[di]
        # The secure range is the last CONTIGUOUS point still secure: once a
        # measurement becomes invalid or insecure, farther points do not count,
        # even if they happen to show a low QBER purely from too few samples.
        if not m["secure"]:
            stop = True
    print(f"\n  Maximum secure distance (contiguous): {last_secure:g} m")
    print(f"  Batch = {args.batch}. {MIN_SIFT} sifted bits suffice only to AVOID a rash conclusion,")
    print("  not to MEASURE. For ~0.4% error around QBER ~1.6%, n_sift ≳ 1000 is needed,")
    print("  i.e. batch ≳ 1000/(P_click·q_basis) — about 2·10⁵ at short range, more at longer ranges.")


def scan_turbulence(src, args, adaptive_pz):
    _hdr(f"TURBULENCE SWEEP — {WATER_TYPES[args.water].name}, "
         f"d = {WATER_TYPES[args.water].d_grid[args.dist]:g} m")
    _row_hdr()
    wi = WATER_ORDER.index(args.water)
    for lv in range(1, 6):
        src.configure(wi, args.dist, lv, args.lam)
        m = src.collect(args.batch, args.pz, args.seed + lv)
        _row(f"L{lv} {OCEAN_TURB_LEVELS[lv]['name']}", m)
    print("\n  Note: MEAN QBER is not very sensitive to the turbulence level (P_click is nearly")
    print("  linear in h, so fading cancels in expectation). Turbulence appears in the")
    print("  VARIABILITY between windows — see verify_window_statistics() in")
    print("  uwoc_channel_model.py and the qber_jitter output of channel_monitor.")


def scan_wavelength(src, args, adaptive_pz):
    _hdr("WAVELENGTH SWEEP — optimal λ depends on water type")
    print("  Expected from [1] §I: clear ocean favors blue-green; turbid water shifts toward red.")
    for water in WATER_ORDER:
        wi = WATER_ORDER.index(water)
        # pick a range around 60% of the secure range so the link survives at every λ
        grid = WATER_TYPES[water].d_grid
        di = max(0, min(len(grid) - 1, len(grid) // 4))
        print(f"\n  ▸ {WATER_TYPES[water].name}, d = {grid[di]:g} m")
        _row_hdr()
        best, best_skr = None, -1.0
        for li, lam in enumerate(LAMBDA_OPTIONS):
            src.configure(wi, di, args.turb, li)
            m = src.collect(args.batch, args.pz, args.seed + li)
            _row(f"λ = {lam} nm", m)
            if m["skr_per_pulse"] > best_skr:
                best_skr, best = m["skr_per_pulse"], lam
        print(f"    → best λ: {best} nm")


def scan_mu(src, args, adaptive_pz):
    _hdr("INTENSITY SWEEP μ — key-rate versus PNS-risk trade-off")
    print("  Higher μ lowers QBER (n̄ ≫ Y₀), but increases the MULTI-PHOTON pulse fraction.")
    print("  Enormous underwater attenuation makes it far easier for Eve to conceal a PNS attack")
    print("  than in the atmosphere ⇒ μ must be capped (MU_CAP in the controller).\n")
    wi = WATER_ORDER.index(args.water)
    print(f"  {'μ':>7}{'P(n≥2 | n≥1)':>16}{'P_click':>12}{'QBER':>9}"
          f"{'SKR/pulse':>12}{'secure':>9}")
    print("  " + "-" * 68)
    base = src.cfg if hasattr(src, "cfg") else LinkConfig()
    for mu in (0.05, 0.1, 0.2, 0.4, 0.8):
        if hasattr(src, "cfg"):
            src.cfg = LinkConfig(**{**base.__dict__, "mu": mu})
        src.configure(wi, args.dist, args.turb, args.lam)
        m = src.collect(args.batch, args.pz, args.seed)
        # P(multi-photon | at least 1) for a weak coherent source (Poisson)
        p_multi = 1 - mu * np.exp(-mu) / (1 - np.exp(-mu))
        print(f"  {mu:>7.2f}{p_multi*100:>15.2f}%{m['p_click']:>12.3e}"
              f"{m['qber']*100:>8.2f}%{m['skr_per_pulse']:>12.3e}"
              f"{'✓' if m['secure'] else '✗':>9}")
    if hasattr(src, "cfg"):
        src.cfg = base


SCANS = {"distance": scan_distance, "turbulence": scan_turbulence,
         "wavelength": scan_wavelength, "mu": scan_mu}


def main():
    ap = argparse.ArgumentParser(description="Measure BB84 over the UWOC channel")
    ap.add_argument("--port", help="FPGA COM port (e.g. COM28)")
    # [Q2] 460800, matching uart_rx/uart_tx in top_module.v. The old 115200
    # default silently decodes nothing on the current bitstream.
    ap.add_argument("--baud", type=int, default=460800)
    ap.add_argument("--simulate", action="store_true")
    ap.add_argument("--scan", default="distance", choices=list(SCANS))
    ap.add_argument("--water", default="clear_ocean", choices=WATER_ORDER)
    ap.add_argument("--dist", type=int, default=2, help="range index 0–15")
    ap.add_argument("--turb", type=int, default=3, choices=range(1, 6))
    ap.add_argument("--lam", type=int, default=0, choices=range(3),
                    help="0=450nm, 1=532nm, 2=650nm")
    ap.add_argument("--pz", type=float, default=0.5, help="P(choosing the Z basis)")
    ap.add_argument("--mu", type=float, default=0.1)
    ap.add_argument("--batch", type=int, default=50000)
    ap.add_argument("--chunk", type=int, default=1,
                    help="number of consecutive qubit commands sent (FPGA mode). The default is 1 "
                         "and is NOT slower — the loop is limited by --timeout, "
                         "not by the FPGA. Chunk > 1 is safe only with "
                         "top_module.v [v11] (64-deep command FIFO; verified up to 64 "
                         "with tb_cmd_fifo.v). Large chunks on an older bitstream "
                         "silently DROP commands and force P_click down to ~0.4 times its value.")
    ap.add_argument("--timeout", type=float, default=0.005,
                    help="time to wait for each qubit report line, in seconds (FPGA "
                         "mode). A 34-byte line takes ≈ 3 ms @115200. Late lines "
                         "are not lost (the next loop iteration reads them, and the "
                         "batch ends with a buffer drain). If the USB-serial port's "
                         "Latency Timer is still at its default 16 ms, lower it to 1 ms to see an effect.")
    ap.add_argument("--step", type=int, default=2, help="range-sweep step")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if not args.port and not args.simulate:
        ap.error("--port <COM> or --simulate is required")

    cfg = LinkConfig(lam_nm=LAMBDA_OPTIONS[args.lam], mu=args.mu)
    src = (FPGASource(args.port, args.baud, chunk=args.chunk,
                      timeout=args.timeout)
           if args.port else SimSource(cfg))
    tag = "FPGA" if args.port else "SIM"

    print("█" * 82)
    print(f"  BB84 OVER THE UWOC CHANNEL — source: {tag}")
    print(f"  μ={cfg.mu}  η_det={cfg.eta_det}  Y₀={cfg.Y0:.2e}  "
          f"f_rep={cfg.f_rep:.0e} Hz  p_z={args.pz}")
    print(f"  QBER upper bound = one-sided Clopper–Pearson {CONF*100:.0f}%; "
          f"'secure' ⟺ upper bound < {QBER_LIMIT*100:.0f}%")
    if args.port:
        print(f"  timeout={args.timeout*1000:.0f} ms  chunk={args.chunk}"
              + ("   ⚠ chunk>1 requires a bitstream with a command FIFO (top_module.v v11)"
                 if args.chunk > 1 else ""))
    print("█" * 82)

    try:
        SCANS[args.scan](src, args, args.pz)
    finally:
        src.close()

    if args.simulate:
        _hdr("ANALYTICAL REFERENCE — maximum secure distance d_max [m]")
        print(f"  {'Water':<17}" + "".join(f"{'L'+str(l):>9}" for l in range(1, 6)))
        print("  " + "-" * 64)
        for w in WATER_ORDER:
            row = [max_secure_distance(w, lv, cfg, criterion="outage")
                   for lv in range(1, 6)]
            print(f"  {WATER_TYPES[w].name:<17}"
                  + "".join(f"{v:>9.2f}" for v in row))
    print()


if __name__ == "__main__":
    main()
