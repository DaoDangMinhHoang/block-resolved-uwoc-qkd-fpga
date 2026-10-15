#!/usr/bin/env python3
"""
send_security_commands.py (rev E) -- configure the DE1 before a capture.

ORDER MATTERS. Fast mode and the rate preset go BEFORE the statistics reset, so
that uwoc_channel's stats_rst restarts qub_index at the same moment
security_block_monitor clears. That is what aligns the security block grid with
the fading block grid, and it is why h_min == h_max is an assertion rather
than a filter. Reset must be the last byte sent.

Real-time run (this is the one the paper's title claims):

    python send_security_commands.py --port COM5 \
        --dist 4 --water 0 --lambda-sel 0 --turb 5 --coh 11 \
        --y0 nominal --fast on --rate 10M --reset

--coh 11 gives 2^(11+5) = 65536 attempts per fading block, matching
ATTEMPT_LOG2 = 16. Any other value breaks the h_min == h_max invariant.
--rate 10M paces attempts at the manuscript's f_rep, so one fading block
takes 65536/1e7 = 6.55 ms of wall clock and the emulator runs in real time.

[rev E, needs the rev B bitstream]
  --turb 6 / 7      matched-variance family controls (Weibull at the L3 index,
                    log-normal at the L4 index). No range offset is applied, so
                    use them at dist_idx 4 (25 m), where the offset is 0.
  --y0 x0.25 ... x16  background-floor dose-response presets (Y0 = 218/2^24 x k)
  --decoy on        interleaved signal / decoy / vacuum attempts
  --decoy-preset P  class probabilities: 0 = 1/2,1/4,1/4  1 = .7,.2,.1
                    2 = .8,.15,.05  3 = 1/3 each
  --decoy-shift S   decoy intensity = signal / 2^S, S = 1..4
  --mu-gain G       signal intensity x G, G = 1, 2, 4, 8 (mu_level stays 8)
"""

import argparse
import time

CMD_RESET = 0x01
CMD_FAST_OFF = 0x70
CMD_FAST_ON = 0x71
CMD_Y0_ZERO = 0x72
CMD_Y0_NOMINAL = 0x73
CMD_DECOY_OFF = 0x74
CMD_DECOY_ON = 0x75

# Rate presets. The single-byte command protocol cannot carry a 32-bit value, so
# top_module holds a ROM indexed by the low three bits.
#
# 0x78-0x7F, NOT 0x80-0x87. A qubit command is
#     0x80 | (data << 2) | (abasis << 1) | bbasis
# (fpga_collect.py), so it occupies exactly 0x80-0x87. Putting the rate presets
# there makes the FPGA swallow every qubit command the host sends: PC mode then
# collects nothing at all and tb_cmd_fifo.v reports 0 of 32 commands run.
RATE_PRESETS = {
    "free": (0x78, "free run, as fast as the channel answers"),
    "20k": (0x79, "20 kattempt/s"),
    "100k": (0x7A, "100 kattempt/s"),
    "1M": (0x7B, "1 Mattempt/s"),
    "2M": (0x7C, "2 Mattempt/s"),
    "5M": (0x7D, "5 Mattempt/s"),
    "10M": (0x7E, "10 Mattempt/s = manuscript f_rep, REAL TIME"),
    "13M1": (0x7F, "13.1 Mattempt/s = 65536 pulses per 5 ms"),
}

# [rev E] background-floor presets 0x18-0x1F: Y0 = 218 x factor.
Y0_PRESETS = {"x0.25": 1, "x0.5": 2, "x1": 3, "x2": 4, "x4": 5, "x8": 6, "x16": 7}


def send(ser, value, label, delay=0.03):
    ser.write(bytes([value & 0xFF]))
    ser.flush()
    print("sent 0x%02X  %s" % (value, label))
    time.sleep(delay)


def main():
    import serial

    ap = argparse.ArgumentParser()
    ap.add_argument("--port", required=True)
    ap.add_argument("--baud", type=int, default=460800)
    ap.add_argument("--dist", type=int, choices=range(16))
    ap.add_argument("--water", type=int, choices=range(4))
    ap.add_argument("--lambda-sel", type=int, choices=range(4), dest="lam")
    ap.add_argument("--turb", type=int, choices=range(8),
                    help="0 off, 1-5 levels, 6/7 ablation rows (rev B)")
    ap.add_argument("--dynamic", action="store_true")
    ap.add_argument("--coh", type=int, choices=range(16))
    ap.add_argument("--y0", choices=["zero", "nominal"] + sorted(Y0_PRESETS))
    ap.add_argument("--decoy", choices=["on", "off"])
    ap.add_argument("--decoy-preset", type=int, choices=range(4))
    ap.add_argument("--decoy-shift", type=int, choices=range(1, 5))
    ap.add_argument("--mu-gain", type=int, choices=(1, 2, 4, 8))
    ap.add_argument("--fast", choices=["on", "off"])
    ap.add_argument("--rate", choices=sorted(RATE_PRESETS))
    ap.add_argument("--reset", action="store_true")
    args = ap.parse_args()

    if (args.water is None) != (args.lam is None):
        raise SystemExit("--water and --lambda-sel must be given together")
    if (args.decoy_preset is None) != (args.decoy_shift is None):
        raise SystemExit("--decoy-preset and --decoy-shift must be given together")
    if args.coh is not None and args.coh != 11:
        print("WARNING: coh_sel = %d gives 2^%d attempts per fading block. "
              "The monitor uses 2^16, so h_min == h_max will fail on every "
              "block." % (args.coh, args.coh + 5))
    if args.turb in (6, 7) and args.dist is not None and args.dist != 4:
        print("WARNING: ablation rows carry no range offset; they are matched to "
              "rows 3/4 only where the offset is 0 (clear ocean 15-25 m, "
              "dist_idx 2-4).")

    with serial.Serial(args.port, args.baud, timeout=1.0) as ser:
        if args.dist is not None:
            send(ser, 0x30 | args.dist, "dist_idx=%d" % args.dist)
        if args.water is not None:
            send(ser, 0x40 | ((args.water & 3) << 2) | (args.lam & 3),
                 "water=%d lambda=%d" % (args.water, args.lam))
        if args.turb is not None:
            send(ser, 0x50 | ((1 if args.dynamic else 0) << 3) | (args.turb & 7),
                 "turb=%d dynamic=%s" % (args.turb, args.dynamic))
        if args.coh is not None:
            send(ser, 0x60 | args.coh, "coh_sel=%d" % args.coh)
        if args.y0 == "zero":
            send(ser, CMD_Y0_ZERO, "Y0 = 0 (radial-null control)")
        elif args.y0 == "nominal":
            send(ser, CMD_Y0_NOMINAL, "Y0 = nominal")
        elif args.y0 is not None:
            send(ser, 0x18 | Y0_PRESETS[args.y0], "Y0 = nominal %s" % args.y0)
        if args.decoy_preset is not None:
            send(ser, 0x20 | (args.decoy_preset << 2) | (args.decoy_shift - 1),
                 "decoy preset %d, decoy = signal / %d"
                 % (args.decoy_preset, 1 << args.decoy_shift))
        if args.mu_gain is not None:
            g = {1: 0, 2: 1, 4: 2, 8: 3}[args.mu_gain]
            send(ser, 0x14 | g, "signal intensity gain x%d" % args.mu_gain)
        if args.decoy == "on":
            send(ser, CMD_DECOY_ON, "decoy intensities ON")
        elif args.decoy == "off":
            send(ser, CMD_DECOY_OFF, "decoy intensities OFF")
        if args.rate:
            code, desc = RATE_PRESETS[args.rate]
            send(ser, code, "rate preset %s: %s" % (args.rate, desc))
        if args.fast == "on":
            send(ser, CMD_FAST_ON, "security fast mode ON")
        elif args.fast == "off":
            send(ser, CMD_FAST_OFF, "security fast mode OFF")
        if args.reset:
            send(ser, CMD_RESET, "statistics + fading-block grid reset (LAST)")

    print("done")


if __name__ == "__main__":
    main()
