#!/usr/bin/env python3
"""
capture_security_blocks.py (rev E) -- read binary block packets from the DE1
and write a CSV that phase2_security_analysis.py / decoy_security_analysis.py
consume.

Packet versions:
  v5   96 bytes  rev A bitstream (Q2 campaign 2026-09)
  v6  144 bytes  rev B bitstream: same offsets 0..92, plus the decoy
                 configuration byte and per-class counters for decoy (class 1)
                 and vacuum (class 2) attempts. The totals still count every
                 attempt; class 0 (signal) = total - class 1 - class 2.

Every packet is validated before it is written. `capture_valid` is what the
analysis filters on; `validation_error` records why a packet failed, so a bad
run is diagnosable afterwards instead of silently producing a flat curve.

The gates that have actually caught bugs in this project:

  trials > clicks STRICTLY
      If evt_bmatch is asserted only on clicks, trials == clicks on every block
      and every gain comes out as 1.0, while the naive ordering chain
      errors <= clicks <= trials <= attempts still passes.

  block_synced == 1
      An earlier monitor revision cleared its internal sync flag on the same
      edge that raised block_pulse, so the reporter latched zero one clock
      later and this gate rejected the entire capture. That is the correct
      outcome: a silent flag is worse than a loud rejection.

  h_min == h_max
      All-or-nothing, not partial. With coh_sel = 11 the security window and the
      fading block are the same grid, so this holds on every block or on none.
      A 60% pass rate means a config or alignment bug, not 60% usable data.

    python capture_security_blocks.py --port COM5 --baud 460800 \
        --out l5_nominal.csv --blocks 2000
"""

import argparse
import csv
import struct
import time

MAGIC = b"\xA5\x5A"
PACKET_LEN = {5: 96, 6: 144}

CLASS_FIELDS = ("attempts", "clicks_all", "trials_z", "clicks_z", "errors_z",
                "trials_x", "clicks_x", "errors_x")
# Class probabilities of realtime_attempt_engine.v rev F, by decoy preset.
DECOY_PRESETS = {0: (0.50, 0.25, 0.25), 1: (0.70, 0.20, 0.10),
                 2: (0.80, 0.15, 0.05), 3: (1 / 3, 1 / 3, 1 / 3)}


def xor_checksum(data):
    v = 0
    for b in data:
        v ^= b
    return v


def parse_packet(p):
    if p[0:2] != MAGIC:
        raise ValueError("magic")
    ver = p[2]
    if ver not in PACKET_LEN:
        raise ValueError("version %d (expected 5 or 6)" % ver)
    n = PACKET_LEN[ver]
    if len(p) != n:
        raise ValueError("length %d for version %d" % (len(p), ver))
    exp = xor_checksum(p[:n - 1])
    if p[n - 1] != exp:
        raise ValueError("checksum got 0x%02X expected 0x%02X" % (p[n - 1], exp))

    def u32(o):
        return struct.unpack_from("<I", p, o)[0]

    def u16(o):
        return struct.unpack_from("<H", p, o)[0]

    def u24(o):
        return p[o] | (p[o + 1] << 8) | (p[o + 2] << 16)

    f = p[3]
    row = {
        "version": ver,
        "fast_mode": (f >> 0) & 1,
        "adaptive": (f >> 1) & 1,
        "has_sift": (f >> 2) & 1,
        "config_changed": (f >> 3) & 1,
        "block_synced": (f >> 4) & 1,
        "realtime_paced": (f >> 5) & 1,
        "rate_valid": (f >> 6) & 1,
        "block_id": u32(4),
        "attempts": u32(8),
        "clicks_all": u32(12),
        "lost_all": u32(16),
        "trials_z": u32(20),
        "clicks_z": u32(24),
        "errors_z": u32(28),
        "trials_x": u32(32),
        "clicks_x": u32(36),
        "errors_x": u32(40),
        "h_mean": u16(44),
        "h_min": u16(46),
        "h_max": u16(48),
        "water": p[50],
        "dist_idx": p[51],
        "turb": p[52],
        "lambda_sel": p[53],
        "mu_level": p[54],
        "adapt_mode": p[55],
        "p_noise": p[56] | (p[57] << 8) | (p[58] << 16),
        "coh_sel": p[59],
        "slot_width": u32(60),
        "attempts_per_sec": u32(64),
        "rate_hz_requested": u32(68),
        "dropped_blocks": u32(72),
        "protocol_errors": u32(76),
        "engine_timeouts": u32(80),
        "engine_overruns": u32(84),
        "watchdog_aborts": u32(88),
        "attempt_log2": p[92],
    }
    if ver >= 6:
        d = p[93]
        row.update(decoy_enable=d & 1, decoy_shift=((d >> 1) & 3) + 1,
                   decoy_preset=(d >> 3) & 3, mu_gain_log2=(d >> 5) & 3)
        for k, name in enumerate(CLASS_FIELDS):
            row["dec_" + name] = u24(95 + 3 * k)
            row["vac_" + name] = u24(119 + 3 * k)
    else:
        row.update(decoy_enable=0, decoy_shift=0, decoy_preset=0, mu_gain_log2=0)
        for name in CLASS_FIELDS:
            row["dec_" + name] = 0
            row["vac_" + name] = 0
    return row


def read_packet(ser):
    """Resynchronise on A5 5A, then read the remainder by version."""
    state = 0
    while True:
        b = ser.read(1)
        if not b:
            return None
        if state == 0:
            if b == b"\xA5":
                state = 1
        else:
            if b == b"\x5A":
                v = ser.read(1)
                if not v:
                    return None
                n = PACKET_LEN.get(v[0])
                if n is None:
                    state = 0
                    continue
                rest = ser.read(n - 3)
                if len(rest) != n - 3:
                    return None
                return MAGIC + v + rest
            state = 1 if b == b"\xA5" else 0


def _ordering_errors(prefix, row, attempts_key):
    errs = []
    t = row[prefix + "trials_z"] + row[prefix + "trials_x"]
    c = row[prefix + "clicks_z"] + row[prefix + "clicks_x"]
    e = row[prefix + "errors_z"] + row[prefix + "errors_x"]
    label = prefix.rstrip("_") or "pooled"
    if not (e <= c <= t <= row[attempts_key]):
        errs.append("%s: E <= C <= A <= N failed" % label)
    for b in ("z", "x"):
        if not (row[prefix + "errors_" + b] <= row[prefix + "clicks_" + b]
                <= row[prefix + "trials_" + b]):
            errs.append("%s: errors_%s <= clicks_%s <= trials_%s failed"
                        % (label, b, b, b))
    if row[prefix + "clicks_all"] > row[attempts_key]:
        errs.append("%s: clicks_all exceeds attempts" % label)
    if c > row[prefix + "clicks_all"]:
        errs.append("%s: sifted clicks exceed clicks_all" % label)
    return errs


def validate(row, rate_tol):
    errs = []
    warns = []
    n = 1 << row["attempt_log2"]
    trials = row["trials_z"] + row["trials_x"]
    clicks = row["clicks_z"] + row["clicks_x"]

    if row["attempts"] != n:
        errs.append("attempts %d != 2^%d" % (row["attempts"], row["attempt_log2"]))
    if row["clicks_all"] + row["lost_all"] != row["attempts"]:
        errs.append("clicks_all + lost_all != attempts")

    # Pooled ordering, then per basis: a pooled check can hide a broken Z or X
    # counter when the other one compensates.
    errs += _ordering_errors("", row, "attempts")
    if trials <= clicks:
        errs.append("FATAL: trials <= clicks -- evt_bmatch is not valid on "
                    "lost attempts, see PATCH item 3")

    if not (row["h_min"] <= row["h_mean"] <= row["h_max"]):
        errs.append("h ordering")
    if row["h_min"] != row["h_max"]:
        errs.append("h_min != h_max -- block spans more than one fading "
                    "realisation; check coh_sel == 11 and reset ordering")

    if not row["block_synced"]:
        errs.append("block_synced = 0")
    if row["config_changed"]:
        errs.append("configuration changed inside block")
    if row["protocol_errors"]:
        errs.append("protocol_errors=%d" % row["protocol_errors"])
    if row["dropped_blocks"]:
        errs.append("dropped_blocks=%d" % row["dropped_blocks"])
    if row["engine_timeouts"]:
        errs.append("engine_timeouts=%d (channel did not answer)"
                    % row["engine_timeouts"])
    if row["watchdog_aborts"]:
        errs.append("watchdog_aborts=%d (UART ownership lost mid-packet)"
                    % row["watchdog_aborts"])

    # [rev E] intensity classes
    n_dec, n_vac = row["dec_attempts"], row["vac_attempts"]
    if row["version"] >= 6:
        for prefix in ("dec_", "vac_"):
            errs += _ordering_errors(prefix, row, prefix + "attempts")
        if n_dec + n_vac > row["attempts"]:
            errs.append("class attempts exceed the block total")
        for name in CLASS_FIELDS:
            if row["dec_" + name] + row["vac_" + name] > row[name]:
                errs.append("class counters exceed total for %s" % name)
        if not row["decoy_enable"] and (n_dec or n_vac):
            errs.append("decoy disabled but class counters are nonzero")
        if row["decoy_enable"]:
            ps, pd_, pv = DECOY_PRESETS[row["decoy_preset"]]
            for label, cnt, p in (("decoy", n_dec, pd_), ("vacuum", n_vac, pv)):
                sd = (n * p * (1 - p)) ** 0.5
                if abs(cnt - n * p) > 6 * sd:
                    warns.append("%s attempts %d far from %.0f (preset %d)"
                                 % (label, cnt, n * p, row["decoy_preset"]))
            if row["vac_clicks_all"] and row["p_noise"] == 0:
                errs.append("vacuum attempts clicked with Y0 = 0")

    # Warnings do not invalidate a block, but they must be visible.
    if row["engine_overruns"]:
        warns.append("engine_overruns=%d: requested rate exceeds what the "
                     "channel sustains; trust attempts_per_sec, not rate_hz"
                     % row["engine_overruns"])
    # With independent unbiased bases P(match) = 1/2, so the TOTAL is N/2 and
    # each basis is N/4. Testing only against N/4 would accept a build in which
    # one basis counter is dead, because the survivor lands exactly there.
    n_att = row["attempts"]
    if not (0.45 * n_att <= trials <= 0.55 * n_att):
        warns.append("total basis matches %d far from N/2 = %d"
                     % (trials, n_att // 2))
    for name in ("trials_z", "trials_x"):
        if not (0.20 * n_att <= row[name] <= 0.30 * n_att):
            warns.append("%s = %d far from N/4 = %d; check the basis bias and "
                         "the Z/X polarity" % (name, row[name], n_att // 4))

    req = row["rate_hz_requested"]
    if not row["rate_valid"]:
        warns.append("rate_valid=0: the first one-second window has not closed, "
                     "attempts_per_sec is meaningless for this block")
    elif req and abs(row["attempts_per_sec"] - req) > rate_tol * req:
        warns.append("measured %d/s vs requested %d/s: pacing not met"
                     % (row["attempts_per_sec"], req))
    return errs, warns


FIELDS = [
    "host_time", "version", "fast_mode", "adaptive", "has_sift",
    "config_changed", "block_synced", "realtime_paced", "rate_valid",
    "block_id", "block_gap",
    "attempts", "clicks_all", "lost_all",
    "trials_z", "clicks_z", "errors_z",
    "trials_x", "clicks_x", "errors_x",
    "h_mean", "h_min", "h_max",
    "water", "dist_idx", "turb", "lambda_sel", "mu_level", "adapt_mode",
    "p_noise", "coh_sel", "slot_width",
    "attempts_per_sec", "rate_hz_requested",
    "dropped_blocks", "protocol_errors", "engine_timeouts",
    "engine_overruns", "watchdog_aborts",
    "attempt_log2",
    "decoy_enable", "decoy_shift", "decoy_preset", "mu_gain_log2",
] + ["dec_" + f for f in CLASS_FIELDS] + ["vac_" + f for f in CLASS_FIELDS] + [
    "capture_valid", "validation_error", "validation_warn",
]


def main():
    import serial

    ap = argparse.ArgumentParser()
    ap.add_argument("--port", required=True)
    ap.add_argument("--baud", type=int, default=460800)
    ap.add_argument("--out", default="fpga_security_blocks.csv")
    ap.add_argument("--blocks", type=int, default=0, help="0 = until Ctrl+C")
    ap.add_argument("--timeout", type=float, default=10.0)
    ap.add_argument("--print-every", type=int, default=1)
    ap.add_argument("--rate-tol", type=float, default=0.02)
    args = ap.parse_args()

    count = good = 0
    prev_id = None

    with serial.Serial(args.port, args.baud, timeout=args.timeout) as ser, \
            open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        print("listening on %s @ %d, writing %s" % (args.port, args.baud, args.out))

        try:
            while args.blocks == 0 or count < args.blocks:
                raw = read_packet(ser)
                if raw is None:
                    print("UART timeout, no complete packet")
                    continue
                try:
                    row = parse_packet(raw)
                except ValueError as exc:
                    print("packet rejected: %s" % exc)
                    continue

                gap = 0 if prev_id is None else \
                    (row["block_id"] - prev_id - 1) & 0xFFFFFFFF
                prev_id = row["block_id"]

                errs, warns = validate(row, args.rate_tol)
                if gap:
                    errs.append("block_id gap=%d" % gap)

                row["block_gap"] = gap
                row["capture_valid"] = 0 if errs else 1
                row["validation_error"] = " | ".join(errs)
                row["validation_warn"] = " | ".join(warns)
                row["host_time"] = time.time()
                w.writerow(row)
                f.flush()
                count += 1
                good += row["capture_valid"]

                if count % max(args.print_every, 1) == 0 or errs or warns:
                    trials = row["trials_z"] + row["trials_x"]
                    clicks = row["clicks_z"] + row["clicks_x"]
                    errors = row["errors_z"] + row["errors_x"]
                    g = clicks / trials if trials else float("nan")
                    q = errors / clicks if clicks else float("nan")
                    extra = ""
                    if row["decoy_enable"]:
                        extra = " dec=%d vac=%d" % (row["dec_attempts"],
                                                    row["vac_attempts"])
                    print("[%s] v%d blk=%7d A=%6d C=%5d E=%4d G=%.4e QBER=%6.2f%% "
                          "h=%.3f rate=%9d/s%s"
                          % ("OK " if row["capture_valid"] else "BAD",
                             row["version"], row["block_id"], trials, clicks,
                             errors, g, 100 * q, row["h_mean"] / 256.0,
                             row["attempts_per_sec"], extra))
                    for m in errs:
                        print("    ERROR " + m)
                    for m in warns:
                        print("    warn  " + m)

        except KeyboardInterrupt:
            print("\nstopped by user")

    print("\npackets=%d valid=%d (%.1f%%)  saved %s"
          % (count, good, 100.0 * good / max(count, 1), args.out))
    if good < count:
        print("Investigate validation_error before using this file.")


if __name__ == "__main__":
    main()
