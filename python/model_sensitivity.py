#!/usr/bin/env python3
"""
model_sensitivity.py -- numbers behind the precision and sensitivity tables.

  (1) fixed-point precision: 8-bit vs 12-bit fading word, 16-bit vs 24-bit
      probabilities, 256 / 1024 / 4096-entry inverse CDF vs continuous.
  (2) error of the first-order QBER expression (7) over the operating grid.
  (3) the range offsets actually stored in hooff_rom.
  (4) sensitivity of the reduced defect to the polarization constant k_s.
  (5) matched-variance family pairs and the background dose-response, as
      predictions for the rev B campaign.

ROM-exact wherever the hardware is concerned (rom_exact_model.py); the
continuous columns are labelled as such. No board is needed.

    python model_sensitivity.py
"""

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gamma as gamma_dist, norm, weibull_min

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rom_exact_model as rm                                   # noqa: E402
from uwoc_channel_model import lognormal_params, weibull_params  # noqa: E402

OUT_DIR = Path("phase2_results")
WATER_NAMES = ("clear ocean", "coastal", "harbor")
D_GRID = {
    0: (5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80),
    1: (2, 4, 6, 8, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43),
    2: (0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10),
}
SIGMA2 = {"L3": 0.2999, "L4": 1.0006, "L5": 2.9996}


def ppf(family, s2, u):
    if family == "gamma":
        return gamma_dist.ppf(u, a=1.0 / s2, scale=s2)
    if family == "lognormal":
        mu_x, sx2 = lognormal_params(s2)
        return np.exp(2.0 * norm.ppf(u, loc=mu_x, scale=math.sqrt(sx2)))
    b1, b2 = weibull_params(s2)
    return weibull_min.ppf(u, c=b1, scale=b2)


def moments(x):
    x = np.asarray(x, dtype=float)
    return x.mean(), x.var() / x.mean() ** 2


def precision_h():
    rows = []
    u = (np.arange(256) + 0.5) / 256
    for name, fam, s2 in (("h_s Gamma", "gamma", 4.5), ("L5 Weibull", "weibull", 2.9996),
                          ("L4 Weibull", "weibull", 1.0006), ("L3 log-normal", "lognormal", 0.2999)):
        x = ppf(fam, s2, u)
        h8 = np.clip(np.round(x * 128), 0, 255) / 128.0
        h12 = np.clip(np.round(x * 256), 0, 4095) / 256.0
        big = ppf(fam, s2, (np.arange(200000) + 0.5) / 200000)
        for label, v in (("continuous", big), ("8-bit (128=1, h<2)", h8),
                         ("12-bit (256=1, h<16)", h12)):
            m, s = moments(v)
            rows.append({"distribution": name, "sigma2_nominal": s2, "word": label,
                         "E_h": m, "sigma2": s})
    return pd.DataFrame(rows)


def precision_p(rom):
    rows = []
    for lam, lam_nm in enumerate((450, 532, 650)):
        for w in range(3):
            for di, d in enumerate(D_GRID[w]):
                code = rom.psig_ref(lam, w, di)
                p = code / rm.S24
                rows.append({"lambda_nm": lam_nm, "water": WATER_NAMES[w], "d_m": d,
                             "P_sig": p, "code24": code,
                             "code16": int(round(p * 65536)), "zero_at_16bit": int(round(p * 65536) == 0),
                             "zero_at_24bit": int(code == 0)})
    return pd.DataFrame(rows)


def lut_depth(rom):
    rows = []
    psig_ref = rom.psig_ref(0, 0, 4)
    epol = rom.epol_code(0, 0, 4)
    fams = {"L3": "lognormal", "L4": "weibull", "L5": "weibull"}
    for lvl, s2 in SIGMA2.items():
        for depth in (256, 1024, 4096, "continuous"):
            if depth == "continuous":
                n = 400000
                h = ppf(fams[lvl], s2, (np.arange(n) + 0.5) / n)
                ps = np.minimum(psig_ref * h, rm.P_MAX)          # no word truncation
                values = h
            else:
                u = (np.arange(depth) + 0.5) / depth
                hq = np.clip(np.round(ppf(fams[lvl], s2, u) * 256), 0, 4095).astype(np.int64)
                ps = np.minimum((psig_ref * hq) >> 8, rm.P_MAX)
                values = hq
            click, err = rm.event_probs(ps, rm.Y0_CODE, epol)
            G, T = 0.5 * click, 0.5 * err
            dfull, V0, G0, _ = rm.defect(G, T, np.ones(len(G)))
            Gb, Tb, cnt = rm.binned_moments(values, G, T, 8)
            d8 = rm.defect(Gb, Tb, cnt)[0]
            rows.append({"level": lvl, "lut": depth, "d8_over_V": d8 / V0,
                         "dfull_over_V": dfull / V0})
    return pd.DataFrame(rows)


def eq7_error(rom):
    """Relative error of eq. (7), which keeps Y0/2 undamped, against the exact
    emulator rule [e_pol P_s + (1 - P_s) Y0/2] / P_click."""
    y0 = rm.Y0_CODE / rm.S24
    worst = []
    for lam in range(3):
        for w in range(3):
            for di in range(16):
                p_nom = rom.psig_ref(lam, w, di) / rm.S24
                if p_nom < 1e-4:
                    continue
                e = rom.epol_code(lam, w, di) / rm.S24
                ps = np.minimum(p_nom * np.linspace(0.0, 16.0, 4001), 1.0)
                pc = 1.0 - (1.0 - ps) * (1.0 - y0)
                exact = (e * ps + (1.0 - ps) * y0 / 2.0) / pc
                eq7 = (e * ps + y0 / 2.0) / pc
                worst.append({"lambda_idx": lam, "water": WATER_NAMES[w],
                              "d_m": D_GRID[w][di],
                              "max_rel_error": float(np.max(np.abs(eq7 - exact) / exact))})
    return pd.DataFrame(worst)


def offsets(rom):
    rows = []
    for w in range(3):
        for off in sorted(set(int(v) for v in rom.hooff[w])):
            ds = [D_GRID[w][i] for i in range(16) if int(rom.hooff[w, i]) == off]
            rows.append({"water": WATER_NAMES[w], "offset": off,
                         "ranges_m": "%g-%g" % (min(ds), max(ds)) if len(ds) > 1 else "%g" % ds[0],
                         "L3_command_runs_row": min(5, max(1, 3 + off))})
    return pd.DataFrame(rows)


def ks_sensitivity(rom):
    rows = []
    b450 = 0.037 * 0.85
    for ks in (0.02, 0.04, 0.08):
        e = min(0.01 + ks * (1.0 - math.exp(-b450 * 25.0)), 0.5)
        code = int(round(e * rm.S24))
        for turb in (1, 3, 4, 5, 6, 7):
            hf = rm.fading_values(rom, 0, 4, turb)
            ps = rm.psig_codes(rom, 0, 0, 4, hf)
            click, err = rm.event_probs(ps, rm.Y0_CODE, code)
            G, T = 0.5 * click, 0.5 * err
            dfull, V0, G0, T0 = rm.defect(G, T, np.ones(len(G)))
            Gb, Tb, cnt = rm.binned_moments(hf, G, T, 8)
            rows.append({"k_s": ks, "e_pol": e, "row": turb,
                         "d8_over_V": rm.defect(Gb, Tb, cnt)[0] / V0,
                         "QBER": T0 / G0})
    return pd.DataFrame(rows)


def dose_response(rom):
    rows = []
    for code in (0, 54, 109, 218, 436, 872, 1744, 3488):
        for turb in (1, 2, 3, 4, 5, 6, 7):
            r = rm.point_defects(rom, turb=turb, p_noise=code)
            rows.append({"p_noise_code": code, "Y0": code / rm.S24, "row": turb,
                         "d8_over_V": r["d_bins_over_V"], "dfull_over_V": r["d_full_over_V"],
                         "bkm_ratio_8": r["bkm_ratio_bins"], "QBER": r["QBER"]})
    return pd.DataFrame(rows)


def main():
    OUT_DIR.mkdir(exist_ok=True)
    rom = rm.load_rom()
    fmt = lambda v: "%.4g" % v
    with pd.option_context("display.width", 220, "display.max_rows", 400):
        t = precision_h()
        t.to_csv(OUT_DIR / "sensitivity_precision_h.csv", index=False)
        print("== (1a) fading word ==\n" + t.to_string(index=False, float_format=fmt))
        t = precision_p(rom)
        t.to_csv(OUT_DIR / "sensitivity_precision_p.csv", index=False)
        print("\n== (1b) probability word: entries rounding to zero ==")
        print(t.groupby(["lambda_nm", "water"])[["zero_at_16bit", "zero_at_24bit"]].sum().to_string())
        print("smallest nonzero 24-bit P_sig: %.3e" % t.loc[t.code24 > 0, "P_sig"].min())
        t = lut_depth(rom)
        t.to_csv(OUT_DIR / "sensitivity_lut_depth.csv", index=False)
        print("\n== (1c) inverse-CDF depth, 25 m clear ocean ==\n"
              + t.pivot(index="lut", columns="level", values=["d8_over_V", "dfull_over_V"]).to_string(float_format=fmt))
        t = eq7_error(rom)
        t.to_csv(OUT_DIR / "eq7_error.csv", index=False)
        print("\n== (2) eq. (7) worst relative error over %d live points: %.2e ==" % (len(t), t.max_rel_error.max()))
        t = offsets(rom)
        t.to_csv(OUT_DIR / "range_offsets.csv", index=False)
        print("\n== (3) hooff_rom ==\n" + t.to_string(index=False))
        t = ks_sensitivity(rom)
        t.to_csv(OUT_DIR / "sensitivity_ks.csv", index=False)
        print("\n== (4) k_s ==\n" + t.pivot(index="row", columns="k_s", values="d8_over_V").to_string(float_format=fmt))
        t = dose_response(rom)
        t.to_csv(OUT_DIR / "dose_response_model.csv", index=False)
        print("\n== (5) background dose-response (8 bins) ==\n"
              + t.pivot(index="p_noise_code", columns="row", values="d8_over_V").to_string(float_format=fmt))


if __name__ == "__main__":
    main()
