#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uwoc_channel_model.py — Underwater wireless optical channel (UWOC) model for BB84
═══════════════════════════════════════════════════════════════════════════════

Replaces the atmospheric Gamma-Gamma model (gamma_gamma_final.v) with an
underwater composite-fading model:

        h = h_ℓ(d, λ; w) · h_f ,        h_f = h_s · h_o

    h_ℓ  : deterministic transmittance (Beer-Lambert + geometry), water type w
    h_s  : SCATTERING fading  ~ Gamma(1/σ_s², σ_s²)      [absent in the atmosphere]
    h_o  : OCEAN TURBULENCE fading ~ Lognormal (weak) / Weibull (moderate-strong)
    h_f  : the two fading terms alone — E[h_s] = E[h_o] = E[h_f] = 1 ⟹ E[h] = h_ℓ

⚠ NOTATION — h CONTAINS h_ℓ, h_f DOES NOT.  Every function below that takes a
  fading argument takes h_f; writing n̄ = μ·h_ℓ·h·η_det counts the path term twice.
  (w is the `water` argument, λ travels inside cfg.lam_nm — see h_ell().)

On top of that sits the PHOTON DETECTION STAGE (mandatory for UWOC, absent from
the current atmospheric model):

        n̄       = μ · h · η_det  =  μ · h_ℓ(d, λ; w) · h_f · η_det
        P_click = 1 − (1 − Y₀)·exp(−n̄)
        QBER    = [ e_pol·(1 − exp(−n̄)) + ½·Y₀ ] / P_click

REFERENCES
  [1] B. Kebapci et al., "FPGA-Based Implementation of an Underwater QKD System
      With BB84 Protocol," IEEE Photonics J., 15(4), 2023.
  [2] P. Salcedo-Serrano et al., "Ergodic capacity analysis of underwater FSO
      systems over scattering-induced fading channels in the presence of
      Weibull oceanic turbulence," IEEE ICC 2022, pp. 3814–3819.
      → eq.(2) h = h_ℓ·h_s·h_o, eq.(3) h_ℓ, eq.(4a/4b) Gamma/Weibull,
        eq.(5) σ²_ho, eq.(6) Nikishov spectrum, Table I/II/IV.
  [3] M. V. Jamali, F. Akhoundi, J. A. Salehi, "Performance Characterization of
      Relay-Assisted Wireless Optical CDMA Networks in Turbulent Underwater
      Channel," IEEE Trans. Wireless Commun., 15(6), 2016. → eq.(3),(4)
  [4] N. G. Nikishov & V. N. Nikishov, "Spectrum of turbulent fluctuations of
      the sea-water refraction index," Int. J. Fluid Mech. Res., 27, 2000.
  [5] L. C. Andrews & R. L. Phillips, "Laser Beam Propagation through Random
      Media," 2nd ed., SPIE Press, 2005.

Run the self-test / numerical checks:
    python python/uwoc_channel_model.py
    python python/uwoc_channel_model.py --plot
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq
from scipy.special import gamma as gamma_fn

__all__ = [
    "WATER_TYPES", "OCEAN_TURB_LEVELS", "LinkConfig",
    "nikishov_spectrum", "sigma2_ho", "sigma2_s",
    "lognormal_params", "weibull_params",
    "h_ell", "extinction", "e_pol",
    "n_bar", "p_click", "qber_from_nbar", "skr",
    "sample_h_s", "sample_h_o", "UWOCChannel",
]

# ═══════════════════════════════════════════════════════════════════════════
# 1. PHYSICAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════

# Nikishov spectrum — constants of δ and the exponential terms  [2] eq.(6), [3] eq.(4)
A_T = 1.863e-2          # decay coefficient of the TEMPERATURE term
A_S = 1.9e-4            # decay coefficient of the SALINITY term (100× smaller → long tail)
A_TS = 9.41e-3          # coefficient of the temperature-salinity cross term
ETA_KOLMOGOROV = 1e-3   # Kolmogorov microscale η [m]  (Table II of [2])


@dataclass(frozen=True)
class WaterType:
    """Water type + optical coefficients at the 532 nm reference wavelength."""
    name: str
    a: float            # absorption coefficient  a(λ) [m⁻¹]
    b: float            # scattering coefficient   b(λ) [m⁻¹]
    d_grid: tuple       # range grid used for the hardware LUT [m]
    source: str

    @property
    def c(self) -> float:
        """Total extinction coefficient c(λ) = a(λ) + b(λ)."""
        return self.a + self.b


# Table I of [2]; the 'harbor' row holds the standard Petzold values widely used
# in the UWOC literature (e.g. [3]) to extend the test range to very turbid water.
WATER_TYPES: dict[str, WaterType] = {
    "clear_ocean": WaterType(
        "Clear ocean", a=0.114, b=0.037,
        d_grid=(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80),
        source="Table I, Salcedo-Serrano ICC'22"),
    "coastal": WaterType(
        "Coastal", a=0.179, b=0.219,
        d_grid=(2, 4, 6, 8, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43),
        source="Table I, Salcedo-Serrano ICC'22"),
    "harbor": WaterType(
        "Harbor (turbid)", a=0.366, b=1.824,
        d_grid=(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10),
        source="Petzold / Jamali TWC'16"),
}

# Multiplier on c(λ) when changing wavelength, relative to the table above (referenced to 532 nm).
# THESE ARE INDICATIVE VALUES, not measurements: they serve to model the
# qualitative conclusion stated in [1] §I and Rosenkrantz & Arnon (SPIE 9224):
#   "blue-green outperform red-yellow-green at open ocean, while
#    red-yellow-green outperform blue-green at coastal turbid waters"
# → if accurate numbers are needed, substitute measured a(λ), b(λ) spectra for the water.
C_LAMBDA_SCALE: dict[str, dict[int, float]] = {
    "clear_ocean": {450: 0.85, 532: 1.00, 650: 2.60},
    "coastal":     {450: 1.10, 532: 1.00, 650: 1.45},
    "harbor":      {450: 1.35, 532: 1.00, 650: 0.90},
}
LAMBDA_OPTIONS = (450, 532, 650)     # nm — 3 options for the λ adaptation knob

# 5 ocean turbulence regimes. The (ε, χ_T, w) triples lie WITHIN the experimental
# Table II [2]: ε ∈ [1e-8, 1e-2] m²/s³, χ_T ∈ [1e-10, 1e-4] K²/s, w ∈ [-5, 0].
#   SMALL ε  → STRONG turbulence
#   LARGE χ_T → STRONG turbulence
#   w → 0   → salinity-dominated (WORST CASE);  w → -5 → temperature-dominated
#
# CALIBRATION: for each level, (ε, w) are spread evenly across the physical range,
# then χ_T is solved linearly (σ²_ho ∝ χ_T at fixed ε, w) so that σ²_ho at the
# 20 m reference range, λ=450 nm, lands exactly on:
#     L1 = 0.02   L2 = 0.08   L3 = 0.30   L4 = 1.00   L5 = 3.00
# The L4 mark of 1.0 is precisely the weak/strong turbulence boundary of [3].
# (To reproduce: see calibrate_turbulence_levels() below.)
OCEAN_TURB_LEVELS: dict[int, dict] = {
    0: dict(name="OFF",      eps=1e-2, chi_T=0.0,      w=-5.0),
    1: dict(name="VERYWEAK", eps=1e-2, chi_T=2.21e-07, w=-5.0),
    2: dict(name="WEAK",     eps=1e-3, chi_T=3.65e-07, w=-4.0),
    3: dict(name="MODERATE", eps=1e-4, chi_T=5.29e-07, w=-3.0),
    4: dict(name="STRONG",   eps=1e-5, chi_T=5.91e-07, w=-2.0),
    5: dict(name="SEVERE",   eps=1e-6, chi_T=3.85e-07, w=-1.0),
}

# Clamp ceiling for σ²_ho.
# REASON: eq.(5) is a WEAK-fluctuation result (Rytov theory, [5] ch.8) — valid
# only when σ² ≲ 1. At long range / strong turbulence the integral yields σ² ≫ 1
# (observed: up to 10⁴ at 70 m), which signals LEAVING THE DOMAIN OF VALIDITY
# rather than real physics — in practice σ²_I saturates around 1–2 through
# saturation of scintillation. We clamp at 5.0 and warn, so the Weibull fit stays stable.
SIGMA2_HO_MAX = 5.0


@dataclass
class LinkConfig:
    """Optoelectronic parameters of the link. Defaults follow the PoC of [1]."""
    lam_nm: int = 450           # wavelength [nm] — [1] uses 405 nm (blue)
    D_rx: float = 0.0508        # receiver aperture diameter [m] (2 inch)
    theta_div: float = 1e-3     # effective half divergence angle [rad].
                                # [1] uses a 1→7.2 mm beam expander; diffraction
                                # divergence is ~4e-5 rad, so the 1e-3 here is mostly
                                # margin for underwater pointing error.
    F: float = 0.85             # correction factor for scattered photons still collected, F ≤ 1
    eta_det: float = 0.18       # SPD efficiency @405 nm — [1] §II
    dark_hz: float = 60.0       # dark count [Hz] — [1] §II: "less than 60 Hz"
    bg_hz: float = 200.0        # background noise (ambient light) [Hz]
    gate_ns: float = 50.0       # detection time window / slot [ns]
    f_rep: float = 10e6         # pulse repetition frequency [Hz] — [1] §III-A: "20 ns
                                # pulses with the repetition rate of 10 MHz"
    tau_coh_ms: float = 5.0     # fading COHERENCE TIME underwater [ms].
                                # Ocean turbulence varies on a ms scale
                                # (~10⁴× slower than the qubit period) → h is
                                # "frozen" across tens of thousands of pulses.
    mu: float = 0.1             # mean photons / pulse (weak coherent)
    e0: float = 0.01            # intrinsic optical error (misalignment, PBS extinction)
    k_s: float = 0.04           # polarisation distortion from scattering — CALIBRATION PARAMETER
                                # (see e_pol(); calibrate against real measurements)
    wave_type: str = "plane"    # 'plane' (Θ=1) or 'spherical' (Θ=0)

    @property
    def k(self) -> float:
        """Wave number k = 2π/λ [rad/m]."""
        return 2.0 * math.pi / (self.lam_nm * 1e-9)

    @property
    def Y0(self) -> float:
        """Probability of a background click (dark + background) within one gate window."""
        return (self.dark_hz + self.bg_hz) * self.gate_ns * 1e-9

    @property
    def coherence_pulses(self) -> int:
        """Pulses fitting entirely inside one coherence block (fading treated as constant)."""
        return max(1, int(self.f_rep * self.tau_coh_ms * 1e-3))


# ═══════════════════════════════════════════════════════════════════════════
# 2. NIKISHOV SPECTRUM & OCEANIC SCINTILLATION INDEX σ²_ho
# ═══════════════════════════════════════════════════════════════════════════

def nikishov_spectrum(kappa, eps: float, chi_T: float, w: float,
                      eta: float = ETA_KOLMOGOROV):
    """
    Refractive-index turbulence power spectrum of seawater Φ_n(κ)  — [2] eq.(6), [4].

        Φ_n(κ) = 0.388e-8 · ε^(-1/3) · κ^(-11/3) · [1 + 2.35(κη)^(2/3)]
                 · (χ_T / w²) · ( w²·e^(-A_T δ) + e^(-A_S δ) - 2w·e^(-A_TS δ) )
        δ = 8.284(κη)^(4/3) + 12.978(κη)²

    kappa : spatial frequency [rad/m], scalar or ndarray
    eps   : turbulent kinetic energy dissipation rate [m²/s³]
    chi_T : temperature variance dissipation rate [K²/s]
    w     : temperature / salinity balance ratio, ∈ [-5, 0)
    """
    kappa = np.asarray(kappa, dtype=float)
    ke = kappa * eta
    delta = 8.284 * ke ** (4.0 / 3.0) + 12.978 * ke ** 2

    # Three exponential terms. The salinity term (A_S) decays ~100× more slowly than
    # the temperature term (A_T) → the spectral tail at large κ is salinity-dominated.
    # That is exactly why w → 0 (salinity-dominated) gives the strongest turbulence.
    term = (w ** 2 * np.exp(-A_T * delta)
            + np.exp(-A_S * delta)
            - 2.0 * w * np.exp(-A_TS * delta))

    return (0.388e-8 * eps ** (-1.0 / 3.0) * kappa ** (-11.0 / 3.0)
            * (1.0 + 2.35 * ke ** (2.0 / 3.0)) * (chi_T / w ** 2) * term)


def sigma2_ho(d: float, eps: float, chi_T: float, w: float,
              lam_nm: float = 450, wave_type: str = "plane",
              n_kappa: int = 4000, n_xi: int = 48) -> float:
    """
    Scintillation index from ocean turbulence  — [2] eq.(5), [3] eq.(3), [5].

        σ²_ho = 8π²k²d ∫₀¹∫₀^∞ κ Φ_n(κ)·[1 − cos( dκ²ξ(1−(1−Θ)ξ)/k )] dκ dξ

    Θ = 1 (plane wave) or Θ = 0 (spherical wave).

    NO closed form → evaluated numerically here, the result going into the FPGA ROM.

    Numerical technique:
      • Integrate κ in log space (κ Φ_n ~ κ^(-8/3) as κ→0, but [1−cos] ~ κ⁴
        so the integral converges; at large κ the exponential terms kill it).
      • The κ range [1e-4, 1e6] rad/m fully covers the salinity tail (A_S·δ ≳ 20).
      • When the cos argument exceeds 50 rad the oscillation is too fast → the
        mean of cos ≈ 0, so [1−cos] is replaced by 1 to avoid sampling error (standard).
    """
    if chi_T <= 0.0 or d <= 0.0:
        return 0.0

    k = 2.0 * math.pi / (lam_nm * 1e-9)
    theta = 1.0 if wave_type == "plane" else 0.0

    # Log-spaced κ grid and integration weights in ln κ  (dκ = κ d(lnκ))
    ln_k = np.linspace(math.log(1e-4), math.log(1e6), n_kappa)
    kap = np.exp(ln_k)
    dlnk = ln_k[1] - ln_k[0]

    phi = nikishov_spectrum(kap, eps, chi_T, w)          # (n_kappa,)
    base = kap ** 2 * phi                                # κ·Φ_n·κ  (change of variable)

    # Gauss-Legendre over ξ ∈ [0, 1]
    xi_nodes, xi_w = np.polynomial.legendre.leggauss(n_xi)
    xi = 0.5 * (xi_nodes + 1.0)
    xi_w = 0.5 * xi_w

    # arg[i, j] = d·κ_j²·ξ_i·(1 − (1−Θ)ξ_i) / k
    shape = xi * (1.0 - (1.0 - theta) * xi)              # (n_xi,)
    arg = (d / k) * np.outer(shape, kap ** 2)            # (n_xi, n_kappa)

    osc = np.where(np.abs(arg) > 50.0, 1.0, 1.0 - np.cos(arg))
    integ = (osc * base[None, :]).sum(axis=1) * dlnk     # κ integral for each ξ
    val = float((integ * xi_w).sum())                    # ξ integral

    return 8.0 * math.pi ** 2 * k ** 2 * d * val


def sigma2_ho_clamped(d: float, eps: float, chi_T: float, w: float,
                      lam_nm: float = 450, wave_type: str = "plane") -> float:
    """sigma2_ho() with the SIGMA2_HO_MAX ceiling applied (see the note at the top of the file)."""
    return min(sigma2_ho(d, eps, chi_T, w, lam_nm, wave_type), SIGMA2_HO_MAX)


def calibrate_turbulence_levels(targets=(0.02, 0.08, 0.30, 1.00, 3.00),
                                seeds=((1e-2, -5.0), (1e-3, -4.0), (1e-4, -3.0),
                                       (1e-5, -2.0), (1e-6, -1.0)),
                                d_ref: float = 20.0, lam_nm: int = 450):
    """
    Reproduce the OCEAN_TURB_LEVELS table: for each given (ε, w), solve for χ_T so
    that σ²_ho(d_ref) equals the target exactly. Valid because σ²_ho ∝ χ_T linearly.
    """
    out = {}
    for lv, ((eps, w), tgt) in enumerate(zip(seeds, targets), start=1):
        chi0 = 1e-8
        s0 = sigma2_ho(d_ref, eps, chi0, w, lam_nm)
        out[lv] = dict(eps=eps, w=w, chi_T=chi0 * tgt / s0)
    return out


# ═══════════════════════════════════════════════════════════════════════════
# 3. FADING DISTRIBUTIONS
# ═══════════════════════════════════════════════════════════════════════════

def lognormal_params(sigma2: float) -> tuple[float, float]:
    """
    Lognormal parameters for WEAK turbulence (σ² < 1)  — [3] eq.(1).

        f(h) = 1/(2h√(2πσ_X²)) · exp( −(ln h − 2μ_X)² / (8σ_X²) )

    with σ_X² = the log-amplitude variance. Relations:  σ²_I = exp(4σ_X²) − 1,
    and normalising E[h] = 1  ⟹  μ_X = −σ_X².
    Returns (mu_X, sigma_X2).
    """
    sigma_X2 = 0.25 * math.log1p(max(sigma2, 1e-12))
    return -sigma_X2, sigma_X2


def weibull_params(sigma2: float) -> tuple[float, float]:
    """
    Weibull parameters for MODERATE-STRONG turbulence  — [2] eq.(4b).

        f(h) = (β₁/β₂)·(h/β₂)^(β₁−1) · exp( −(h/β₂)^β₁ )

    Solve β₁ from    σ² = Γ(1+2/β₁)/Γ²(1+1/β₁) − 1    (decreasing in β₁),
    then   β₂ = 1/Γ(1+1/β₁)   to normalise E[h] = 1.
    """
    def scint(b1: float) -> float:
        return gamma_fn(1.0 + 2.0 / b1) / gamma_fn(1.0 + 1.0 / b1) ** 2 - 1.0

    sigma2 = float(np.clip(sigma2, 1e-6, 50.0))
    b1 = brentq(lambda x: scint(x) - sigma2, 0.12, 200.0, xtol=1e-10)
    b2 = 1.0 / gamma_fn(1.0 + 1.0 / b1)
    return b1, b2


# ---- Scattering fading: σ_s²(d) -----------------------------------------
# Table IV of [2] (at w = −3). The paper calls σ_s² "exponentially distance-
# dependent"; we fit ln σ_s² = B·(d − d₁) with d₁ = the range where σ_s² = 1.
TABLE_IV = {
    "clear_ocean": [(55, 0.153), (65, 1.157), (67, 1.758), (70, 3.291)],
    "coastal":     [(27, 0.142), (33, 0.916), (35, 1.705), (37, 3.172)],
}


def _fit_sigma_s() -> dict[str, tuple[float, float]]:
    """Fit (B, d₁) for each water type that has data; extrapolate for harbor."""
    fits: dict[str, tuple[float, float]] = {}
    for wt, pts in TABLE_IV.items():
        d = np.array([p[0] for p in pts], dtype=float)
        y = np.log(np.array([p[1] for p in pts], dtype=float))
        B, intercept = np.polyfit(d, y, 1)      # ln σ_s² = B·d + intercept
        d1 = -intercept / B                     # where ln σ_s² = 0
        fits[wt] = (float(B), float(d1))

    # Harbor extrapolation: power law in c(λ) through 2 known anchor points.
    #   B(c)  = B₀ · c^s ,   d₁(c) = d₁₀ · c^t
    c_lo, c_hi = WATER_TYPES["clear_ocean"].c, WATER_TYPES["coastal"].c
    (B_lo, d1_lo), (B_hi, d1_hi) = fits["clear_ocean"], fits["coastal"]
    s = math.log(B_hi / B_lo) / math.log(c_hi / c_lo)
    t = math.log(d1_hi / d1_lo) / math.log(c_hi / c_lo)
    c_h = WATER_TYPES["harbor"].c
    fits["harbor"] = (B_lo * (c_h / c_lo) ** s, d1_lo * (c_h / c_lo) ** t)
    return fits


SIGMA_S_FIT = _fit_sigma_s()


def sigma2_s(d: float, water: str) -> float:
    """
    Scattering fading strength σ_s² at range d [m]  — calibrated against Table IV [2].

    ⚠ For 'harbor' this is an EXTRAPOLATION (power law in c), not measured
      data — [2] only publishes clear ocean and coastal.
    """
    B, d1 = SIGMA_S_FIT[water]
    return float(np.clip(math.exp(B * (d - d1)), 1e-4, 50.0))


def sample_h_s(sigma2: float, size, rng: np.random.Generator) -> np.ndarray:
    """h_s ~ Gamma(shape=1/σ_s², scale=σ_s²), E[h_s] = 1  — [2] eq.(4a)."""
    if sigma2 <= 1e-6:
        return np.ones(size)
    return rng.gamma(shape=1.0 / sigma2, scale=sigma2, size=size)


def sample_h_o(sigma2: float, size, rng: np.random.Generator) -> np.ndarray:
    """h_o ~ Lognormal (σ² < 1) or Weibull (σ² ≥ 1), E[h_o] = 1."""
    if sigma2 <= 1e-6:
        return np.ones(size)
    if sigma2 < 1.0:
        mu_X, sX2 = lognormal_params(sigma2)
        # h = exp(2X), X ~ N(mu_X, sX2)  → per the log-amplitude convention of [3]
        return np.exp(2.0 * rng.normal(mu_X, math.sqrt(sX2), size=size))
    b1, b2 = weibull_params(sigma2)
    return b2 * rng.weibull(b1, size=size)


# ═══════════════════════════════════════════════════════════════════════════
# 4. LOSS, POLARISATION ERROR, DETECTION STAGE
# ═══════════════════════════════════════════════════════════════════════════

def extinction(water: str, lam_nm: int = 532) -> float:
    """c(λ) = a + b, scaled with wavelength (see the C_LAMBDA_SCALE note)."""
    wt = WATER_TYPES[water]
    nearest = min(C_LAMBDA_SCALE[water], key=lambda L: abs(L - lam_nm))
    return wt.c * C_LAMBDA_SCALE[water][nearest]


def scattering_coef(water: str, lam_nm: int = 532) -> float:
    """b(λ) — used for the scattering optical depth τ_s = b·d."""
    wt = WATER_TYPES[water]
    nearest = min(C_LAMBDA_SCALE[water], key=lambda L: abs(L - lam_nm))
    return wt.b * C_LAMBDA_SCALE[water][nearest]


def h_ell(d: float, water: str, cfg: LinkConfig) -> float:
    """
    Deterministic transmittance h_ℓ — [2] eq.(3), kept in the paper's form
    (including the π):

        h_ℓ(d, λ; w) = D_rx² / (π·(d·tanθ_div)²) · exp( −F·c_w(λ)·d )
                        └─── geometry ────┘         └─ Beer-Lambert ─┘

    Arguments vs parameter: d and λ (= cfg.lam_nm) are what a sweep varies, the
    water type w (= `water`) is a categorical scenario label — hence the ";" in
    the symbol.  NOTE that w and λ enter ONLY through c_w(λ): the geometric term
    depends on neither.  That is why the wavelength knob reduces to exactly
    argmin_λ c_w(λ), which is the whole λ-selection claim.

    It is a TRANSMITTANCE (≤ 1, multiplies), not a loss — do not read the ℓ as
    "loss in dB", which would carry the opposite sign.

    The geometric term is clamped to ≤ 1 (you cannot collect more than you send).
    """
    geo = cfg.D_rx ** 2 / (math.pi * (d * math.tan(cfg.theta_div)) ** 2)
    geo = min(geo, 1.0)
    return geo * math.exp(-cfg.F * extinction(water, cfg.lam_nm) * d)


def e_pol(d: float, water: str, cfg: LinkConfig) -> float:
    """
    Polarisation error from MULTIPLE SCATTERING — a UWOC-specific mechanism, absent in air.

        e_det(d) = e₀ + k_s·( 1 − exp(−b(λ)·d) )

    Paper [1] attributes the QBER rise with distance to "polarization distortion
    and low SNR", and fits it LINEARLY in d (Fig. 19). For small b·d,
    1 − e^(−bd) ≈ b·d ⟹ this model is linear in that region, matching Fig. 19.
    k_s is a CALIBRATION PARAMETER to be fitted per system from measurements.
    Clamped at 0.5 (complete loss of polarisation information).
    """
    tau_s = scattering_coef(water, cfg.lam_nm) * d
    return float(min(0.5, cfg.e0 + cfg.k_s * (1.0 - math.exp(-tau_s))))


def n_bar(h_f: float | np.ndarray, d: float, water: str,
          cfg: LinkConfig) -> float | np.ndarray:
    """
    Mean photo-electrons per pulse:  n̄ = μ · h_ℓ(d, λ; w) · h_f · η_det.

    ⚠ TAKES h_f — the fading pair h_s·h_o with unit mean — NOT the composite h.
      h already contains h_ℓ, so passing it here applies the path term twice.
      Pass h_f = 1.0 to get the mean n̄ (since E[h_f] = 1).
    """
    return cfg.mu * h_ell(d, water, cfg) * h_f * cfg.eta_det


def p_click(nbar, Y0: float):
    """P_click = 1 − (1 − Y₀)·exp(−n̄)."""
    return 1.0 - (1.0 - Y0) * np.exp(-nbar)


def qber_from_nbar(nbar, Y0: float, e_det: float):
    """
        QBER = [ e_det·(1 − e^(−n̄)) + ½·Y₀ ] / P_click

    A signal click is wrong with probability e_det; a background / dark-count click
    is fully random and therefore wrong 50% of the time.
    """
    pc = p_click(nbar, Y0)
    num = e_det * (1.0 - np.exp(-nbar)) + 0.5 * Y0
    return np.where(pc > 0, num / np.maximum(pc, 1e-300), 0.5)


def h2(x):
    """Binary entropy H₂(x)."""
    x = np.clip(x, 1e-12, 1 - 1e-12)
    return -x * np.log2(x) - (1 - x) * np.log2(1 - x)


def skr(pclick: float, qber: float, cfg: LinkConfig, q: float = 0.5) -> float:
    """
    Asymptotic secret key rate (BB84, no decoy states):

        R = q · f_rep · P_click · max(0, 1 − 2·H₂(QBER))

    q = 1/2 for uniform bases; q = p_z² + p_x² for biased bases (the basis_prob_z knob).
    """
    return max(0.0, q * cfg.f_rep * pclick * (1.0 - 2.0 * h2(qber)))


# ═══════════════════════════════════════════════════════════════════════════
# 5. COMPOSITE CHANNEL CLASS (used for Monte-Carlo and by lut_gen)
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class UWOCChannel:
    """UWOC channel for one configuration (water type, range, turbulence level)."""
    water: str
    distance: float
    turb_level: int
    cfg: LinkConfig = field(default_factory=LinkConfig)

    def __post_init__(self):
        t = OCEAN_TURB_LEVELS[self.turb_level]
        self.sigma2_ho = sigma2_ho_clamped(self.distance, t["eps"], t["chi_T"],
                                           t["w"], self.cfg.lam_nm,
                                           self.cfg.wave_type)
        self.sigma2_s = sigma2_s(self.distance, self.water)
        self.h_ell = h_ell(self.distance, self.water, self.cfg)
        self.e_det = e_pol(self.distance, self.water, self.cfg)
        self.Y0 = self.cfg.Y0
        self.nbar_mean = n_bar(1.0, self.distance, self.water, self.cfg)

    # ---- mean statistics (no Monte-Carlo needed) ------------------------
    def mean_metrics(self) -> dict:
        pc = p_click(self.nbar_mean, self.Y0)
        q = float(qber_from_nbar(self.nbar_mean, self.Y0, self.e_det))
        return dict(sigma2_ho=self.sigma2_ho, sigma2_s=self.sigma2_s,
                    # key kept as L_dB: it is a COLUMN NAME of data/sim_table.csv,
                    # renaming it would invalidate the closed campaign dataset
                    L_dB=10 * math.log10(max(self.h_ell, 1e-300)),
                    nbar=self.nbar_mean, p_click=pc, qber=q,
                    e_det=self.e_det, skr=skr(pc, q, self.cfg))

    # ---- WINDOWED statistics (where turbulence actually shows up) --------
    def window_statistics(self, n_windows: int = 4000, window: int | None = None,
                          seed: int = 5, p_basis_z: float = 0.5,
                          qber_limit: float = 0.11) -> dict:
        """
        Per-sliding-window QBER / click statistics — the very quantities that
        channel_monitor.v and adaptive_controller.v actually see.

        ⚠ THE KEY PHYSICAL POINT — COHERENCE TIME:
          Ocean turbulence varies on a ms scale, while the qubit rate is
          10⁵–10⁶ /s.  So within ONE 256-qubit window (~256 µs @1 MHz) the fading
          coefficient h_f is effectively FROZEN, not resampled per qubit.

          Consequence: the LONG-TERM mean of QBER barely depends on the turbulence
          level (P_click ≈ linear in h_f, so E[·] cancels the fading).
          Turbulence only shows up through the VARIANCE BETWEEN WINDOWS and through
          the OUTAGE PROBABILITY.  That is the signal the adaptive controller exploits.
          → Therefore the RTL MUST sample h on a coherence-time counter, NOT
            per qubit (see uwoc_channel.v §3).

        window = None → use exactly one coherence block (cfg.coherence_pulses).
        That is the RIGHT window size: longer and the fading is averaged away,
        shorter and the QBER estimate drowns in shot noise.
        """
        if window is None:
            window = self.cfg.coherence_pulses
        rng = np.random.default_rng(seed)

        # 1 fading sample per window (fading frozen within the window)
        h_f = (sample_h_s(self.sigma2_s, n_windows, rng)
               * sample_h_o(self.sigma2_ho, n_windows, rng))
        nb = n_bar(h_f, self.distance, self.water, self.cfg)

        p_sig = 1.0 - np.exp(-nb)                      # (n_windows,)
        q_basis = p_basis_z ** 2 + (1 - p_basis_z) ** 2

        # Within each window: binomial counting instead of per-qubit simulation
        n_sift_att = rng.binomial(window, q_basis, size=n_windows)
        sig = rng.binomial(n_sift_att, np.clip(p_sig, 0, 1))
        rest = n_sift_att - sig
        noise = rng.binomial(rest, min(self.Y0, 1.0))
        clicks = sig + noise
        errs = (rng.binomial(sig, min(self.e_det, 1.0))
                + rng.binomial(noise, 0.5))

        valid = clicks > 0
        qw = np.where(valid, errs / np.maximum(clicks, 1), 0.5)

        outage = float((qw > qber_limit).mean())
        return dict(
            qber_win=qw, clicks_win=clicks, window=window,
            clicks_per_window=float(clicks.mean()),
            qber_mean=float(qw.mean()), qber_p90=float(np.percentile(qw, 90)),
            qber_std=float(qw.std()),
            click_rate=float(clicks.sum() / (n_windows * window)),
            outage=outage,
            # SKR accounts for outage: an insecure window is discarded entirely
            skr_per_pulse=float(np.mean(np.where(
                qw <= qber_limit,
                (clicks / window) * np.maximum(0.0, 1 - 2 * h2(qw)), 0.0))),
        )

    # ---- Qubit-level Monte-Carlo ----------------------------------------
    def simulate(self, n: int, seed: int = 42, p_basis_z: float = 0.5) -> dict:
        """
        Simulate n BB84 qubits through the channel. Returns statistics + per-qubit arrays.

        The decision structure is IDENTICAL to the RTL (uwoc_channel.v) for cross-checking:
          sig_det   ~ Bernoulli(1 − e^(−n̄))       → click from the signal
          noise_det ~ Bernoulli(Y₀)               → click do dark/background
          click     = sig_det OR noise_det
          error     = sig_det ? Bern(e_det) : Bern(0.5)
        This split is exactly equivalent to the QBER formula above but needs
        NO DIVISION → it is synthesisable on the FPGA.
        """
        rng = np.random.default_rng(seed)

        h_s = sample_h_s(self.sigma2_s, n, rng)
        h_o = sample_h_o(self.sigma2_ho, n, rng)
        nb = n_bar(h_s * h_o, self.distance, self.water, self.cfg)

        sig_det = rng.random(n) < (1.0 - np.exp(-nb))
        noise_det = rng.random(n) < self.Y0
        click = sig_det | noise_det

        err_p = np.where(sig_det, self.e_det, 0.5)
        error = rng.random(n) < err_p

        a_basis = rng.random(n) >= p_basis_z          # False = Z, True = X
        b_basis = rng.random(n) >= p_basis_z
        basis_match = (a_basis == b_basis)

        sifted = click & basis_match
        n_sift = int(sifted.sum())
        n_err = int((sifted & error).sum())
        q = n_err / n_sift if n_sift else 0.5
        pc = float(click.mean())
        q_eff = p_basis_z ** 2 + (1 - p_basis_z) ** 2

        return dict(
            n=n, h_s=h_s, h_o=h_o, nbar=nb, click=click, error=error,
            basis_match=basis_match, sifted=sifted,
            p_click=pc, sift_eff=n_sift / n, qber=q,
            skr_per_pulse=max(0.0, (n_sift / n) * (1 - 2 * h2(q))),
            skr_bps=skr(pc, q, self.cfg, q=q_eff),
        )


# ═══════════════════════════════════════════════════════════════════════════
# 6. SELF-TEST / NUMERICAL CHECKS
# ═══════════════════════════════════════════════════════════════════════════

def _hr(t=""):
    print("\n" + "═" * 78)
    if t:
        print("  " + t)
        print("═" * 78)


def verify_sigma_s_fit():
    """Cross-check the σ_s² fit against Table IV of [2]."""
    _hr("KIỂM CHỨNG 1 — fit σ_s²(d) so với Table IV (Salcedo-Serrano ICC'22)")
    print(f"{'Nước':<13}{'d [m]':>7}{'σ_s² (paper)':>15}{'σ_s² (fit)':>13}"
          f"{'sai số':>10}")
    print("-" * 78)
    worst = 0.0
    for wt, pts in TABLE_IV.items():
        B, d1 = SIGMA_S_FIT[wt]
        for d, ref in pts:
            got = sigma2_s(d, wt)
            e = abs(got - ref) / ref * 100
            worst = max(worst, e)
            print(f"{wt:<13}{d:>7.0f}{ref:>15.3f}{got:>13.3f}{e:>9.1f}%")
        print(f"{'':13}→ σ_s²(d) = exp({B:.4f}·(d − {d1:.2f}))")
    B, d1 = SIGMA_S_FIT["harbor"]
    print(f"{'harbor':<13}  (NGOẠI SUY, không có số liệu đo)"
          f"   σ_s²(d) = exp({B:.4f}·(d − {d1:.2f}))")
    print(f"\n  Sai số lớn nhất: {worst:.2f}%  "
          f"→ {'ĐẠT' if worst < 5 else 'CẦN XEM LẠI'} (ngưỡng 5%)")
    return worst < 5.0


def verify_weibull_inversion():
    """Check the σ² → (β₁, β₂) inversion and the E[h] = 1 normalisation."""
    _hr("KIỂM CHỨNG 2 — nghịch đảo tham số Weibull / Lognormal")
    print(f"{'σ² mục tiêu':>13}{'β₁':>10}{'β₂':>10}"
          f"{'E[h] MC':>11}{'σ² MC':>11}{'sai số σ²':>12}")
    print("-" * 78)
    rng = np.random.default_rng(7)
    ok = True
    for s2 in (0.05, 0.2, 0.5, 0.9, 1.0, 1.5, 2.5, 4.0):
        h = sample_h_o(s2, 400_000, rng)
        m, v = h.mean(), h.var() / h.mean() ** 2
        e = abs(v - s2) / s2 * 100
        ok &= (e < 3.0 and abs(m - 1) < 0.02)
        if s2 >= 1.0:
            b1, b2 = weibull_params(s2)
            tag = f"{b1:>10.3f}{b2:>10.3f}"
        else:
            _, sX2 = lognormal_params(s2)
            tag = f"{'LN':>10}{math.sqrt(sX2):>10.3f}"
        print(f"{s2:>13.2f}{tag}{m:>11.4f}{v:>11.4f}{e:>11.2f}%")
    print(f"\n  {'ĐẠT' if ok else 'CẦN XEM LẠI'} — E[h]=1 (±2%) và σ² khớp (±3%)")
    return ok


def verify_sigma_ho():
    """Print σ²_ho versus range for the 5 turbulence levels — check monotonicity."""
    _hr("KIỂM CHỨNG 3 — σ²_ho(d) từ phổ Nikishov (sóng phẳng, λ = 450 nm)")
    dists = [5, 10, 20, 30, 50, 70]
    print(f"{'Level':<11}{'ε [m²/s³]':>12}{'χ_T [K²/s]':>13}{'w':>5}"
          + "".join(f"{d:>9}m" for d in dists))
    print("-" * 78)
    prev = None
    mono = True
    for lv in range(1, 6):
        t = OCEAN_TURB_LEVELS[lv]
        vals = [sigma2_ho(d, t["eps"], t["chi_T"], t["w"], 450) for d in dists]
        print(f"{t['name']:<11}{t['eps']:>12.0e}{t['chi_T']:>13.0e}"
              f"{t['w']:>5.0f}" + "".join(f"{v:>10.3f}" for v in vals))
        if prev is not None:
            mono &= all(a <= b * 1.001 for a, b in zip(prev, vals))
        prev = vals
        # monotonically increasing with d within one level
        mono &= all(vals[i] <= vals[i + 1] * 1.001 for i in range(len(vals) - 1))
    print(f"\n  Đơn điệu tăng theo d và theo level: "
          f"{'ĐẠT' if mono else 'CẦN XEM LẠI'}")
    print("  Đối chiếu định tính [3]: nhiễu loạn MẠNH (σ² > 1) xuất hiện ở cự ly")
    print("  chỉ 10–100 m — khác hẳn khí quyển (hàng km). Kết quả trên phù hợp.")
    return mono


def verify_qber_vs_distance(cfg: LinkConfig):
    """Reproduce the shape of the QBER-vs-distance curve of [1] Fig. 19."""
    _hr("KIỂM CHỨNG 4 — QBER & SKR theo cự ly (tái tạo dạng Fig. 19 của [1])")
    for water in ("clear_ocean", "coastal", "harbor"):
        wt = WATER_TYPES[water]
        print(f"\n  ▸ {wt.name}  (a={wt.a}, b={wt.b}, c={wt.c:.3f} m⁻¹)  "
              f"turb level 3 (MODERATE)")
        print(f"    {'d [m]':>7}{'L [dB]':>10}{'σ²_ho':>9}{'σ_s²':>9}"
              f"{'e_pol':>9}{'n̄':>11}{'P_click':>10}{'QBER':>9}{'SKR [bps]':>12}")
        print("    " + "-" * 82)
        for d in wt.d_grid[::2]:
            ch = UWOCChannel(water, float(d), 3, cfg)
            m = ch.mean_metrics()
            flag = "" if m["qber"] < 0.11 else "  ✗>11%"
            print(f"    {d:>7.1f}{m['L_dB']:>10.1f}{m['sigma2_ho']:>9.3f}"
                  f"{m['sigma2_s']:>9.3f}{m['e_det']*100:>8.2f}%"
                  f"{m['nbar']:>11.2e}{m['p_click']:>10.2e}"
                  f"{m['qber']*100:>8.2f}%{m['skr']:>12.1f}{flag}")
    print("\n  Nhận xét kỳ vọng: QBER tăng gần TUYẾN TÍNH ở cự ly ngắn (vì")
    print("  1−e^(−bd) ≈ bd), rồi tăng vọt khi P_click chạm sàn dark count.")


def max_secure_distance(water: str, turb_level: int, cfg: LinkConfig,
                        qber_limit: float = 0.11, criterion: str = "mean",
                        outage_limit: float = 0.10,
                        window: int | None = None) -> float:
    """
    Largest range still secure (found by binary search).

    criterion = 'mean'   : long-term mean QBER < qber_limit
                           → NOT sensitive to the turbulence level (see window_statistics)
    criterion = 'outage' : P(window QBER > qber_limit) < outage_limit
                           → the practical criterion, and it DOES separate turbulence levels
    """
    def secure(d: float) -> bool:
        ch = UWOCChannel(water, d, turb_level, cfg)
        if criterion == "mean":
            return ch.mean_metrics()["qber"] < qber_limit
        st = ch.window_statistics(n_windows=1200, window=window,
                                  seed=int(d * 977) & 0xFFFF,
                                  qber_limit=qber_limit)
        return st["outage"] < outage_limit

    lo, hi = 0.2, WATER_TYPES[water].d_grid[-1] * 1.5
    if not secure(lo):
        return 0.0
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if secure(mid):
            lo = mid
        else:
            hi = mid
    return lo


def verify_wavelength_choice(cfg: LinkConfig):
    """Check the conclusion: the optimal λ depends on the water type."""
    _hr("KIỂM CHỨNG 5 — bước sóng tối ưu theo loại nước (knob thích nghi λ)")
    print("  So sánh tại cự ly bảo mật tối đa (d_max ứng với λ=532 nm),")
    print("  chỉ tiêu = SKR [bps]; ô có QBER ≥ 11% ghi 'insec'.\n")
    print(f"  {'Nước':<17}{'d [m]':>7}"
          + "".join(f"{L:>15}nm" for L in LAMBDA_OPTIONS) + f"{'λ tốt nhất':>14}")
    print("  " + "-" * 79)
    verdicts = {}
    for water in ("clear_ocean", "coastal", "harbor"):
        ref = LinkConfig(**{**cfg.__dict__, "lam_nm": 532})
        d = max(0.5, 0.7 * max_secure_distance(water, 3, ref))
        cells, skrs = [], []
        for lam in LAMBDA_OPTIONS:
            c2 = LinkConfig(**{**cfg.__dict__, "lam_nm": lam})
            m = UWOCChannel(water, d, 3, c2).mean_metrics()
            skrs.append(m["skr"])
            cells.append(f"{m['skr']:>17.1f}" if m["skr"] > 0 else f"{'insec':>17}")
        best = LAMBDA_OPTIONS[int(np.argmax(skrs))]
        verdicts[water] = best
        print(f"  {WATER_TYPES[water].name:<17}{d:>7.1f}" + "".join(cells)
              + f"{best:>12} nm")
    ok = (verdicts["clear_ocean"] <= verdicts["coastal"] <= verdicts["harbor"])
    print("\n  Kỳ vọng theo [1] §I: clear ocean ưu blue-green; nước đục dịch dần")
    print("  về phía red-yellow.  λ*(clear) ≤ λ*(coastal) ≤ λ*(harbor): "
          f"{'ĐẠT' if ok else 'CẦN XEM LẠI'}")
    return ok


def verify_monte_carlo(cfg: LinkConfig):
    """
    Cross-check the qubit-level Monte-Carlo against the analytical formulas.

    Split into 2 INDEPENDENT tests so the two error sources do not mix:
      (a) LINEARISATION ERROR (Jensen gap): compare E[1−e^(−n̄(h_f))] averaged over
          the distribution of h_f  with  1−e^(−n̄(E[h_f]))  at the centre. This is
          exactly the assumption the RTL relies on (P_click linear in h_f) → it
          must be small.
      (b) Bernoulli SAMPLING ERROR: run at an operating point with a large enough
          P_click (short range) to gather enough clicks to test for bias.
    """
    _hr("KIỂM CHỨNG 6 — Monte-Carlo vs giải tích")
    rng = np.random.default_rng(11)

    print("  (a) Sai số tuyến tính hoá P_click ≈ μ·h_ℓ·h_f·η  (giả thiết của RTL)")
    print(f"      {'Nước':<13}{'d':>6}{'lv':>4}{'n̄(E[h_f])':>12}"
          f"{'P_click tâm':>14}{'E[P_click]':>13}{'khe hở':>10}")
    print("      " + "-" * 72)
    ok_lin = True
    for water, d in (("clear_ocean", 20.0), ("coastal", 8.0), ("harbor", 1.5)):
        for lv in (1, 3, 5):
            ch = UWOCChannel(water, d, lv, cfg)
            h_f = (sample_h_s(ch.sigma2_s, 400_000, rng)
                   * sample_h_o(ch.sigma2_ho, 400_000, rng))
            nb = n_bar(h_f, d, water, cfg)
            exact = float(p_click(nb, ch.Y0).mean())
            centre = p_click(ch.nbar_mean, ch.Y0)
            gap = abs(exact - centre) / max(exact, 1e-300) * 100
            ok_lin &= gap < 5.0
            print(f"      {water:<13}{d:>6.1f}{lv:>4}{ch.nbar_mean:>12.2e}"
                  f"{centre:>14.4e}{exact:>13.4e}{gap:>9.2f}%")

    print("\n  (b) Tính không chệch của tầng quyết định Bernoulli (RTL-style)")
    print(f"      {'Nước':<13}{'d':>6}{'lv':>4}{'P_click lt':>13}{'P_click MC':>13}"
          f"{'QBER lt':>10}{'QBER MC':>10}{'#sift':>8}")
    print("      " + "-" * 72)
    ok_mc = True
    for water, d in (("clear_ocean", 8.0), ("coastal", 3.0), ("harbor", 0.5)):
        for lv in (1, 3, 5):
            ch = UWOCChannel(water, d, lv, cfg)
            an = ch.mean_metrics()
            mc = ch.simulate(300_000, seed=100 + lv)
            n_sift = int(mc["sifted"].sum())
            rel = abs(mc["p_click"] - an["p_click"]) / max(an["p_click"], 1e-12)
            ok_mc &= rel < 0.10
            print(f"      {water:<13}{d:>6.1f}{lv:>4}{an['p_click']:>13.4e}"
                  f"{mc['p_click']:>13.4e}{an['qber']*100:>9.2f}%"
                  f"{mc['qber']*100:>9.2f}%{n_sift:>8d}")

    print(f"\n  (a) khe hở Jensen < 5%: {'ĐẠT' if ok_lin else 'CẦN XEM LẠI'}"
          f"   → RTL được phép tuyến tính hoá")
    print(f"  (b) P_click MC khớp lý thuyết < 10%: "
          f"{'ĐẠT' if ok_mc else 'CẦN XEM LẠI'}")
    return ok_lin and ok_mc


def verify_window_design(cfg: LinkConfig, water="clear_ocean", d=15.0):
    """
    THE MOST IMPORTANT DESIGN RESULT for channel_monitor.v.

    The underwater QBER estimation window is SQUEEZED FROM BOTH SIDES:
      • too SHORT → too few clicks, QBER drowns in shot noise;
                     the controller reacts to noise rather than to the channel.
      • too LONG  → beyond the coherence time, fading is averaged away;
                     the controller no longer "sees" the turbulence.
    The optimal window ≈ one coherence block, and it is only usable if the click
    count in that block is large enough (≳ 100).  With the 256 attempts of the
    current atmospheric RTL this FAILS → channel_monitor.v must be changed.
    """
    _hr("KIỂM CHỨNG 7 — thiết kế cửa sổ ước lượng QBER (kẹp giữa 2 phía)")
    print(f"  {WATER_TYPES[water].name}, d = {d} m, f_rep = {cfg.f_rep:.0e} Hz, "
          f"τ_coh = {cfg.tau_coh_ms} ms")
    print(f"  → 1 khối kết hợp = {cfg.coherence_pulses:,} xung\n")

    print(f"  {'cửa sổ':>10}{'clicks/cs':>12}{'QBER mean':>12}"
          f"{'std(L1)':>10}{'std(L5)':>10}{'Δstd':>9}   nhận xét")
    print("  " + "-" * 76)
    MIN_CLICKS = 50          # below this, the QBER estimation error exceeds ~7% absolute
    windows = [256, 1024, 4096, 16384, cfg.coherence_pulses, 500_000]
    best_w = None
    for w in windows:
        s1 = UWOCChannel(water, d, 1, cfg).window_statistics(2500, w, seed=3)
        s5 = UWOCChannel(water, d, 5, cfg).window_statistics(2500, w, seed=3)
        gap = s5["qber_std"] - s1["qber_std"]
        usable = (s1["clicks_per_window"] >= MIN_CLICKS
                  and w <= 2 * cfg.coherence_pulses)
        if usable:
            best_w = w                      # keep the LARGEST usable window
        note = ("shot-noise chi phối" if s1["clicks_per_window"] < MIN_CLICKS
                else "fading bị trung bình hoá" if w > 2 * cfg.coherence_pulses
                else "VÙNG DÙNG ĐƯỢC")
        tag = " ← khối kết hợp" if w == cfg.coherence_pulses else ""
        print(f"  {w:>10,}{s1['clicks_per_window']:>12.1f}"
              f"{s1['qber_mean']*100:>11.2f}%{s1['qber_std']*100:>9.2f}%"
              f"{s5['qber_std']*100:>9.2f}%{gap*100:>8.2f}%   {note}{tag}")

    ok = best_w is not None
    if ok:
        print(f"\n  Cửa sổ khuyến nghị: {best_w:,} attempt "
              f"(≈ {best_w / cfg.f_rep * 1e3:.2f} ms, "
              f"{best_w / cfg.coherence_pulses:.2f}× khối kết hợp)")
        print(f"  RTL hiện tại (channel_monitor.v) dùng 256 attempt "
              f"→ nhỏ hơn {best_w // 256}× ⇒ PHẢI SỬA (bước 3).")
    else:
        print("\n  KHÔNG có cửa sổ dùng được ở điểm làm việc này:")
        print("  cần tăng f_rep, tăng μ, hoặc rút ngắn cự ly.")
    print(f"  {'ĐẠT' if ok else 'CẦN XEM LẠI'} — tồn tại vùng cửa sổ dùng được")
    return ok


def verify_window_statistics(cfg: LinkConfig, water="clear_ocean", d=15.0):
    """Turbulence shows up in the window-to-window variability, not in the mean QBER."""
    _hr("KIỂM CHỨNG 8 — nhiễu loạn lộ ra ở ĐÂU? (mean vs biến động cửa sổ)")
    w = cfg.coherence_pulses
    print(f"  {WATER_TYPES[water].name}, d = {d} m, cửa sổ = 1 khối kết hợp "
          f"({w:,} xung)\n")
    print(f"  {'Level':<11}{'σ²_ho':>9}{'clicks/cs':>11}{'QBER mean':>12}"
          f"{'QBER std':>11}{'QBER p90':>11}{'outage':>9}{'SKR/pulse':>12}")
    print("  " + "-" * 78)
    means, stds = [], []
    for lv in range(1, 6):
        ch = UWOCChannel(water, d, lv, cfg)
        st = ch.window_statistics(n_windows=6000, window=w, seed=31 + lv)
        means.append(st["qber_mean"]); stds.append(st["qber_std"])
        print(f"  {OCEAN_TURB_LEVELS[lv]['name']:<11}{ch.sigma2_ho:>9.3f}"
              f"{st['clicks_per_window']:>11.0f}{st['qber_mean']*100:>11.2f}%"
              f"{st['qber_std']*100:>10.2f}%{st['qber_p90']*100:>10.2f}%"
              f"{st['outage']*100:>8.1f}%{st['skr_per_pulse']:>12.2e}")
    spread_mean = (max(means) - min(means)) / max(np.mean(means), 1e-9)
    spread_std = (max(stds) - min(stds)) / max(np.mean(stds), 1e-9)
    ok = spread_std > 2 * spread_mean
    print(f"\n  Biến thiên tương đối theo level:  QBER mean = {spread_mean*100:.1f}%"
          f"   |   QBER std = {spread_std*100:.1f}%")
    print(f"  std nhạy hơn mean ≥2×: {'ĐẠT' if ok else 'CẦN XEM LẠI'}")
    print("\n  ⇒ KẾT LUẬN THIẾT KẾ (dùng cho bước 3):")
    print("    channel_monitor.v phải xuất thêm ĐỘ BIẾN ĐỘNG giữa các cửa sổ")
    print("    (qber_var) và LOSS RATE, chứ không chỉ QBER trung bình — nếu")
    print("    không, adaptive_controller sẽ 'mù' trước nhiễu loạn dưới nước.")
    return ok


def report_secure_range(cfg: LinkConfig):
    """Table of maximum secure range — used directly in the paper."""
    _hr("BẢNG A — cự ly bảo mật tối đa d_max [m]")
    print("  Hai tiêu chí:")
    print("    (i)  mean   : QBER trung bình dài hạn < 11%")
    print(f"    (ii) outage : P(QBER cửa sổ {cfg.coherence_pulses:,}-xung > 11%)"
          f" < 10%  ← thực tế hơn\n")
    for crit in ("mean", "outage"):
        print(f"  ── tiêu chí '{crit}' " + "─" * 46)
        print(f"  {'Nước':<17}" + "".join(f"{'L' + str(l):>10}" for l in range(1, 6)))
        for water in ("clear_ocean", "coastal", "harbor"):
            row = [max_secure_distance(water, lv, cfg, criterion=crit)
                   for lv in range(1, 6)]
            print(f"  {WATER_TYPES[water].name:<17}"
                  + "".join(f"{v:>10.2f}" for v in row))
        print()
    print("  Cột 'outage' GIẢM dần theo level — đây mới là 'secure operating")
    print("  range' đúng nghĩa, và là chỉ tiêu để so sánh Fixed vs Adaptive")
    print("  (tương ứng dòng cuối bảng Key Results trong README).")


def make_plots(cfg: LinkConfig, outdir: str):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)

    # --- Fig A: QBER & SKR vs distance, 3 water types ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    for water, col in (("clear_ocean", "#00DDFF"), ("coastal", "#FFB000"),
                       ("harbor", "#FF3333")):
        wt = WATER_TYPES[water]
        ds = np.linspace(wt.d_grid[0], wt.d_grid[-1], 60)
        qs, rs = [], []
        for d in ds:
            m = UWOCChannel(water, float(d), 3, cfg).mean_metrics()
            qs.append(m["qber"] * 100)
            rs.append(max(m["skr"], 1e-3))
        ax1.plot(ds, qs, color=col, lw=2.2, label=wt.name)
        ax2.semilogy(ds, rs, color=col, lw=2.2, label=wt.name)
    ax1.axhline(11, color="w", ls="--", lw=1.5, label="BB84 limit 11%")
    ax1.set_xlabel("Link distance d [m]"); ax1.set_ylabel("QBER [%]")
    ax1.set_ylim(0, 30); ax1.set_title("QBER vs distance (turb level 3)")
    ax2.set_xlabel("Link distance d [m]"); ax2.set_ylabel("SKR [bps]")
    ax2.set_title("Secure key rate vs distance")
    for a in (ax1, ax2):
        a.grid(alpha=.3); a.legend()
    fig.tight_layout()
    p = os.path.join(outdir, "uwoc_qber_skr_vs_distance.png")
    fig.savefig(p, dpi=160); plt.close(fig)
    print(f"  ✔ {p}")

    # --- Fig B: σ²_ho vs distance cho 5 level ---
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ds = np.linspace(2, 80, 40)
    for lv in range(1, 6):
        t = OCEAN_TURB_LEVELS[lv]
        v = [sigma2_ho(float(d), t["eps"], t["chi_T"], t["w"], cfg.lam_nm)
             for d in ds]
        ax.semilogy(ds, v, lw=2, label=f"L{lv} {t['name']}")
    ax.axhline(1.0, color="k", ls="--", lw=1, label="weak/strong boundary")
    ax.set_xlabel("d [m]"); ax.set_ylabel(r"$\sigma^2_{h_o}$")
    ax.set_title("Oceanic turbulence scintillation index (Nikishov spectrum)")
    ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(outdir, "uwoc_sigma2_ho.png")
    fig.savefig(p, dpi=160); plt.close(fig)
    print(f"  ✔ {p}")


def main():
    ap = argparse.ArgumentParser(description="UWOC channel model — self-test")
    ap.add_argument("--plot", action="store_true", help="Xuất hình kiểm chứng")
    ap.add_argument("--outdir", default="python/Images")
    ap.add_argument("--lam", type=int, default=450, choices=LAMBDA_OPTIONS)
    ap.add_argument("--mu", type=float, default=0.1)
    args = ap.parse_args()

    cfg = LinkConfig(lam_nm=args.lam, mu=args.mu)

    print("█" * 78)
    print("  UWOC CHANNEL MODEL — SELF-TEST")
    print(f"  λ={cfg.lam_nm} nm  μ={cfg.mu}  η_det={cfg.eta_det}  "
          f"Y₀={cfg.Y0:.2e}  f_rep={cfg.f_rep:.0e} Hz")
    print("█" * 78)

    r = [verify_sigma_s_fit(), verify_weibull_inversion(), verify_sigma_ho()]
    verify_qber_vs_distance(cfg)
    r.append(verify_wavelength_choice(cfg))
    r.append(verify_monte_carlo(cfg))
    r.append(verify_window_design(cfg))
    r.append(verify_window_statistics(cfg))
    report_secure_range(cfg)

    _hr("TỔNG KẾT")
    print(f"  {sum(r)}/{len(r)} kiểm chứng định lượng ĐẠT")
    print("  Bước tiếp theo:  python python/uwoc_lut_gen.py")

    if args.plot:
        _hr("XUẤT HÌNH")
        make_plots(cfg, args.outdir)


if __name__ == "__main__":
    main()
