"""
utils.py — shared helpers for the kidney donation model pipeline.
"""

import json
import numpy as np
from pathlib import Path

# ── Repo-relative paths ───────────────────────────────────────────────────────
REPO_ROOT  = Path(__file__).resolve().parent.parent
DATA_RAW   = REPO_ROOT / "data" / "raw"
DATA_PROC  = REPO_ROOT / "data" / "processed"
RESULTS    = REPO_ROOT / "results"

for _p in (DATA_RAW, DATA_PROC, RESULTS):
    _p.mkdir(parents=True, exist_ok=True)


# ── Life table helpers ────────────────────────────────────────────────────────
def load_life_table(path: Path = None) -> np.ndarray:
    """
    Load age-specific annual probability of death (qx) from a CSV produced
    by 01_download_lifetables.py.

    Returns array of shape (101,) for ages 0–100.
    If path is None, falls back to the Gompertz-Makeham approximation used
    during development.
    """
    if path is not None and path.exists():
        import pandas as pd
        df = pd.read_csv(path)
        # Expects columns: age, qx
        df = df[df["age"] <= 100].sort_values("age")
        return df["qx"].values[:101]
    else:
        # Gompertz-Makeham fit calibrated to 2021 CDC life tables
        # A=accident hazard, B=baseline aging, c=aging rate
        A, B, c = 0.0007, 0.00005, 0.095
        ages = np.arange(101)
        return np.clip(A + B * np.exp(c * ages), 0.0, 1.0)


# ── Parameter I/O ─────────────────────────────────────────────────────────────
def load_params(path: Path = None) -> dict:
    """Load assembled parameters from params.json, or return hard-coded
    base-case values if the file does not yet exist."""
    if path is None:
        path = DATA_PROC / "params.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return _hardcoded_base_params()


def save_params(params: dict, path: Path = None) -> None:
    if path is None:
        path = DATA_PROC / "params.json"
    with open(path, "w") as f:
        json.dump(params, f, indent=2)
    print(f"Parameters saved → {path}")


def _hardcoded_base_params() -> dict:
    """
    Hard-coded base-case parameters, confirmed against primary sources.
    Used as fallback when data pipeline has not been run.

    Sources documented inline — see README and src/04_literature_params.py.
    """
    return {
        # ── ESRD risk ─────────────────────────────────────────────────────
        # Muzaale 2014, JAMA 311:579 (Table 2 / PMC full text)
        "esrd_15yr_donor_overall":  0.0031,   # 30.8/10,000 at 15 yr
        "esrd_15yr_donor_black":    0.00747,  # 74.7/10,000 at 15 yr
        "esrd_15yr_donor_white":    0.00227,  # 22.7/10,000 at 15 yr
        "esrd_15yr_nondonor":       0.00039,  # 3.9/10,000 matched controls

        # Grams 2016, NEJM 374:411 — donor-candidate baseline by race/sex
        "esrd_15yr_nondonor_black": 0.00200,  # sex-avg Black (0.24% M, 0.15% F)
        "esrd_15yr_nondonor_white": 0.00050,  # sex-avg White (0.06% M, 0.04% F)

        # Weibull shape for ESRD hazard accumulation (>1 = accelerating)
        # Calibrated so integration over 15 yr reproduces Muzaale CIF
        "weibull_shape":            1.5,

        # Massie 2017, JASN 28:2749 — hazard ratios within donor cohort
        "hr_black_race":            2.96,   # 95% CI 2.25–3.89
        "hr_male_sex":              1.88,   # 95% CI 1.50–2.35
        "hr_age_per_decade_nonblack": 1.40, # 95% CI 1.23–1.59 (non-Black only)

        # ── Dialysis / ESRD survival ──────────────────────────────────────
        # USRDS ADR 2023 / Renal Fellow Network USRDS summary
        "dialysis_1yr_mort":        0.22,   # first year on HD
        "dialysis_annual_mort":     0.17,   # subsequent years (from 5-yr surv 42%)

        # ── Waitlist outcomes ─────────────────────────────────────────────
        # SRTR 2022 ADR, Figure KI 24
        "wl_mort_per_100py":        5.4,    # deaths/100 patient-years, standard
        "wl_mort_black_per_100py":  6.5,    # approximate from SRTR Figure KI 25
        "wl_mort_white_per_100py":  5.0,

        # SRTR 2022 ADR / Schold AJT 2023 — post-KAS250 (March 2021)
        "wl_std_median_days":       985,    # ~32.8 months overall
        # Wainright 2017 AJT 17:1103 + UNOS conference abstract
        "wl_pld_median_days":       100,    # prior living donors post-KAS

        # SRTR 2022 ADR Figure KI 22: 3-yr cohort listed 2017–2019
        "wl_removal_rate_yr":       0.064,  # ~6.4%/yr competing removal

        # ── Post-transplant outcomes ──────────────────────────────────────
        # SRTR 2022 ADR, Table KI 11/12 (DDKT — used for ESRD arm)
        "graft_5yr_age1834":        0.814,
        "graft_5yr_age3564":        0.760,
        "graft_5yr_age65p":         0.678,

        # SRTR 2022 ADR approximate patient survival
        "posttx_annual_mort":       0.052,  # ~1-(0.85^(1/3)) from 3-yr surv 85%
        "posttx_annual_mort_black": 0.065,  # approximate from SRTR race data
        "posttx_annual_mort_white": 0.048,

        # Annual graft failure rate post year-1 (approximate)
        "graft_annual_fail_postyear1": 0.025,

        # ── Background mortality ──────────────────────────────────────────
        # CDC NVSR Vol 72 No 12 (Nov 2023) — 2021 US life tables
        # Full table loaded separately via load_life_table()
        # HR for donor all-cause mortality vs matched controls
        # Muzaale 2014 / Segev 2010 / Grams 2018 meta-analysis: HR ≈ 1.0
        "donor_mort_hr":            1.0,

        # Mjøen 2014 (Kidney Int 86:162) — upper-bound / sensitivity scenario
        "donor_mort_hr_mjoeen":     1.30,
    }


# ── Statistical helpers ───────────────────────────────────────────────────────
def beta_params_from_mean_se(mean: float, se_frac: float = 0.20):
    """Return (alpha, beta) for a Beta distribution with given mean and
    fractional SE. Clamps shape parameters to >= 0.5."""
    se = mean * se_frac
    denom = se ** 2
    a = mean * (mean * (1 - mean) / denom - 1)
    b = (1 - mean) * (mean * (1 - mean) / denom - 1)
    return max(float(a), 0.5), max(float(b), 0.5)


def weibull_annual_hazard(t: float, lam: float, k: float) -> float:
    """Annual hazard at time t for Weibull(lambda, k). Returns 0 at t=0."""
    if t <= 0:
        return 0.0
    return float(np.clip((k / lam) * (t / lam) ** (k - 1), 0, 0.5))


def weibull_scale_from_cumrisk(cum_risk_15: float, k: float) -> float:
    """
    Return Weibull scale lambda such that:
      1 - exp(-(15/lambda)^k) = cum_risk_15
    """
    return 15.0 / (-np.log(1.0 - cum_risk_15)) ** (1.0 / k)


def median_to_annual_tx_prob(median_days: float) -> float:
    """Convert median waiting time (days) to annual transplant probability
    under an exponential waiting-time model."""
    median_yrs = median_days / 365.25
    rate = np.log(2) / median_yrs
    return float(1.0 - np.exp(-rate))
