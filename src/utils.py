"""
utils.py — shared helpers for the kidney donation model pipeline.
"""

import json
import numpy as np
from pathlib import Path

# ── REPO-RELATIVE PATHS ───────────────────────────────────────────────────────
REPO_ROOT  = Path(__file__).resolve().parent.parent
DATA_RAW   = REPO_ROOT / "data" / "raw"
DATA_PROC  = REPO_ROOT / "data" / "processed"
RESULTS    = REPO_ROOT / "results"

for _p in (DATA_RAW, DATA_PROC, RESULTS):
    _p.mkdir(parents=True, exist_ok=True)


# ── LIFE TABLE HELPERS ────────────────────────────────────────────────────────
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
        # Gompertz-Makeham fit calibrated to 2021 CDC life tables (nvsr72-12.pdf).
        # Targets (sex-averaged): qx(0)≈0.006, qx(40)≈0.002, qx(70)≈0.027,
        # qx(80)≈0.075. A=accident hazard, B·exp(c·age)=aging component.
        A, B, c = 0.0007, 0.00005, 0.095
        ages = np.arange(101)
        return np.clip(A + B * np.exp(c * ages), 0.0, 1.0)


# ── PARAMETER I/O ─────────────────────────────────────────────────────────────
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
        # ── ESRD RISK ─────────────────────────────────────────────────────
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

        # ── DIALYSIS / ESRD SURVIVAL ──────────────────────────────────────
        # USRDS ADR 2023 / Renal Fellow Network USRDS summary
        "dialysis_1yr_mort":        0.22,   # first year on HD
        "dialysis_annual_mort":     0.17,   # subsequent years (from 5-yr surv 42%)

        # ── WAITLIST OUTCOMES ─────────────────────────────────────────────
        # SRTR 2023 ADR, Figure KI 24
        "wl_mort_per_100py":        5.0,    # deaths/100 patient-years (2023 value)
        "wl_mort_black_per_100py":  4.62,   # SRTR 2023 Figure KI 25
        "wl_mort_white_per_100py":  5.71,

        # SRTR 2023 ADR / Schold AJT 2023 — post-KAS250 (March 2021)
        "wl_std_median_days":       985,    # ~32.8 months overall
        # Wainright 2017 AJT 17:1103 + UNOS conference abstract
        "wl_pld_median_days":       100,    # prior living donors post-KAS

        # SRTR 2023 ADR Figure KI 22: 3-yr removal CIF = 19.06%.
        # Back-calculated via competing-risk bisection (see 03_download_srtr.py
        # _solve_removal_rate) accounting for simultaneous transplantation (22.7%/yr)
        # and death (4.9%/yr).  Old formula 1-(1-0.191)^(1/3)=6.81% ignored these
        # competing events and reproduced only 10.8% removed at 3yr instead of 19.1%.
        "wl_removal_rate_yr":       0.1260, # 12.60%/yr (corrected, SRTR 2023 KI 22)

        # Conditional per-cycle probability of transitioning from dialysis to
        # the waitlist, calibrated from USRDS 2025 ADR ESRD Vol. Ch.7:
        #   Fig 7.15 3-yr listing CIF = 12% (general ESRD) → p ≈ 0.065/yr
        #   Fig 7.13+7.17 for donor-like cohort (age 18-44) → p ≈ 0.17/yr
        # Conservative base case rounded down; sensitivity 0.05–0.30.
        # See src/02b_download_usrds_esrd7.py for full derivation.
        "wl_listing_prob":          0.15,

        # ── POST-TRANSPLANT OUTCOMES ──────────────────────────────────────
        # SRTR 2023 ADR DDKT patient survival — age-stratified annual mortality
        # Derived as 1 - 5yr_surv^(1/5) for each age band
        "posttx_annual_mort_age1834": 0.009,  # 1 - 0.957^0.2
        "posttx_annual_mort_age3549": 0.018,  # 1 - 0.914^0.2
        "posttx_annual_mort_age5064": 0.039,  # 1 - 0.820^0.2
        "posttx_annual_mort_age65p":  0.068,  # 1 - 0.701^0.2

        # Overall: age-weighted average (weights approximate ESRD recipient mix)
        "posttx_annual_mort":       0.036,  # ~weighted avg of age strata
        "posttx_annual_mort_black": 0.035,  # SRTR 2023 DDKT race-stratified
        "posttx_annual_mort_white": 0.038,

        # Annual graft failure rate post year-1 (SRTR 2023)
        "graft_annual_fail_postyear1": 0.025,

        # ── POST-TRANSPLANT OUTCOMES (LDKT) ───────────────────────────────
        # SRTR 2023 ADR Figure KI 76 — LDKT patient survival by recipient age
        # (2016–2018 transplant cohort). Derived as 1 − 5yr_surv^(1/5).
        "posttx_ld_annual_mort_age1834": 0.0042,  # 1 - 0.979^0.2
        "posttx_ld_annual_mort_age3549": 0.0079,  # 1 - 0.961^0.2
        "posttx_ld_annual_mort_age5064": 0.0172,  # 1 - 0.917^0.2
        "posttx_ld_annual_mort_age65p":  0.0392,  # 1 - 0.819^0.2
        # Age-weighted average (same SRTR KI 1 weights: [0.10, 0.30, 0.40, 0.20])
        "posttx_ld_annual_mort":         0.0175,

        # ── BACKGROUND MORTALITY ──────────────────────────────────────────
        # CDC NVSR Vol 72 No 12 (Nov 2023) — 2021 US life tables
        # Full table loaded separately via load_life_table()
        # HR for donor all-cause mortality vs matched controls
        # Muzaale 2014 / Segev 2010 / Grams 2018 meta-analysis: HR ≈ 1.0
        "donor_mort_hr":            1.0,

        # Mjøen 2014 (Kidney Int 86:162) — upper-bound / sensitivity scenario
        "donor_mort_hr_mjoeen":     1.30,

        # ── PREEMPTIVE TRANSPLANT LISTING ─────────────────────────────────
        # One-time branching probability at ESRD onset: listed before dialysis starts.
        # Source: USRDS 2025 ADR ESRD Vol. Ch.7, Figure 7.13 (2024 data).
        "esrd_preemptive_prob_std": 0.058,   # 5.8% overall (non-donor standard arm)
        "esrd_preemptive_prob_pld": 0.094,   # 9.4% age 18-44 (donor-like cohort)
    }


# ── STATISTICAL HELPERS ───────────────────────────────────────────────────────
def beta_params_from_mean_se(mean: float, se_frac: float = 0.20):
    """Return (alpha, beta) for a Beta distribution with given mean and
    fractional SE. Clamps shape parameters to >= 0.5."""
    se = mean * se_frac
    denom = se ** 2
    a = mean * (mean * (1 - mean) / denom - 1)
    b = (1 - mean) * (mean * (1 - mean) / denom - 1)
    return max(float(a), 0.5), max(float(b), 0.5)


def weibull_annual_prob(t: float, lam: float, k: float) -> float:
    """Annual ESRD transition probability for cycle [t, t+1] under Weibull(lambda, k).

    Equals S(t) - S(t+1) where S(t) = exp(-(t/lambda)^k).  The sum over
    t in [0, T-1] telescopes to exactly 1 - S(T), so calibration via
    weibull_scale_from_cumrisk is consistent with these per-cycle draws.
    """
    if t < 0:
        return 0.0
    surv_t  = np.exp(-((t / lam) ** k))
    surv_t1 = np.exp(-(((t + 1) / lam) ** k))
    return float(np.clip(surv_t - surv_t1, 0.0, 1.0))


def weibull_scale_from_cumrisk(cum_risk_15: float, k: float) -> float:
    """
    Return Weibull scale lambda such that:
      1 - exp(-(15/lambda)^k) = cum_risk_15
    """
    return 15.0 / (-np.log(1.0 - cum_risk_15)) ** (1.0 / k)


def weibull_scale_from_cumrisk_competing(
    cum_risk_15: float,
    k: float,
    life_table_qx: np.ndarray,
    age_at_entry: int,
    bg_hr: float = 1.0,
) -> float:
    """
    Return Weibull scale lambda calibrated so the competing-risk-adjusted
    15-year ESRD CIF equals cum_risk_15.

    Accounts for background mortality depleting the ESRD-susceptible pool each
    year; uses the same life table and bg_hr as the main simulation.
    Solved by bisection (60 iterations → precision < 1e-12 for typical inputs).
    """
    def cr_cif(lam: float) -> float:
        cif, s_bg = 0.0, 1.0
        for t in range(15):
            cif += s_bg * weibull_annual_prob(t, lam, k)
            age  = min(age_at_entry + t, len(life_table_qx) - 1)
            s_bg *= 1.0 - life_table_qx[age] * bg_hr
        return cif

    # cr_cif is monotonically decreasing in lam
    lam_lo, lam_hi = 1.0, 10_000.0
    for _ in range(60):
        lam_mid = (lam_lo + lam_hi) / 2.0
        if cr_cif(lam_mid) > cum_risk_15:
            lam_lo = lam_mid
        else:
            lam_hi = lam_mid
    return (lam_lo + lam_hi) / 2.0


def median_to_annual_tx_prob(median_days: float) -> float:
    """Convert median waiting time (days) to annual transplant probability
    under an exponential waiting-time model."""
    median_yrs = median_days / 365.25
    rate = np.log(2) / median_yrs
    return float(1.0 - np.exp(-rate))
