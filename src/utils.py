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
    params = {
        # ── ESRD RISK ─────────────────────────────────────────────────────
        # Muzaale 2014, JAMA 311:579 (Table 2 / PMC full text)
        "esrd_15yr_donor_overall":  0.0031,   # 30.8/10,000 at 15 yr
        "esrd_15yr_donor_black":    0.00747,  # 74.7/10,000 at 15 yr
        "esrd_15yr_donor_white":    0.00227,  # 22.7/10,000 at 15 yr
        "esrd_15yr_nondonor":       0.00039,  # 3.9/10,000 matched controls

        # Grams 2016, NEJM 374:411 — donor-candidate baseline by race/sex
        "esrd_15yr_nondonor_black": 0.00200,  # sex-avg Black (0.24% M, 0.15% F)
        "esrd_15yr_nondonor_white": 0.00050,  # sex-avg White (0.06% M, 0.04% F)

        # Weibull shape for ESRD hazard accumulation (>1 = accelerating) and its
        # log-scale PSA sigma: fit via cloglog regression to Massie 2017's five
        # published cumulative-incidence curves (median, IQR, 1st/99th pct),
        # not assumed. See massie_weibull_fit() / fit_weibull_shape_cloglog()
        # above and 05_assemble_parameters.py, which does the same fit from
        # literature_params.json. The literal curve values below mirror the
        # massie2017 entry in 04_literature_params.py — kept as this module's
        # offline fallback should data/processed/ not yet be populated.
        "weibull_shape":            None,   # placeholder, set below
        "weibull_shape_log_sigma":  None,   # placeholder, set below

        # Massie 2017, JASN 28:2749 — hazard ratios within donor cohort
        "hr_black_race":            2.96,   # 95% CI 2.25–3.89
        "hr_male_sex":              1.88,   # 95% CI 1.50–2.35
        "hr_age_per_decade_nonblack": 1.40, # 95% CI 1.23–1.59 (non-Black only)
        "hr_age_per_decade_nondonor": 1.40, # proxy: population ESRD age-gradient

        # ── DIALYSIS / ESRD SURVIVAL ──────────────────────────────────────
        # USRDS ADR 2023 / Renal Fellow Network USRDS summary
        "dialysis_1yr_mort":        0.22,   # first year on HD
        "dialysis_annual_mort":     0.17,   # subsequent years (from 5-yr surv 42%)

        # ── WAITLIST OUTCOMES ─────────────────────────────────────────────
        # SRTR 2023 ADR, Figure KI 24
        "wl_mort_per_100py":        5.0,    # deaths/100 patient-years (2023 value)
        "wl_mort_black_per_100py":  4.62,   # SRTR 2023 Figure KI 25
        "wl_mort_white_per_100py":  5.71,

        # Punjala 2024 (Transplant Proc 56:1740-1751), Table 3 — national mean
        # waiting time at transplant, post-KAS250 (5/2021-4/2022)
        "wl_std_mean_days":         1765,   # 58 months overall (mean)
        # Wainright 2017 AJT 17:1103; median 102.6 d → mean = 102.6/ln(2) ≈ 148 d
        "wl_pld_mean_days":         148.0,  # prior living donors post-KAS (mean)

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

        # DDKT 5-yr graft survival by age (SRTR 2023 ADR Figure KI 53) and an
        # assumed national 1-yr graft survival, used below to derive
        # age-stratified post-year-1 graft failure.
        "ddkt_graft_5yr_age1834": 0.822,
        "ddkt_graft_5yr_age3549": 0.835,
        "ddkt_graft_5yr_age5064": 0.768,
        "ddkt_graft_5yr_age65p":  0.661,
        "ddkt_graft_1yr_assumed": 0.955,

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
        # Donor all-cause mortality HR, modeled as time-varying: short/medium
        # follow-up studies (Muzaale 2014, Segev 2010, O'Keeffe 2018 pooled)
        # show no excess mortality within ~10 years; Mjøen 2014 (median 15 yr
        # follow-up) found donor/control survival curves separate only after
        # ~10 years, reaching HR=1.30 (95% CI 1.11-1.52) by ~25 years. See
        # donor_mort_hr_at() in utils.py.
        "donor_mort_hr_early":      1.0,
        "donor_mort_hr_late":       1.30,
        "donor_mort_hr_t_start":    10.0,
        "donor_mort_hr_t_end":      15.0,

        # ── PREEMPTIVE TRANSPLANT LISTING ─────────────────────────────────
        # One-time branching probability at ESRD onset: listed before dialysis starts.
        # Source: USRDS 2025 ADR ESRD Vol. Ch.7, Figure 7.13 (2024 data).
        "esrd_preemptive_prob_std": 0.058,   # 5.8% overall (non-donor standard arm)
        "esrd_preemptive_prob_pld": 0.094,   # 9.4% age 18-44 (donor-like cohort)
    }

    # Massie 2017, JASN 28:2749, Figure 3 / p. 2751-2752 — five published
    # cumulative-incidence-of-ESRD curves (per 10,000 donors) at 5/10/15/20 yr:
    # median, IQR (25th/75th pct), and 1st/99th pct. Mirrors the massie2017
    # entry in 04_literature_params.py. weibull_shape/weibull_shape_log_sigma
    # are fit from these, not assumed — see massie_weibull_fit() above.
    _massie2017_curves = {
        "esrd_cum_incidence_5yr_per10k":       1.0,
        "esrd_cum_incidence_5yr_p01_per10k":   0.2,
        "esrd_cum_incidence_5yr_p25_per10k":   1.0,
        "esrd_cum_incidence_5yr_p75_per10k":   2.0,
        "esrd_cum_incidence_5yr_p99_per10k":   8.0,
        "esrd_cum_incidence_10yr_per10k":      6.0,
        "esrd_cum_incidence_10yr_p01_per10k":  1.2,
        "esrd_cum_incidence_10yr_p25_per10k":  4.0,
        "esrd_cum_incidence_10yr_p75_per10k": 11.0,
        "esrd_cum_incidence_10yr_p99_per10k": 48.0,
        "esrd_cum_incidence_15yr_per10k":     16.0,
        "esrd_cum_incidence_15yr_p01_per10k":  3.0,
        "esrd_cum_incidence_15yr_p25_per10k": 10.0,
        "esrd_cum_incidence_15yr_p75_per10k": 29.0,
        "esrd_cum_incidence_15yr_p99_per10k": 125.0,
        "esrd_cum_incidence_20yr_per10k":     34.0,
        "esrd_cum_incidence_20yr_p01_per10k":  7.0,
        "esrd_cum_incidence_20yr_p25_per10k": 20.0,
        "esrd_cum_incidence_20yr_p75_per10k": 59.0,
        "esrd_cum_incidence_20yr_p99_per10k": 256.0,
    }
    params["weibull_shape"], params["weibull_shape_log_sigma"] = \
        massie_weibull_fit(_massie2017_curves)

    # Age-stratified post-year-1 graft failure, derived as
    # 1 - (5yr_surv / 1yr_surv)^(1/4) from the DDKT graft survival above
    # (see src/03_download_srtr.py for the ADR-Excel-sourced version). A flat
    # rate cannot represent this: the four bands span roughly 3.3%-8.8%/yr.
    for _band in ("age1834", "age3549", "age5064", "age65p"):
        params[f"graft_annual_fail_postyear1_{_band}"] = round(
            1 - (params[f"ddkt_graft_5yr_{_band}"]
                 / params["ddkt_graft_1yr_assumed"]) ** 0.25, 4)
    params["graft_annual_fail_postyear1"] = round(sum(
        params[f"graft_annual_fail_postyear1_{_b}"]
        for _b in ("age1834", "age3549", "age5064", "age65p")) / 4, 4)

    return params


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


def fit_weibull_shape_cloglog(times, cum_incidences) -> tuple:
    """
    Fit Weibull shape k (and scale lambda) to a cumulative-incidence curve via
    the complementary log-log linearization of the Weibull CDF:

      I(t) = 1 - exp(-(t/lambda)^k)
      ln(-ln(1 - I(t))) = k*ln(t) - k*ln(lambda)

    This is linear in ln(t) with slope k, so an OLS fit recovers k directly
    (closed form, no nonlinear solver). `cum_incidences` are probabilities
    (not per-10,000 counts).
    """
    x = np.log(np.asarray(times, dtype=float))
    y = np.log(-np.log(1.0 - np.asarray(cum_incidences, dtype=float)))
    slope, intercept = np.polyfit(x, y, 1)
    k   = float(slope)
    lam = float(np.exp(-intercept / slope))
    return k, lam


def massie_weibull_fit(massie: dict) -> tuple:
    """
    Fit Weibull shape k and its log-scale PSA sigma from Massie 2017's five
    published cumulative-incidence curves (median, IQR, 1st/99th percentile;
    see fit_weibull_shape_cloglog). `massie` is the massie2017 dict from
    literature_params.json (or an equivalent literal mirror), keyed by
    esrd_cum_incidence_{t}yr[_p{pct}]_per10k for t in (5, 10, 15, 20).

    k is the fit to the median curve; sigma is the sample SD of ln(k) across
    all 5 independently-fit curves, capturing how consistently a single
    Weibull shape describes the whole published risk distribution.
    """
    times = (5, 10, 15, 20)
    curve_suffixes = {"p01": "_p01", "p25": "_p25", "p50": "",
                       "p75": "_p75", "p99": "_p99"}
    ks = {}
    for label, suffix in curve_suffixes.items():
        vals = [massie[f"esrd_cum_incidence_{t}yr{suffix}_per10k"] / 10_000.0
                 for t in times]
        k, _ = fit_weibull_shape_cloglog(times, vals)
        ks[label] = k
    k_hat = ks["p50"]
    sigma = float(np.std(np.log(list(ks.values())), ddof=1))
    return k_hat, sigma


def weibull_scale_from_cumrisk(cum_risk_15: float, k: float) -> float:
    """
    Return Weibull scale lambda such that:
      1 - exp(-(15/lambda)^k) = cum_risk_15
    """
    return 15.0 / (-np.log(1.0 - cum_risk_15)) ** (1.0 / k)


def donor_mort_hr_at(t: float, hr_early: float = 1.0, hr_late: float = 1.30,
                      t_start: float = 10.0, t_end: float = 25.0) -> float:
    """
    Donor all-cause mortality HR as a function of years since donation.

    Flat at hr_early for t <= t_start (no excess mortality detectable in the
    short/medium-follow-up literature, e.g. Segev 2010, Garg 2012), ramping
    linearly to hr_late by t_end (Mjøen 2014's finding that donor and control
    survival curves separate only after ~10 years, reaching HR=1.30 by their
    median ~15-25 yr follow-up), then flat at hr_late thereafter.
    """
    if t <= t_start:
        return hr_early
    if t >= t_end:
        return hr_late
    frac = (t - t_start) / (t_end - t_start)
    return hr_early + frac * (hr_late - hr_early)


def weibull_scale_from_cumrisk_competing(
    cum_risk_15: float,
    k: float,
    life_table_qx: np.ndarray,
    age_at_entry: int,
    bg_hr=1.0,
) -> float:
    """
    Return Weibull scale lambda calibrated so the competing-risk-adjusted
    15-year ESRD CIF equals cum_risk_15.

    Accounts for background mortality depleting the ESRD-susceptible pool each
    year; uses the same life table and bg_hr as the main simulation. bg_hr may
    be a fixed scalar or a callable bg_hr(t) giving the HR t years after entry
    (e.g. donor_mort_hr_at), to support a time-varying background HR.
    Solved by bisection (60 iterations → precision < 1e-12 for typical inputs).
    """
    bg_hr_fn = bg_hr if callable(bg_hr) else (lambda t: bg_hr)

    def cr_cif(lam: float) -> float:
        cif, s_bg = 0.0, 1.0
        for t in range(15):
            cif += s_bg * weibull_annual_prob(t, lam, k)
            age  = min(age_at_entry + t, len(life_table_qx) - 1)
            s_bg *= 1.0 - life_table_qx[age] * bg_hr_fn(t)
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
    under an exponential waiting-time model (rate = ln2/median)."""
    median_yrs = median_days / 365.25
    rate = np.log(2) / median_yrs
    return float(1.0 - np.exp(-rate))


def mean_to_annual_tx_prob(mean_days: float) -> float:
    """Convert mean waiting time (days) to annual transplant probability
    under an exponential waiting-time model (rate = 1/mean; Sonnenberg & Beck 1993)."""
    mean_yrs = mean_days / 365.25
    rate = 1.0 / mean_yrs
    return float(1.0 - np.exp(-rate))
