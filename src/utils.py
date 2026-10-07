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
        # Grams 2016 projections by age (20/40/60) — non-donor age gradient
        **{f"grams_nondonor_15yr_{cell}_age{age}": pct / 100
           for cell, row in {"black_male":   (0.08, 0.24, 0.32),
                             "black_female": (0.05, 0.15, 0.18),
                             "white_male":   (0.02, 0.06, 0.13),
                             "white_female": (0.01, 0.04, 0.08)}.items()
           for age, pct in zip((20, 40, 60), row)},

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

        # ── DIALYSIS / ESRD SURVIVAL ──────────────────────────────────────
        # USRDS 2023 ADR ESRD Ch.6 Fig 6.5b — dialysis mortality per 1,000 PY,
        # 2019 (pre-COVID), adjusted; annual prob = 1 − exp(−rate/1000).
        # See src/02c_usrds_mortality.py.
        "dialysis_mort_age1844":    0.06819,  # 70.6/1000 PY
        "dialysis_mort_age4564":    0.11438,  # 121.5/1000 PY
        "dialysis_mort_age6574":    0.17765,  # 195.6/1000 PY
        "dialysis_mort_age75p":     0.25142,  # 289.6/1000 PY
        # Fig 6.4 — 2019 incident HD cohort, year-1 ÷ year-2 rate (214.1/201.9)
        "dialysis_yr1_mult":        1.0606,
        "dialysis_mort_scale":      1.0,      # uniform sensitivity scale

        # ── WAITLIST OUTCOMES ─────────────────────────────────────────────
        # SRTR 2023 ADR, Figure KI 24
        "wl_mort_per_100py":        5.0,    # deaths/100 patient-years (2023 value)
        "wl_mort_black_per_100py":  4.62,   # SRTR 2023 Figure KI 26
        "wl_mort_white_per_100py":  5.71,
        # By age (Figure KI 25, 2019 pre-COVID) — used by the model
        "wl_mort_per_100py_age1834": 1.98,
        "wl_mort_per_100py_age3549": 2.95,
        "wl_mort_per_100py_age5064": 4.97,
        "wl_mort_per_100py_age65p":  7.43,
        # Race ÷ overall, 2019 (KI 26 / KI 24); applied only in race subgroups
        "wl_mort_race_mult":        1.0,
        "wl_mort_race_mult_black":  4.47 / 4.78,
        "wl_mort_race_mult_white":  5.57 / 4.78,

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

        # ── PRIOR-DONOR POST-ESRD PATHWAY ─────────────────────────────────
        # Muzaale 2016, Transplantation 100:1306 — 99 donors who developed ESRD.
        # Donor arm only; voucher / ESRD-conditional priority arms do not use these.
        "donor_esrd_preemptive_prob": 21 / 99,  # 20 preemptive listings + 1 preemptive LDKT
        "donor_posttx_mort_mult":     0.7,      # adjusted post-Tx mortality HR (0.2–2.4)
        # Calibrated in 05_assemble_parameters.py (calibrate_donor_esrd_pathway):
        # listing to 43.4% listed by 12 mo; multiplier to 1/3/5/10-yr mortality
        "donor_wl_listing_prob":      0.3081,
        "donor_dialysis_mort_mult":   0.78,
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


# ── POST-ESRD STATE MORTALITY ─────────────────────────────────────────────────
# Shared by every model implementation (06 Monte Carlo + analytic, 07, 09, 10)
# so the age dependence, donor adjustments, and background floor stay in sync.
#
# Registry inputs come in age bands. Rates are placed at each band's midpoint
# and interpolated log-linearly between midpoints (held flat beyond the outer
# ones). Step functions created artifacts at band edges: e.g. post-Tx
# mortality jumping from 3.9% to 6.9% at exactly 65 made the ESRD-conditional
# priority benefit rise between ages 65 and 73. Open-ended bands use an
# assumed representative age: 31 for 18-44, 26 for 18-34, 80 for 75+, 70 for
# 65+ (most recipients and candidates aged 65+ are under 75).
#
# Every post-ESRD mortality is floored at the general-population life-table qx
# for the same age (q_floor). Without the floor, fixed state-specific rates
# fall below background mortality at old ages (e.g. 65+ post-Tx 6.85%/yr is
# exceeded by qx from age 83), making ESRD paradoxically protective.

USRDS_BAND_MIDS = ((31.0, "age1844"), (54.5, "age4564"), (69.5, "age6574"), (80.0, "age75p"))
SRTR_BAND_MIDS  = ((26.0, "age1834"), (42.0, "age3549"), (57.0, "age5064"), (70.0, "age65p"))


def interp_age_bands(p: dict, prefix: str, mids, age: float) -> float:
    """Log-linear interpolation of p[f"{prefix}{band}"] between band midpoints."""
    xs = [m for m, _ in mids]
    ys = [np.log(float(p[f"{prefix}{band}"])) for _, band in mids]
    return float(np.exp(np.interp(age, xs, ys)))


def dialysis_annual_mort(p: dict, age: int, first_year: bool = False,
                         donor: bool = False, q_floor: float = 0.0) -> float:
    """
    Annual all-cause mortality on dialysis at a given age.

    Age bands 18-44 / 45-64 / 65-74 / 75+ from USRDS 2023 ADR Fig 6.5b (2019,
    adjusted), interpolated between band midpoints; first-year excess from
    Fig 6.4. donor=True applies the prior-donor multiplier calibrated to
    Muzaale 2016 (05_assemble_parameters.py). dialysis_mort_scale is a uniform
    sensitivity scale (base 1.0).
    """
    q = interp_age_bands(p, "dialysis_mort_", USRDS_BAND_MIDS, age)
    q *= p.get("dialysis_mort_scale", 1.0)
    if first_year:
        q *= p.get("dialysis_yr1_mult", 1.0)
    if donor:
        q *= p.get("donor_dialysis_mort_mult", 1.0)
    return float(min(max(q, q_floor), 1.0))


def posttx_annual_mort(p: dict, age: int, ldkt: bool = False,
                       donor: bool = False, q_floor: float = 0.0) -> float:
    """
    Age-stratified annual post-transplant mortality (SRTR 2023 ADR; DDKT base,
    LDKT for sensitivity), interpolated between band midpoints. donor=True
    applies the Muzaale 2016 adjusted post-transplant mortality HR for prior
    donors vs matched non-donors.
    """
    prefix = "posttx_ld_annual_mort_" if ldkt else "posttx_annual_mort_"
    q = interp_age_bands(p, prefix, SRTR_BAND_MIDS, age)
    if donor:
        q *= p.get("donor_posttx_mort_mult", 1.0)
    return float(min(max(q, q_floor), 1.0))


def waitlist_annual_mort_at(p: dict, age: int, q_floor: float = 0.0) -> float:
    """
    Annual waitlist mortality at a given age, floored.

    SRTR 2023 ADR Figure KI 25, 2019 (pre-COVID, matching the USRDS dialysis
    year): deaths per 100 patient-years of waiting by age band, interpolated
    between band midpoints. wl_mort_race_mult (base 1.0) applies the KI 26
    race ratio in race subgroups, assuming race and age act multiplicatively.
    """
    rate = interp_age_bands(p, "wl_mort_per_100py_", SRTR_BAND_MIDS, age) / 100.0
    rate *= p.get("wl_mort_race_mult", 1.0)
    return float(max(1.0 - np.exp(-rate), q_floor))


def grams_nondonor_age_scale(p: dict, age: float, race: str = "Overall",
                             sex: str = "Overall", reference_age: float = 40) -> float:
    """
    Ratio of Grams 2016 15-yr no-donation ESRD risk at `age` to that at
    `reference_age`, log-linearly interpolated between Grams' age 20/40/60
    projections (held flat outside). Multiplying an age-40 non-donor rate by
    this reproduces Grams exactly for Grams cells and keeps the overall arm's
    Muzaale matched-control anchor at 40.

    Race "Overall" uses the White gradient (White donors dominate the donor
    pool); sex "Overall" averages female and male. Black risk flattens after
    40 in Grams (≈1.13×/decade) while White keeps rising (≈1.45×/decade).
    """
    r = "black" if race == "Black" else "white"
    sexes = ("female", "male") if sex == "Overall" else (sex.lower(),)
    vals = [np.mean([p[f"grams_nondonor_15yr_{r}_{s}_age{a}"] for s in sexes])
            for a in (20, 40, 60)]
    curve = lambda x: np.exp(np.interp(x, [20, 40, 60], np.log(vals)))
    return float(curve(age) / curve(reference_age))


def esrd_onset_mortality_curve(p: dict, age_at_esrd: int, donor: bool,
                               years, life_table_qx: np.ndarray) -> dict:
    """
    Cumulative all-cause mortality t years after ESRD onset, for a cohort
    entering ESRD at age_at_esrd. Same D1/D2/WL/PT transitions as the main
    model (see run_arm_analytic in 06_markov_simulation.py).

    donor=True  → prior-donor pathway: priority wait, donor-specific preemptive
                  listing, listing rate, dialysis and post-Tx multipliers.
    donor=False → general ESRD pathway: standard wait and inputs.

    Used to calibrate the donor pathway to, and validate the general pathway
    against, Muzaale 2016 (Transplantation 100:1306) mortality after ESRD.
    """
    max_age = len(life_table_qx) - 1
    mean_days = p["wl_pld_mean_days"] if donor else p["wl_std_mean_days"]
    wl_tx      = mean_to_annual_tx_prob(float(mean_days))
    wl_remove  = float(p["wl_removal_rate_yr"])
    wl_listing = float(p["donor_wl_listing_prob"] if donor else p["wl_listing_prob"])
    pre = float(p["donor_esrd_preemptive_prob"] if donor else p["esrd_preemptive_prob_std"])

    D1, D2, WL, PT = 1.0 - pre, 0.0, pre, 0.0
    horizon = max(years)
    curve = {}
    for yr in range(horizon):
        age  = age_at_esrd + yr
        qx   = float(life_table_qx[min(age, max_age)])
        dm1  = dialysis_annual_mort(p, age, first_year=True,  donor=donor, q_floor=qx)
        dm   = dialysis_annual_mort(p, age, first_year=False, donor=donor, q_floor=qx)
        wlm  = waitlist_annual_mort_at(p, age, q_floor=qx)
        ptm  = posttx_annual_mort(p, age, donor=donor, q_floor=qx)
        gf   = graft_annual_fail(p, age)

        D1s = D1 * (1 - dm1); D1l = D1s * wl_listing
        D2s = D2 * (1 - dm);  D2l = D2s * wl_listing
        WLs = WL * (1 - wlm); WLt = WLs * wl_tx; WLr = (WLs - WLt) * wl_remove
        PTs = PT * (1 - ptm); PTf = PTs * gf

        D1 = 0.0
        D2 = (D1s - D1l) + (D2s - D2l) + WLr + PTf
        WL = D1l + D2l + (WLs - WLt - WLr)
        PT = WLt + (PTs - PTf)
        if yr + 1 in years:
            curve[yr + 1] = 1.0 - (D1 + D2 + WL + PT)
    return curve


def graft_annual_fail(p: dict, age: int) -> float:
    """Age-stratified post-year-1 graft failure (SRTR 2023 ADR Figure KI 53),
    interpolated between band midpoints."""
    return interp_age_bands(p, "graft_annual_fail_postyear1_", SRTR_BAND_MIDS, age)
