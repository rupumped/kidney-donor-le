"""
05_assemble_parameters.py
──────────────────────────
Combines outputs from scripts 01–04 into a single params.json used by
the simulation. Applies any cross-source reconciliation logic documented
in the parameter confirmation notes.

Key reconciliation decisions:
  1. White non-donor ESRD baseline: use Grams 2016 (0.05%) not Muzaale (0%)
  2. Wait times: national mean waiting time at transplant, Punjala 2024 (Table 3)
     — post-KAS250 (5/2021-4/2022): 58 months (1765 d) standard, 102.6 d PLD;
     pre-KAS250 (8/2018-7/2019) sensitivity: 61 months (1857 d) standard
  3. Donor all-cause mortality HR: base case 1.0 (Muzaale/US), sensitivity 1.30 (Mjøen)
  4. ESRD Weibull shape k and PSA sigma fit via cloglog regression to Massie 2017's
     5 published cumulative-incidence curves (median, IQR, 1st/99th pct)
  5. Race-specific waitlist/post-tx parameters from SRTR 2023 ADR
  6. wl_listing_prob=0.15 calibrated from USRDS 2025 ADR ESRD Ch.7 Figs 7.13/7.15/7.17
  7. Dialysis mortality age-stratified from USRDS 2023 ADR Ch.6 (Fig 6.5b, 2019)
  8. Prior-donor post-ESRD pathway from Muzaale 2016 (direct + calibrated inputs)
  9. Waitlist mortality age-stratified from SRTR 2023 ADR KI 25 (2019)
 10. Non-donor ESRD age gradient from Grams 2016 by race and sex (ages 20/40/60)

Output:
  data/processed/params.json
"""

import sys
import json
import math
import numpy as np
from pathlib import Path

# Allow running from repo root or src/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (DATA_PROC, save_params, _hardcoded_base_params, massie_weibull_fit,
                   load_life_table, dialysis_annual_mort, esrd_onset_mortality_curve)


def load_json(path: Path) -> dict:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def calibrate_donor_esrd_pathway(params: dict, m16: dict) -> dict:
    """
    Fit donor_dialysis_mort_mult and donor_wl_listing_prob to Muzaale 2016.

    For each candidate multiplier on a grid, the listing probability is solved
    in closed form so the cohort listed by 12 months matches the observed
    43.4% (preemptive + dialysis-year-1 survivors × listing prob); the
    multiplier minimizing squared error against the observed 1/3/5/10-yr
    mortality after ESRD is kept. Cohort enters ESRD at the study's median
    age (50). An absolute fit is used rather than one relative to Muzaale's
    matched non-donors: their comparator had 39% LDKT and a 21.5-month median
    wait, so the general-path mismatch reflects transplant access, not a
    dialysis-era effect that would cancel (see 08_calibration.py, CHECK 3).
    """
    lt  = load_life_table(DATA_PROC / "lifetable_combined_2021.csv")
    age = int(m16.get("median_age_at_esrd", 50))
    obs = {int(k): v for k, v in m16.get(
        "mort_after_esrd_donor", {"1": 0.091, "3": 0.202, "5": 0.255, "10": 0.324}).items()}
    listed_12mo = m16.get("listed_12mo_donor", 0.434)
    pre = params["donor_esrd_preemptive_prob"]
    years = sorted(obs)

    best = None
    for mult in np.arange(0.05, 2.5001, 0.005):
        p = dict(params, donor_dialysis_mort_mult=float(mult))
        dm1 = dialysis_annual_mort(p, age, first_year=True, donor=True, q_floor=lt[age])
        p["donor_wl_listing_prob"] = float(np.clip(
            (listed_12mo - pre) / ((1 - pre) * (1 - dm1)), 0.0, 1.0))
        curve = esrd_onset_mortality_curve(p, age, True, years, lt)
        sse = sum((curve[t] - obs[t]) ** 2 for t in years)
        if best is None or sse < best[0]:
            best = (sse, p["donor_dialysis_mort_mult"], p["donor_wl_listing_prob"])

    _, mult, listing = best
    return {"donor_dialysis_mort_mult": round(mult, 3),
            "donor_wl_listing_prob":    round(listing, 4)}


def main():
    print("=== 05_assemble_parameters.py ===\n")

    # Load component files
    lit      = load_json(DATA_PROC / "literature_params.json")
    usrds    = load_json(DATA_PROC / "usrds_params.json")
    usrds_e7 = load_json(DATA_PROC / "usrds_esrd7_params.json")
    usrds_mort = load_json(DATA_PROC / "usrds_mortality_params.json")
    srtr     = load_json(DATA_PROC / "srtr_params.json")

    muz = lit.get("muzaale2014", {})
    gr  = lit.get("grams2016", {})
    mas = lit.get("massie2017", {})
    mj  = lit.get("mjoeen2014", {})
    w   = lit.get("wainright2017", {})

    params = {}
    _fallback = _hardcoded_base_params()

    # ── ESRD RISK ─────────────────────────────────────────────────────────
    # Donor values: Muzaale 2014
    params["esrd_15yr_donor_overall"] = (
        muz.get("esrd_15yr_donor_overall_per10k", 30.8) / 10_000
    )
    params["esrd_15yr_donor_black"] = (
        muz.get("esrd_15yr_donor_black_per10k", 74.7) / 10_000
    )
    params["esrd_15yr_donor_white"] = (
        muz.get("esrd_15yr_donor_white_per10k", 22.7) / 10_000
    )

    # Non-donor baseline:
    # Overall: Muzaale 2014 matched controls
    params["esrd_15yr_nondonor"] = (
        muz.get("esrd_15yr_nondonor_matched_per10k", 3.9) / 10_000
    )
    # Race-stratified: Grams 2016 sex-averaged (RECONCILIATION DECISION 1)
    params["esrd_15yr_nondonor_black"] = (
        gr.get("esrd_15yr_nodonate_black_sexavg_pct", 0.195) / 100
    )
    params["esrd_15yr_nondonor_white"] = (
        gr.get("esrd_15yr_nodonate_white_sexavg_pct", 0.050) / 100
    )

    # Weibull shape k and its PSA log-sigma: fit (not assumed) via cloglog
    # regression to Massie 2017's five published cumulative-incidence curves
    # (median, IQR, 1st/99th pct) — see massie_weibull_fit() in utils.py.
    # (RECONCILIATION DECISION 4)
    params["weibull_shape"], params["weibull_shape_log_sigma"] = massie_weibull_fit(mas)

    # Hazard ratios within donor cohort (Massie 2017)
    params["hr_black_race"]             = mas.get("hr_black_race", 2.96)
    params["hr_male_sex"]               = mas.get("hr_male_sex", 1.88)
    params["hr_age_per_decade_nonblack"] = mas.get("hr_age_per_10yr_nonblack", 1.40)
    # Non-donor age gradient: Grams 2016 projections at ages 20/40/60 by race
    # and sex (RECONCILIATION DECISION 10), replacing a borrowed 1.40/decade
    # proxy that overstated Black non-donor risk away from age 40.
    for cell, by_age in gr.get("esrd_15yr_nodonate_by_age_pct", {}).items():
        for a, pct in by_age.items():
            params[f"grams_nondonor_15yr_{cell}_age{a}"] = pct / 100
    for k, v in _fallback.items():
        if k.startswith("grams_nondonor_15yr_"):
            params.setdefault(k, v)

    # ── DIALYSIS / ESRD SURVIVAL ──────────────────────────────────────────
    # Age-stratified (RECONCILIATION DECISION 7): USRDS 2023 ADR Fig 6.5b,
    # 2019 adjusted; first-year multiplier from Fig 6.4. Replaces the former
    # all-ages 22%/17% rates, which are dominated by older patients.
    for band in ("age1844", "age4564", "age6574", "age75p"):
        key = f"dialysis_mort_{band}"
        params[key] = usrds_mort.get(key, _fallback[key])
    params["dialysis_yr1_mult"]   = usrds_mort.get("dialysis_yr1_mult",
                                                   _fallback["dialysis_yr1_mult"])
    params["dialysis_mort_scale"] = 1.0   # uniform sensitivity scale

    # ── WAITLIST OUTCOMES ─────────────────────────────────────────────────
    # Key uses 2023 suffix; 2022 fallback retained for backwards compatibility
    params["wl_mort_per_100py"] = srtr.get(
        "pretx_mort_per_100py_overall_2023",
        srtr.get("pretx_mort_per_100py_overall_2022", 5.0)
    )
    params["wl_mort_black_per_100py"]    = srtr.get("pretx_mort_per_100py_black", 4.62)
    params["wl_mort_white_per_100py"]    = srtr.get("pretx_mort_per_100py_white", 5.71)
    # By age (KI 25), 2019 to match the USRDS dialysis year — used by the model
    # (RECONCILIATION DECISION 9); the all-ages rates above remain for the
    # SRTR KI 22 calibration check in 08_calibration.py.
    for band in ("age1834", "age3549", "age5064", "age65p"):
        key = f"wl_mort_per_100py_{band}"
        params[key] = srtr.get(f"pretx_mort_per_100py_{band}_2019", _fallback[key])
    # Race ÷ overall (KI 26 / KI 24, 2019), applied to the age table in race
    # subgroups, assuming race and age act multiplicatively
    _overall_19 = srtr.get("pretx_mort_per_100py_overall_2019", 4.78)
    params["wl_mort_race_mult"] = 1.0
    for race, fb in (("black", 4.47), ("white", 5.57)):
        params[f"wl_mort_race_mult_{race}"] = round(
            srtr.get(f"pretx_mort_per_100py_{race}_2019", fb) / _overall_19, 4)
    params["wl_removal_rate_yr"]         = srtr.get("wl_annual_removal_competing", 0.1260)

    # Wait times: post-KAS250 figures (RECONCILIATION DECISION 2)
    params["wl_std_mean_days"] = srtr.get("wl_std_mean_days", 1765)
    # Punjala 2024 reports PLD median; convert to mean for exponential model (mean = median/ln2)
    params["wl_pld_mean_days"] = w.get("pld_mwt_days_post_kas", 102.6) / math.log(2)
    # Sensitivity scenarios
    params["wl_std_mean_days_prekas250"] = srtr.get("wl_std_mean_days_prekas250", 1857)
    params["wl_pld_mean_days_from_activation"] = w.get("pld_mwt_from_activation", 23.0) / math.log(2)

    # Conditional per-cycle listing probability, calibrated from
    # USRDS 2025 ADR ESRD Ch.7 Figs 7.13/7.15/7.17 (RECONCILIATION DECISION 6)
    params["wl_listing_prob"] = usrds_e7.get("wl_listing_prob", 0.15)

    # Preemptive transplant listing probability at ESRD onset (USRDS 2025 Fig 7.13)
    params["esrd_preemptive_prob_std"] = usrds_e7.get(
        "usrds_esrd7_fig713_preemptive_listing_overall_2024", 0.058)
    params["esrd_preemptive_prob_pld"] = usrds_e7.get(
        "usrds_esrd7_fig713_preemptive_listing_age1844_2024", 0.094)
    params["wl_listing_prob_sens_low"]  = usrds_e7.get("wl_listing_prob_sens_low",  0.05)
    params["wl_listing_prob_sens_high"] = usrds_e7.get("wl_listing_prob_sens_high", 0.30)

    # ── POST-TRANSPLANT OUTCOMES ──────────────────────────────────────────
    # Age-stratified annual mortality from SRTR 2023 DDKT patient survival
    params["posttx_annual_mort_age1834"] = srtr.get("posttx_dd_annual_mort_age1834", 0.009)
    params["posttx_annual_mort_age3549"] = srtr.get("posttx_dd_annual_mort_age3549", 0.018)
    params["posttx_annual_mort_age5064"] = srtr.get("posttx_dd_annual_mort_age5064", 0.039)
    params["posttx_annual_mort_age65p"]  = srtr.get("posttx_dd_annual_mort_age65p",  0.068)

    # Overall: age-weighted average of SRTR age-strata
    # Weights from SRTR 2023 ADR Figure KI 1 (incident kidney transplant
    # recipients by age): 18-34 ≈ 9%, 35-49 ≈ 28%, 50-64 ≈ 41%, 65+ ≈ 22%.
    # Rounded to [0.10, 0.30, 0.40, 0.20] for conservatism (shifts weight to
    # younger, healthier recipients, slightly understating average mortality).
    _age_weights = [0.10, 0.30, 0.40, 0.20]
    _age_keys    = ["posttx_annual_mort_age1834", "posttx_annual_mort_age3549",
                    "posttx_annual_mort_age5064", "posttx_annual_mort_age65p"]
    params["posttx_annual_mort"] = sum(
        params[k] * w for k, w in zip(_age_keys, _age_weights)
    )

    # Race-stratified: SRTR 2023 DDKT patient survival by race
    params["posttx_annual_mort_black"] = srtr.get(
        "posttx_dd_annual_mort_black",
        srtr.get("posttx_annual_mort_black", 0.035)
    )
    params["posttx_annual_mort_white"] = srtr.get(
        "posttx_dd_annual_mort_white",
        srtr.get("posttx_annual_mort_white", 0.038)
    )
    # Age-stratified post-year-1 graft failure (SRTR 2023 ADR Figure KI 53).
    # Prefer the value already computed in 03_download_srtr.py; if a key is
    # missing (e.g. stale/partial srtr_params.json), recompute it here from
    # the same 5-yr/1-yr graft survival inputs rather than falling back to
    # an arbitrary flat number. A single flat rate materially understated
    # every age band (previously 0.025 vs. a true ~3.3%-8.8%/yr range).
    _gf_one_yr = srtr.get("ddkt_graft_1yr_assumed", 0.955)
    _gf_fallback_5yr = {"age1834": 0.822, "age3549": 0.835,
                         "age5064": 0.768, "age65p":  0.661}
    _gf_keys = []
    for _band, _fb5 in _gf_fallback_5yr.items():
        _key = f"graft_annual_fail_postyear1_{_band}"
        _default = round(1 - (srtr.get(f"ddkt_graft_5yr_{_band}", _fb5)
                               / _gf_one_yr) ** 0.25, 4)
        params[_key] = srtr.get(_key, _default)
        _gf_keys.append(_key)
    params["graft_annual_fail_postyear1"] = round(
        sum(params[k] * w for k, w in zip(_gf_keys, _age_weights)), 4
    )

    # Age-stratified LDKT annual mortality (SRTR 2023 ADR Figure KI 76)
    for band, fallback in [("age1834", 0.979), ("age3549", 0.961),
                           ("age5064", 0.917), ("age65p",  0.819)]:
        surv5 = srtr.get(f"posttx_ld_5yr_patient_surv_{band}", fallback)
        params[f"posttx_ld_annual_mort_{band}"] = round(1 - surv5 ** (1 / 5), 4)
    _ld_keys = ["posttx_ld_annual_mort_age1834", "posttx_ld_annual_mort_age3549",
                "posttx_ld_annual_mort_age5064", "posttx_ld_annual_mort_age65p"]
    params["posttx_ld_annual_mort"] = round(
        sum(params[k] * w for k, w in zip(_ld_keys, _age_weights)), 4
    )

    # Graft survival by age/donor type
    params["ddkt_graft_5yr_age1834"] = srtr.get("ddkt_graft_5yr_age1834", 0.814)
    params["ddkt_graft_5yr_age3564"] = srtr.get("ddkt_graft_5yr_age3564", 0.760)
    params["ddkt_graft_5yr_age65p"]  = srtr.get("ddkt_graft_5yr_age65p",  0.678)
    params["ldkt_graft_5yr_age1834"] = srtr.get("ldkt_graft_5yr_age1834", 0.900)
    params["ldkt_graft_5yr_age3564"] = srtr.get("ldkt_graft_5yr_age3564", 0.848)
    params["ldkt_graft_5yr_age65p"]  = srtr.get("ldkt_graft_5yr_age65p",  0.808)

    # ── BACKGROUND MORTALITY ──────────────────────────────────────────────
    # Donor all-cause mortality HR, modeled as time-varying rather than a
    # single fixed number: short/medium follow-up studies (Segev 2010, Garg
    # 2012, Berger 2011, and the O'Keeffe 2018 pooled estimate that combines
    # them) show no excess donor mortality, while Mjøen 2014 (median ~15 yr
    # follow-up) found donor/control survival curves separate only after
    # ~10 years, reaching HR=1.30 (95% CI 1.11-1.52) by ~25 years. These are
    # not conflicting estimates of one static HR but evidence of a
    # non-proportional (time-varying) hazard, so the disagreement is resolved
    # structurally (donor_mort_hr_at() in utils.py) rather than by pooling.
    # (RECONCILIATION DECISION 3)
    params["donor_mort_hr_early"]   = 1.0
    params["donor_mort_hr_late"]    = mj.get("hr_all_cause_mortality", 1.30)
    params["donor_mort_hr_t_start"] = 10.0
    params["donor_mort_hr_t_end"]   = 15.0

    # ── PRIOR-DONOR POST-ESRD PATHWAY (RECONCILIATION DECISION 8) ─────────
    # Donors who develop ESRD are followed closely, healthier, and listed
    # earlier than the general ESRD population (Muzaale 2016). Direct inputs:
    # preemptive fraction and post-Tx mortality HR. Calibrated inputs: listing
    # rate (to 43.4% listed by 12 mo) and dialysis-mortality multiplier (to the
    # 1/3/5/10-yr mortality-after-ESRD curve).
    m16 = lit.get("muzaale2016", {})
    params["donor_esrd_preemptive_prob"] = (
        m16.get("n_preemptive", 21) / m16.get("n_donors_esrd", 99))
    params["donor_posttx_mort_mult"] = m16.get("hr_posttx_mort_adjusted", 0.7)
    params.update(calibrate_donor_esrd_pathway(params, m16))

    # ── METADATA ──────────────────────────────────────────────────────────
    params["_sources"] = {
        "esrd_donor_risk":        "Muzaale 2014 JAMA 311:579",
        "esrd_nondonor_baseline": "Grams 2016 NEJM 374:411 (race-stratified); "
                                  "Muzaale 2014 (overall)",
        "esrd_hr_within_donors":  "Massie 2017 JASN 28:2749",
        "weibull_shape":          "Fit via cloglog regression to Massie 2017's published "
                                  "cumulative-incidence curves (median, IQR, 1st/99th pct)",
        "dialysis_survival":      "USRDS 2023 ADR ESRD Ch.6 Fig 6.5b (2019 adjusted, by age); "
                                  "first-year multiplier Fig 6.4",
        "donor_esrd_pathway":     "Muzaale 2016 Transplantation 100:1306 (preemptive 21/99; "
                                  "post-Tx HR 0.7; listing and dialysis multiplier calibrated)",
        "waitlist_outcomes":      "SRTR 2023 ADR",
        "wl_listing_prob":        "USRDS 2025 ADR ESRD Vol. Ch.7 Figs 7.13/7.15/7.17",
        "pld_wait_time":          "Wainright 2017 AJT 17:1103; UNOS ATC abstract 2015",
        "posttx_survival":        "SRTR 2023 ADR DDKT patient survival (age-stratified)",
        "donor_mort_hr_early":    "Muzaale 2014; Segev 2010 JAMA 303:959; "
                                  "O'Keeffe 2018 meta-analysis (PubMed 29379948)",
        "donor_mort_hr_late":     "Mjøen 2014 Kidney Int 86:162 (HR=1.30, 95% CI "
                                  "1.11-1.52; survival curves separate after ~10 yr)",
        "life_tables":            "CDC NVSR Vol 72 No 12 (Nov 2023) — 2021 US life tables",
    }
    params["_reconciliation_notes"] = [
        "White non-donor baseline uses Grams 2016 (0.05%) not Muzaale (0 events/unstable)",
        "Wait times use Punjala 2024 national mean waiting time at transplant: "
        "1765d standard post-KAS250 (100d PLD, Wainright 2017); 1857d pre-KAS250 sensitivity",
        "Donor all-cause mortality HR=1.0 base case per US evidence; "
        "1.30 as sensitivity (Mjøen)",
        "Weibull shape k and PSA sigma fit via cloglog regression to Massie 2017's 5 "
        "published cumulative-incidence curves (median, IQR, 1st/99th pct), not assumed "
        "(see massie_weibull_fit in utils.py)",
        "Post-tx mortality age-stratified from SRTR 2023 DDKT 5-yr patient survival",
        "wl_listing_prob=0.15: back-calculated from USRDS 2025 Fig 7.15 3-yr CIF=12% "
        "(general ESRD p≈0.065/yr), scaled for donor-like 18-44 cohort via Figs 7.13+7.17 "
        "(post-dialysis yr1 listing ≈9.2%, conditional p≈0.17/yr); conservative base 0.15; "
        "sensitivity 0.05-0.30. Replaces prior placeholder 0.75.",
        "Dialysis mortality age-stratified (USRDS 2023 Fig 6.5b, 2019 pre-COVID) instead of "
        "all-ages 22%/17%; all post-ESRD state mortality floored at life-table qx.",
        "Waitlist mortality by age from SRTR KI 25 (2019), race ratio from KI 26; all "
        "age-banded inputs interpolated log-linearly between band midpoints.",
        "Donor arm uses Muzaale 2016 prior-donor inputs after ESRD; voucher and ESRD-"
        "conditional priority arms do not (they isolate priority access, not donor health).",
    ]

    save_params(params)

    # Print summary for verification
    print("\nAssembled parameters:")
    print(f"  Donor ESRD 15yr (overall):    {params['esrd_15yr_donor_overall']:.4%}")
    print(f"  Nondonor ESRD 15yr (overall): {params['esrd_15yr_nondonor']:.4%}")
    print(f"  RR:                           {params['esrd_15yr_donor_overall']/params['esrd_15yr_nondonor']:.1f}×")
    print(f"  Nondonor ESRD 15yr (black):   {params['esrd_15yr_nondonor_black']:.4%}")
    print(f"  Nondonor ESRD 15yr (white):   {params['esrd_15yr_nondonor_white']:.4%}")
    print(f"  Weibull shape (fit):          {params['weibull_shape']:.4f}")
    print(f"  Weibull shape log-sigma (fit):{params['weibull_shape_log_sigma']:.4f}")
    print(f"  Dialysis mortality by age:    "
          + " / ".join(f"{params[f'dialysis_mort_{b}']:.1%}"
                       for b in ("age1844", "age4564", "age6574", "age75p"))
          + f"  (yr-1 ×{params['dialysis_yr1_mult']:.2f})")
    print(f"  Donor | ESRD: preemptive {params['donor_esrd_preemptive_prob']:.1%}, "
          f"listing {params['donor_wl_listing_prob']:.3f}/yr, "
          f"dialysis mort ×{params['donor_dialysis_mort_mult']:.2f}, "
          f"post-Tx mort ×{params['donor_posttx_mort_mult']:.2f}")
    print(f"  Waitlist mort (overall):      {params['wl_mort_per_100py']:.1f}/100 PY")
    print(f"  Std wait (mean days):         {params['wl_std_mean_days']}")
    print(f"  PLD wait (mean days):         {params['wl_pld_mean_days']:.1f}")
    print(f"  WL listing prob/yr:           {params['wl_listing_prob']:.0%}")
    print(f"  Post-tx annual mort (overall):{params['posttx_annual_mort']:.1%}")
    print(f"    age 18-34:                  {params['posttx_annual_mort_age1834']:.1%}")
    print(f"    age 35-49:                  {params['posttx_annual_mort_age3549']:.1%}")
    print(f"    age 50-64:                  {params['posttx_annual_mort_age5064']:.1%}")
    print(f"    age 65+:                    {params['posttx_annual_mort_age65p']:.1%}")
    print(f"  Donor mort HR (early, <={params['donor_mort_hr_t_start']:.0f}yr): {params['donor_mort_hr_early']:.2f}")
    print(f"  Donor mort HR (late, >={params['donor_mort_hr_t_end']:.0f}yr):  {params['donor_mort_hr_late']:.2f}")
    print("\nDone.")


if __name__ == "__main__":
    main()
