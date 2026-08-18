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
  4. ESRD Weibull shape k=1.5 (supported by Massie 2017 20-yr data)
  5. Race-specific waitlist/post-tx parameters from SRTR 2023 ADR
  6. wl_listing_prob=0.15 calibrated from USRDS 2025 ADR ESRD Ch.7 Figs 7.13/7.15/7.17

Output:
  data/processed/params.json
"""

import sys
import json
from pathlib import Path

# Allow running from repo root or src/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import DATA_PROC, save_params, _hardcoded_base_params


def load_json(path: Path) -> dict:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def main():
    print("=== 05_assemble_parameters.py ===\n")

    # Load component files
    lit      = load_json(DATA_PROC / "literature_params.json")
    usrds    = load_json(DATA_PROC / "usrds_params.json")
    usrds_e7 = load_json(DATA_PROC / "usrds_esrd7_params.json")
    srtr     = load_json(DATA_PROC / "srtr_params.json")

    muz = lit.get("muzaale2014", {})
    gr  = lit.get("grams2016", {})
    mas = lit.get("massie2017", {})
    mj  = lit.get("mjoeen2014", {})
    w   = lit.get("wainright2017", {})

    params = {}

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

    # Weibull shape: k=1.5 supported by Massie 2017 20-yr data (RECONCILIATION DECISION 4)
    params["weibull_shape"] = 1.5

    # Hazard ratios within donor cohort (Massie 2017)
    params["hr_black_race"]             = mas.get("hr_black_race", 2.96)
    params["hr_male_sex"]               = mas.get("hr_male_sex", 1.88)
    params["hr_age_per_decade_nonblack"] = mas.get("hr_age_per_10yr_nonblack", 1.40)
    params["hr_age_per_decade_nondonor"] = mas.get("hr_age_per_10yr_nondonor", 1.40)

    # ── DIALYSIS / ESRD SURVIVAL ──────────────────────────────────────────
    params["dialysis_1yr_mort"]        = usrds.get("hd_1yr_mortality", 0.22)
    params["dialysis_annual_mort"]     = usrds.get("hd_annual_mort_postyear1", 0.17)

    # ── WAITLIST OUTCOMES ─────────────────────────────────────────────────
    # Key uses 2023 suffix; 2022 fallback retained for backwards compatibility
    params["wl_mort_per_100py"] = srtr.get(
        "pretx_mort_per_100py_overall_2023",
        srtr.get("pretx_mort_per_100py_overall_2022", 5.0)
    )
    params["wl_mort_black_per_100py"]    = srtr.get("pretx_mort_per_100py_black", 4.62)
    params["wl_mort_white_per_100py"]    = srtr.get("pretx_mort_per_100py_white", 5.71)
    params["wl_removal_rate_yr"]         = srtr.get("wl_annual_removal_competing", 0.1260)

    # Wait times: post-KAS250 figures (RECONCILIATION DECISION 2)
    params["wl_std_median_days"] = srtr.get("wl_std_median_days", 1765)
    params["wl_pld_median_days"] = w.get("pld_mwt_days_post_kas", 102.6)
    # Sensitivity scenarios
    params["wl_std_median_days_prekas250"] = srtr.get("wl_std_median_days_prekas250", 1857)
    params["wl_pld_median_days_from_activation"] = w.get("pld_mwt_from_activation", 23.0)

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
    # HR base case = 1.0 (Muzaale/Segev/Grams meta-analysis)
    # (RECONCILIATION DECISION 3)
    params["donor_mort_hr"]         = 1.0
    params["donor_mort_hr_mjoeen"]  = mj.get("hr_all_cause_mortality", 1.30)

    # ── METADATA ──────────────────────────────────────────────────────────
    params["_sources"] = {
        "esrd_donor_risk":        "Muzaale 2014 JAMA 311:579",
        "esrd_nondonor_baseline": "Grams 2016 NEJM 374:411 (race-stratified); "
                                  "Muzaale 2014 (overall)",
        "esrd_hr_within_donors":  "Massie 2017 JASN 28:2749",
        "weibull_shape":          "Calibrated to Massie 2017 20-yr cumulative incidence",
        "dialysis_survival":      "USRDS 2023/2024 ADR (hardcoded from published tables)",
        "waitlist_outcomes":      "SRTR 2023 ADR",
        "wl_listing_prob":        "USRDS 2025 ADR ESRD Vol. Ch.7 Figs 7.13/7.15/7.17",
        "pld_wait_time":          "Wainright 2017 AJT 17:1103; UNOS ATC abstract 2015",
        "posttx_survival":        "SRTR 2023 ADR DDKT patient survival (age-stratified)",
        "donor_mort_hr_base":     "Muzaale 2014; Segev 2010 JAMA 303:959; "
                                  "Grams 2018 meta-analysis (PubMed 29379948)",
        "donor_mort_hr_sens":     "Mjøen 2014 Kidney Int 86:162",
        "life_tables":            "CDC NVSR Vol 72 No 12 (Nov 2023) — 2021 US life tables",
    }
    params["_reconciliation_notes"] = [
        "White non-donor baseline uses Grams 2016 (0.05%) not Muzaale (0 events/unstable)",
        "Wait times use Punjala 2024 national mean waiting time at transplant: "
        "1765d standard post-KAS250 (100d PLD, Wainright 2017); 1857d pre-KAS250 sensitivity",
        "Donor all-cause mortality HR=1.0 base case per US evidence; "
        "1.30 as sensitivity (Mjøen)",
        "Weibull k=1.5 supported by Massie 2017 20-yr data showing accelerating hazard",
        "Post-tx mortality age-stratified from SRTR 2023 DDKT 5-yr patient survival",
        "wl_listing_prob=0.15: back-calculated from USRDS 2025 Fig 7.15 3-yr CIF=12% "
        "(general ESRD p≈0.065/yr), scaled for donor-like 18-44 cohort via Figs 7.13+7.17 "
        "(post-dialysis yr1 listing ≈9.2%, conditional p≈0.17/yr); conservative base 0.15; "
        "sensitivity 0.05-0.30. Replaces prior placeholder 0.75.",
    ]

    save_params(params)

    # Print summary for verification
    print("\nAssembled parameters:")
    print(f"  Donor ESRD 15yr (overall):    {params['esrd_15yr_donor_overall']:.4%}")
    print(f"  Nondonor ESRD 15yr (overall): {params['esrd_15yr_nondonor']:.4%}")
    print(f"  RR:                           {params['esrd_15yr_donor_overall']/params['esrd_15yr_nondonor']:.1f}×")
    print(f"  Nondonor ESRD 15yr (black):   {params['esrd_15yr_nondonor_black']:.4%}")
    print(f"  Nondonor ESRD 15yr (white):   {params['esrd_15yr_nondonor_white']:.4%}")
    print(f"  Weibull shape:                {params['weibull_shape']}")
    print(f"  Dialysis 1yr mortality:       {params['dialysis_1yr_mort']:.0%}")
    print(f"  Waitlist mort (overall):      {params['wl_mort_per_100py']:.1f}/100 PY")
    print(f"  Std wait (median days):       {params['wl_std_median_days']}")
    print(f"  PLD wait (median days):       {params['wl_pld_median_days']:.1f}")
    print(f"  WL listing prob/yr:           {params['wl_listing_prob']:.0%}")
    print(f"  Post-tx annual mort (overall):{params['posttx_annual_mort']:.1%}")
    print(f"    age 18-34:                  {params['posttx_annual_mort_age1834']:.1%}")
    print(f"    age 35-49:                  {params['posttx_annual_mort_age3549']:.1%}")
    print(f"    age 50-64:                  {params['posttx_annual_mort_age5064']:.1%}")
    print(f"    age 65+:                    {params['posttx_annual_mort_age65p']:.1%}")
    print(f"  Donor mort HR (base):         {params['donor_mort_hr']:.2f}")
    print(f"  Donor mort HR (Mjøen):        {params['donor_mort_hr_mjoeen']:.2f}")
    print("\nDone.")


if __name__ == "__main__":
    main()
