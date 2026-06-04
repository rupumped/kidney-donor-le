"""
04_literature_params.py
───────────────────────
Encodes parameters extracted from primary literature during the parameter
confirmation process. Each value is documented with its source, table/figure
reference, and any caveats noted during confirmation.

Papers covered:
  - Muzaale 2014  (JAMA 311:579)       — primary ESRD risk study
  - Grams 2016    (NEJM 374:411)       — donor-candidate ESRD baseline tool
  - Massie 2017   (JASN 28:2749)       — post-donation ESRD hazard ratios
  - Mjøen 2014    (Kidney Int 86:162)  — Norwegian cohort (sensitivity upper bound)
  - Wainright 2017 (AJT 17:1103)      — PLD wait time post-KAS
  - Ibrahim 2009  (NEJM 360:459)       — long-term donor outcomes
  - Segev 2010    (JAMA 303:959)       — perioperative mortality and long-term survival
  - Grams 2018    (meta-analysis, PubMed 29379948) — all-cause mortality HR

Output:
  data/processed/literature_params.json
"""

import sys
import json
from pathlib import Path

# Allow running from repo root or src/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import DATA_PROC

LITERATURE_PARAMS = {

    # ════════════════════════════════════════════════════════════════════════
    # MUZAALE 2014 — Risk of ESRD Following Live Kidney Donation
    # JAMA 2014;311(6):579–586  PMCID: PMC4411956
    # Design: 96,217 US donors vs 20,024 matched healthy NHANES controls
    # Median follow-up: 7.6 yr donors, 15.0 yr controls
    # ════════════════════════════════════════════════════════════════════════
    "muzaale2014": {

        # Table 2 / PMC full text — 15-yr cumulative incidence
        "esrd_15yr_donor_overall_per10k":      30.8,   # 95% CI 24.3–38.5
        "esrd_15yr_nondonor_matched_per10k":    3.9,   # 95% CI 0.8–8.9
        "esrd_15yr_donor_black_per10k":        74.7,   # 95% CI 47.8–105.8
        "esrd_15yr_nondonor_black_per10k":     23.9,   # 95% CI 1.6–62.4
        "esrd_15yr_donor_white_per10k":        22.7,   # 95% CI 15.6–30.1
        "esrd_15yr_nondonor_white_per10k":      0.0,   # CAVEAT: zero events in white controls
                                                        # → use Grams 2016 for white baseline

        # Lifetime risk (to age 80) from Table 2
        "esrd_lifetime_donor_per10k":          90.0,
        "esrd_lifetime_nondonor_matched_per10k": 14.0,
        "esrd_lifetime_general_pop_per10k":   326.0,  # unscreened controls

        # Relative risk at 15 years
        "rr_15yr_overall":                      7.9,   # 30.8 / 3.9
        "rr_lifetime_overall":                  6.4,   # 90 / 14

        # CAVEAT on Muzaale RR:
        # The 8× RR vs MATCHED controls is the appropriate figure for modelling.
        # The unscreened general population is not the right comparator.

        # Study notes
        "_note_white_controls": (
            "Zero ESRD events in white non-donors at 15 yr, making the white "
            "non-donor figure (0/n ≈ 0.0) an unstable estimate. "
            "Use Grams 2016 NEJM for white donor-candidate baseline."
        ),
    },

    # ════════════════════════════════════════════════════════════════════════
    # GRAMS 2016 — Kidney-Failure Risk Projection for the Living Kidney-Donor Candidate
    # NEJM 2016;374(5):411–421
    # Design: meta-analysis of 7 general-population cohorts (n=4,933,314)
    #         calibrated to US ESRD incidence; compared to 52,998 US living donors
    # This is the appropriate source for pre-donation donor-CANDIDATE baseline risk
    # ════════════════════════════════════════════════════════════════════════
    "grams2016": {

        # 15-year ESRD risk projections in the ABSENCE of donation
        # for a 40-year-old with health characteristics similar to age-matched donors
        # NEJM Table 2 / text
        "esrd_15yr_nodonate_black_male_pct":   0.24,
        "esrd_15yr_nodonate_black_female_pct": 0.15,
        "esrd_15yr_nodonate_white_male_pct":   0.06,
        "esrd_15yr_nodonate_white_female_pct": 0.04,

        # Sex-averaged (simple mean) — used as model base case
        "esrd_15yr_nodonate_black_sexavg_pct": 0.195,  # (0.24+0.15)/2
        "esrd_15yr_nodonate_white_sexavg_pct": 0.050,  # (0.06+0.04)/2

        # Key finding: lifetime ESRD risk highest at youngest ages, esp. young Black
        # Post-donation 15-yr risk was 3.5–5.3× higher than these projections
        "donation_multiplier_vs_projection_range": [3.5, 5.3],

        "_note": (
            "This paper provides the correct donor-candidate baseline — not USRDS "
            "general population rates. The Grams tool is available at "
            "https://www.transplantmodels.com/esrdrisk/"
        ),
    },

    # ════════════════════════════════════════════════════════════════════════
    # MASSIE 2017 — Quantifying Postdonation Risk of ESRD in Living Kidney Donors
    # JASN 2017;28(9):2749–2755
    # Design: 133,824 US living donors 1987–2015, CMS linkage, Cox regression
    # ════════════════════════════════════════════════════════════════════════
    "massie2017": {

        # Hazard ratios for post-donation ESRD within the donor cohort
        "hr_black_race":                    2.96,   # 95% CI 2.25–3.89
        "hr_male_sex":                      1.88,   # 95% CI 1.50–2.35
        "hr_age_per_10yr_nonblack":         1.40,   # 95% CI 1.23–1.59

        # Cumulative incidence by time point (per 10,000 donors)
        "esrd_cum_incidence_5yr_per10k":    1.0,    # range 1–2
        "esrd_cum_incidence_10yr_per10k":   6.0,    # range 4–11
        "esrd_cum_incidence_15yr_per10k":  16.0,    # range 10–29
        "esrd_cum_incidence_20yr_per10k":  34.0,    # range 20–59

        "_note": (
            "The 20-year figure (34/10,000) supports using a Weibull shape > 1 "
            "(accelerating hazard). Under a constant-hazard model, 15-yr rate of "
            "16/10,000 implies 20-yr rate of ~21/10,000 — well below the observed 34. "
            "Use k=1.5 Weibull as base case."
        ),
    },

    # ════════════════════════════════════════════════════════════════════════
    # MJØEN 2014 — Long-term Risks for Kidney Donors
    # Kidney Int 2014;86(1):162–167
    # Design: 1,901 Norwegian donors (1963–2007) vs 32,621 matched controls
    # Median follow-up: 15.1 yr donors, 24.9 yr controls
    # ════════════════════════════════════════════════════════════════════════
    "mjoeen2014": {

        # Primary outcomes
        "hr_all_cause_mortality":           1.30,   # 95% CI 1.11–1.52
        "hr_cardiovascular_mortality":      1.40,   # 95% CI 1.03–1.91
        "hr_esrd":                         11.38,   # 95% CI 4.37–29.6

        # Absolute ESRD: 9 donors (0.47%) vs 22 non-donors (0.068%)
        "esrd_rate_donors_pct":             0.47,
        "esrd_rate_controls_pct":           0.068,

        # CAVEATS — this is the upper bound / contested estimate
        "_caveats": [
            "80% of donors were first-degree relatives of recipient, sharing genetic "
            "ESRD risk. This inflates the HR vs a matched-donor control design.",
            "Janki et al. (Eur J Epidemiol 2017) identified methodological issues "
            "suggesting risk overestimation in matched controls.",
            "US studies using NHANES-matched healthy controls (Muzaale 2014, Segev 2010) "
            "show HR for all-cause mortality ≈ 1.0.",
            "Mjøen HR=1.30 is used as the pessimistic scenario in sensitivity analyses, "
            "not as the base case.",
        ],
    },

    # ════════════════════════════════════════════════════════════════════════
    # WAINRIGHT 2017 — Impact of KAS on Prior Living Kidney Donors' Access
    # AJT 2017;17(4):1103–1111
    # Design: OPTN data, incident/prevalent PLDs pre-KAS vs post-KAS (2013–2015)
    # ════════════════════════════════════════════════════════════════════════
    "wainright2017": {

        "pld_mwt_days_post_kas":        102.6,
        "pld_mwt_days_pre_kas":          82.3,
        "pld_mwt_pvalue":                0.98,   # not significantly different

        # UNOS ATC abstract (2015) — similar time period
        "pld_mwt_days_ats_abstract":     93.0,

        # CJASN 2016 study (n=210 PLDs, Jan 2010–Jul 2015)
        # ~50% received kidney within 98 days of listing
        # ~50% from activation: 23 days median
        "pld_mwt_overall_listing":       98.0,
        "pld_mwt_from_activation":       23.0,

        # Key finding on inactive listing
        "pct_inactive_at_listing":       0.68,   # 68% listed as inactive initially
        "pct_inactive_90d_to_1yr":       0.18,
        "pct_inactive_over_1yr":         0.14,

        "_note": (
            "The 145-day figure in the original research plan was pre-KAS. "
            "Post-KAS confirmed values are 93–103 days from listing, ~23 days "
            "from activation. The administrative friction of inactive listing is "
            "a real-world parameter that should be captured in sensitivity analyses."
        ),
    },

    # ════════════════════════════════════════════════════════════════════════
    # SEGEV 2010 — Perioperative Mortality and Long-term Survival
    # JAMA 2010;303(10):959–966
    # Design: 80,347 US living kidney donors 1994–2009 vs NHANES controls
    # ════════════════════════════════════════════════════════════════════════
    "segev2010": {
        "perioperative_mortality_per_10k":  3.1,    # 31 deaths / 80,347 donors
        "hr_long_term_all_cause_mort":      0.85,   # donors BETTER than controls
                                                    # (healthy donor effect)
        "_note": (
            "Donors had LOWER long-term all-cause mortality than matched controls "
            "in this US study, reflecting healthy-donor selection. HR=1.0 is the "
            "appropriate base case; the Mjøen HR=1.30 is the pessimistic scenario."
        ),
    },

    # ════════════════════════════════════════════════════════════════════════
    # GRAMS 2018 META-ANALYSIS — Mid- and Long-term Health Risks
    # PubMed 29379948; Ann Intern Med 2018;168(4):276–284
    # 52 studies, 118,426 donors, 117,656 nondonors; avg follow-up 1–24 yr
    # ════════════════════════════════════════════════════════════════════════
    "grams2018_meta": {
        # Pooled all-cause mortality HR from Table 2, pooled estimate across
        # 9 studies with mortality outcome (subset of the 52 total studies).
        # HR = 0.984, 95% CI 0.743–1.302 (p=0.91 for heterogeneity).
        "hr_all_cause_mortality":       0.984,
        "hr_all_cause_mortality_ci_lo": 0.743,
        "hr_all_cause_mortality_ci_hi": 1.302,
        "hr_cvd":                  None,    # no significant increase found
        "hr_hypertension":         None,    # no significant increase found
        "hr_type2_diabetes":       None,    # no significant increase found
        "_finding": (
            "All-cause mortality HR 0.984 (95% CI 0.743–1.302) — not significantly "
            "elevated vs nondonor populations. No evidence of higher risk for CVD, "
            "hypertension, T2DM, or adverse psychosocial outcomes. Supports HR=1.0 "
            "base case; CI used to parameterise PSA log-normal distribution."
        ),
    },
}


def main():
    print("=== 04_literature_params.py ===\n")
    print("Encoding literature-derived parameters...\n")

    out = DATA_PROC / "literature_params.json"
    with open(out, "w") as f:
        json.dump(LITERATURE_PARAMS, f, indent=2)
    print(f"Saved → {out}")

    # Print summary
    print()
    muz = LITERATURE_PARAMS["muzaale2014"]
    gr  = LITERATURE_PARAMS["grams2016"]
    mas = LITERATURE_PARAMS["massie2017"]
    mj  = LITERATURE_PARAMS["mjoeen2014"]
    w   = LITERATURE_PARAMS["wainright2017"]

    print("Summary of key literature values:")
    print(f"  Muzaale 2014 — 15-yr ESRD per 10k: donor={muz['esrd_15yr_donor_overall_per10k']}, "
          f"matched control={muz['esrd_15yr_nondonor_matched_per10k']}, "
          f"RR={muz['rr_15yr_overall']:.1f}×")
    print(f"  Grams 2016   — 15-yr baseline (40yr, no donate): "
          f"Black M={gr['esrd_15yr_nodonate_black_male_pct']}%, "
          f"White M={gr['esrd_15yr_nodonate_white_male_pct']}%")
    print(f"  Massie 2017  — 20-yr donor ESRD per 10k: {mas['esrd_cum_incidence_20yr_per10k']} "
          f"(supports accelerating Weibull hazard)")
    print(f"  Mjøen 2014   — HR all-cause mort: {mj['hr_all_cause_mortality']} "
          f"(upper bound / contested)")
    print(f"  Wainright 2017 — PLD median wait post-KAS: {w['pld_mwt_days_post_kas']} days")
    print()
    print("Note: Muzaale white non-donor rate is 0 (zero events).")
    print("      Use Grams 2016 for white donor-candidate baseline in race-stratified runs.")
    print("\nDone.")


if __name__ == "__main__":
    main()
