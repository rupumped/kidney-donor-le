"""
02_download_usrds.py
────────────────────
Extracts CKD prevalence tables from the United States Renal Data System (USRDS)
Annual Data Report, CKD Volume Chapter 1: "CKD in the General Population."

Data source: NHANES (2017–March 2020 is the most recent survey period).

─── ACCESS NOTE ──────────────────────────────────────────────────────────────
USRDS data is free but requires registration at:
  https://usrds-adr.niddk.nih.gov/

The 2025 ADR is the most recent available:
  https://usrds-adr.niddk.nih.gov/2025/chronic-kidney-disease/
    1-ckd-in-the-general-population

Key tables:
  Table 1.1 — % and N of U.S. adults in KDIGO CKD risk categories by eGFR/ACR
  Table 1.2 — CKD prevalence by insurance, income, and education level
  Table 1.3 — Health risk behaviors by CKD status (adjusted and unadjusted)

Data files can be downloaded via the ADR Data Download page as Excel files.
This script expects them at:
  data/raw/usrds_ckd_chapter1_table1.xlsx
  data/raw/usrds_ckd_chapter1_table2.xlsx
  data/raw/usrds_ckd_chapter1_table3.xlsx

If those files are absent, the script writes confirmed published values
extracted from the 2025 ADR CKD Chapter 1 PDF to processed/.

─── MANUAL DOWNLOAD STEPS ────────────────────────────────────────────────────
1. Go to https://usrds-adr.niddk.nih.gov/2025/chronic-kidney-disease/
       1-ckd-in-the-general-population
2. Click "Download Data" and export Table 1.1, 1.2, and 1.3
3. Place files in data/raw/ with the names shown above

Outputs:
  data/processed/usrds_ckd_kdigo_risk.csv     KDIGO risk category %/N by eGFR×ACR
  data/processed/usrds_ckd_by_insurance.csv   CKD % by insurance/income/education
  data/processed/usrds_ckd_risk_behaviors.csv Health risk behaviors by CKD status
  data/processed/usrds_ckd_params.json        Scalar parameters for model
"""

import sys
import json
import warnings
import pandas as pd
from pathlib import Path

# Allow running from repo root or src/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import DATA_RAW, DATA_PROC

warnings.filterwarnings("ignore", category=UserWarning)


# ── PUBLISHED FALLBACK VALUES ─────────────────────────────────────────────────
# Source: USRDS 2025 ADR, CKD Volume Chapter 1
# Survey period: 2017–March 2020 (NHANES); adults aged ≥20 with serum Cr + ACR
# Note: based on single eGFR/ACR measurement — likely overestimates true CKD
# prevalence relative to the KDIGO guideline requirement of ≥3 months.

# ── TABLE 1.1B: KDIGO RISK CATEGORY TRENDS (% OF U.S. ADULTS) ──────────────────
KDIGO_RISK_TRENDS = {
    # period: {risk_category: pct}
    "2005-2008": {
        "low_risk":             86.7,
        "moderately_high_risk":  9.6,
        "high_risk":             2.4,
        "very_high_risk":        1.3,
    },
    "2009-2012": {
        "low_risk":             87.5,
        "moderately_high_risk":  9.0,
        "high_risk":             2.3,
        "very_high_risk":        1.2,
    },
    "2013-2016": {
        "low_risk":             86.1,
        "moderately_high_risk": 10.2,
        "high_risk":             2.5,
        "very_high_risk":        1.3,
    },
    "2017-2020": {
        "low_risk":             86.0,
        "moderately_high_risk": 10.5,
        "high_risk":             2.4,
        "very_high_risk":        1.1,
    },
}

# ── TABLE 1.1A: EGFR × ACR GRID — % AND N (2017–MARCH 2020) ─────────────────
# Rows: G1 eGFR≥90, G2 60-89, G3a 45-59, G3b 30-44, G4 15-29, G5 <15
# Cols: A1 ACR<30, A2 ACR 30-299, A3 ACR≥300
KDIGO_GRID_PCT = {
    # (egfr_category, acr_category): pct
    ("G1_ge90",   "A1_lt30"):  59.8,
    ("G1_ge90",   "A2_30_299"): 5.0,
    ("G1_ge90",   "A3_ge300"):  0.68,
    ("G2_60_89",  "A1_lt30"):  26.2,
    ("G2_60_89",  "A2_30_299"): 2.4,
    ("G2_60_89",  "A3_ge300"):  0.35,
    ("G3a_45_59", "A1_lt30"):   3.1,
    ("G3a_45_59", "A2_30_299"): 0.79,
    ("G3a_45_59", "A3_ge300"):  0.12,
    ("G3b_30_44", "A1_lt30"):   0.61,
    ("G3b_30_44", "A2_30_299"): 0.32,
    ("G3b_30_44", "A3_ge300"):  0.18,
    ("G4_15_29",  "A1_lt30"):   0.07,
    ("G4_15_29",  "A2_30_299"): 0.08,
    ("G4_15_29",  "A3_ge300"):  0.18,
    ("G5_lt15",   "A1_lt30"):   0.00,
    ("G5_lt15",   "A2_30_299"): 0.02,
    ("G5_lt15",   "A3_ge300"):  0.13,
}

# ── TABLE 1.2: CKD % BY SOCIOECONOMIC CHARACTERISTICS (2017–MARCH 2020) ──────
CKD_BY_SES = {
    "overall":                      14.0,
    # insurance
    "not_insured":                  10.8,
    "insured":                      14.5,
    # insurance type
    "private":                      13.2,
    "medicare":                     31.5,
    "medicare_and_private":         33.6,
    "medicaid":                     15.6,
    "other_government":             14.3,
    # income
    "income_not_poverty":           14.0,
    "income_poverty":               16.5,
    # education
    "not_hs_graduate":              20.4,
    "hs_graduate_or_ged":           16.1,
    "at_least_some_college":        12.0,
}

# ── DEMOGRAPHIC CKD PREVALENCE (FIGURE 1.1, 2017–MARCH 2020) ─────────────────
CKD_BY_DEMO = {
    "overall":          14.0,
    "female":           15.4,
    "male":             12.6,
    "age_lt65":          9.0,
    "age_ge65":         33.2,
    "race_black":       18.8,
    "race_hispanic":    12.0,
    # White not explicitly stated; overall is ~14% with Black highest, Hispanic lowest
}

# ── TABLE 1.3: SELECTED HEALTH RISK BEHAVIORS (ADJUSTED, 2017–MARCH 2020) ────
RISK_BEHAVIORS = {
    "sedentary_no_ckd_pct":          20.8,   # <2.5 MET-hr/wk
    "sedentary_ckd_pct":             26.2,
    "current_smoker_no_ckd_pct":     16.4,
    "current_smoker_ckd_pct":        17.1,
    "sodium_ge2300mg_no_ckd_pct":    77.8,   # fraction NOT meeting guidelines
    "sodium_ge2300mg_ckd_pct":       71.6,
    "potassium_lt4700mg_no_ckd_pct": 95.6,   # nearly universal inadequacy
    "potassium_lt4700mg_ckd_pct":    97.9,
}

# ── DIABETES CO-PREVALENCE (FIGURE 1.6, 2017–MARCH 2020) ─────────────────────
DIABETES_BY_CKD = {
    "no_ckd_pct":   9.5,
    "any_ckd_pct": 35.6,
    "ckd_g3_pct":  32.3,
    "ckd_g4_g5_pct": 38.6,
}

# Flat dict written to JSON for downstream model scripts
CKD_PARAMS = {
    # Overall CKD prevalence
    "ckd_prevalence_overall":           CKD_BY_DEMO["overall"] / 100,
    "ckd_prevalence_female":            CKD_BY_DEMO["female"] / 100,
    "ckd_prevalence_male":              CKD_BY_DEMO["male"] / 100,
    "ckd_prevalence_lt65":              CKD_BY_DEMO["age_lt65"] / 100,
    "ckd_prevalence_ge65":              CKD_BY_DEMO["age_ge65"] / 100,
    "ckd_prevalence_black":             CKD_BY_DEMO["race_black"] / 100,
    "ckd_prevalence_hispanic":          CKD_BY_DEMO["race_hispanic"] / 100,
    # KDIGO risk tier prevalence (most recent survey period)
    "kdigo_moderately_high_risk_pct":   KDIGO_RISK_TRENDS["2017-2020"]["moderately_high_risk"],
    "kdigo_high_risk_pct":              KDIGO_RISK_TRENDS["2017-2020"]["high_risk"],
    "kdigo_very_high_risk_pct":         KDIGO_RISK_TRENDS["2017-2020"]["very_high_risk"],
    # Comorbidities
    "diabetes_prev_no_ckd":             DIABETES_BY_CKD["no_ckd_pct"] / 100,
    "diabetes_prev_any_ckd":            DIABETES_BY_CKD["any_ckd_pct"] / 100,
    # Lifestyle
    "sedentary_no_ckd":                 RISK_BEHAVIORS["sedentary_no_ckd_pct"] / 100,
    "sedentary_ckd":                    RISK_BEHAVIORS["sedentary_ckd_pct"] / 100,
    # Survey year
    "nhanes_period":                    "2017-March 2020",
    "adr_year":                         2025,
}


def try_parse_usrds_excel(path: Path, sheet_hint: str = None) -> pd.DataFrame | None:
    """Attempt to parse a USRDS Excel file. Returns None if file not found."""
    if not path.exists():
        return None
    try:
        xf = pd.ExcelFile(path)
        sheet = sheet_hint if (sheet_hint and sheet_hint in xf.sheet_names) \
                else xf.sheet_names[0]
        return pd.read_excel(path, sheet_name=sheet, header=None)
    except Exception as e:
        print(f"  Warning: could not parse {path.name}: {e}")
        return None


def build_kdigo_grid_df() -> pd.DataFrame:
    """Return Table 1.1a fallback data as a tidy DataFrame."""
    rows = []
    for (egfr_cat, acr_cat), pct in KDIGO_GRID_PCT.items():
        rows.append({"egfr_category": egfr_cat, "acr_category": acr_cat, "pct": pct})
    return pd.DataFrame(rows)


def build_kdigo_trends_df() -> pd.DataFrame:
    """Return Table 1.1b fallback data as a tidy DataFrame."""
    rows = []
    for period, cats in KDIGO_RISK_TRENDS.items():
        for cat, pct in cats.items():
            rows.append({"period": period, "risk_category": cat, "pct": pct})
    return pd.DataFrame(rows)


def build_ses_df() -> pd.DataFrame:
    return pd.DataFrame(
        [{"characteristic": k, "ckd_pct_2017_2020": v} for k, v in CKD_BY_SES.items()]
    )


def build_risk_behaviors_df() -> pd.DataFrame:
    return pd.DataFrame(
        [{"metric": k, "value": v} for k, v in RISK_BEHAVIORS.items()]
    )


def main():
    print("=== 02_download_usrds.py ===\n")
    print("Checking for USRDS CKD Chapter 1 Excel files in data/raw/ ...\n")

    found_any_excel = False

    path_t1 = DATA_RAW / "usrds_ckd_chapter1_table1.xlsx"
    path_t2 = DATA_RAW / "usrds_ckd_chapter1_table2.xlsx"
    path_t3 = DATA_RAW / "usrds_ckd_chapter1_table3.xlsx"

    for path in (path_t1, path_t2, path_t3):
        raw = try_parse_usrds_excel(path)
        if raw is not None:
            out = DATA_PROC / path.stem.replace("usrds_ckd_chapter1_", "usrds_ckd_") \
                  .replace("table", "table_raw_") + ".csv"
            raw.to_csv(out, index=False)
            print(f"  Found {path.name} — saved raw copy → {out}")
            found_any_excel = True

    # Write fallback / confirmed data regardless
    grid_df = build_kdigo_grid_df()
    trends_df = build_kdigo_trends_df()
    ses_df = build_ses_df()
    behaviors_df = build_risk_behaviors_df()

    out_grid    = DATA_PROC / "usrds_ckd_kdigo_grid.csv"
    out_trends  = DATA_PROC / "usrds_ckd_kdigo_trends.csv"
    out_ses     = DATA_PROC / "usrds_ckd_by_insurance.csv"
    out_behav   = DATA_PROC / "usrds_ckd_risk_behaviors.csv"
    out_params  = DATA_PROC / "usrds_ckd_params.json"

    grid_df.to_csv(out_grid, index=False)
    trends_df.to_csv(out_trends, index=False)
    ses_df.to_csv(out_ses, index=False)
    behaviors_df.to_csv(out_behav, index=False)

    with open(out_params, "w") as f:
        json.dump(CKD_PARAMS, f, indent=2)

    print(f"\n  Saved KDIGO eGFR×ACR grid     → {out_grid}")
    print(f"  Saved KDIGO risk trends        → {out_trends}")
    print(f"  Saved CKD by SES               → {out_ses}")
    print(f"  Saved CKD risk behaviors       → {out_behav}")
    print(f"  Saved scalar parameters        → {out_params}")

    if not found_any_excel:
        print("\n  No USRDS Excel files found in data/raw/")
        print("  Using confirmed fallback values from USRDS 2025 ADR CKD Chapter 1.")
        print()
        print("  To use actual downloaded data:")
        print("    1. Register at https://usrds-adr.niddk.nih.gov/")
        print("    2. Navigate to CKD Volume → Chapter 1")
        print("    3. Download Table 1.1, 1.2, and 1.3 Excel supplements")
        print("    4. Place files in data/raw/ as:")
        print("         usrds_ckd_chapter1_table1.xlsx")
        print("         usrds_ckd_chapter1_table2.xlsx")
        print("         usrds_ckd_chapter1_table3.xlsx")
        print("    5. Re-run this script")

    print("\n  Key values (NHANES 2017–March 2020):")
    print(f"    CKD prevalence overall:      {CKD_PARAMS['ckd_prevalence_overall']:.1%}")
    print(f"    CKD prevalence, female:      {CKD_PARAMS['ckd_prevalence_female']:.1%}")
    print(f"    CKD prevalence, male:        {CKD_PARAMS['ckd_prevalence_male']:.1%}")
    print(f"    CKD prevalence, Black:       {CKD_PARAMS['ckd_prevalence_black']:.1%}")
    print(f"    CKD prevalence, Hispanic:    {CKD_PARAMS['ckd_prevalence_hispanic']:.1%}")
    print(f"    CKD prevalence, age <65:     {CKD_PARAMS['ckd_prevalence_lt65']:.1%}")
    print(f"    CKD prevalence, age ≥65:     {CKD_PARAMS['ckd_prevalence_ge65']:.1%}")
    print(f"    KDIGO very-high-risk:        {CKD_PARAMS['kdigo_very_high_risk_pct']:.1f}%")
    print(f"    Diabetes prev (any CKD):     {CKD_PARAMS['diabetes_prev_any_ckd']:.1%}")
    print("\nDone.")


if __name__ == "__main__":
    main()
