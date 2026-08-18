"""
03_download_srtr.py
───────────────────
Extracts waitlist mortality, wait times, and post-transplant survival parameters
from the SRTR 2023 Annual Data Report supporting data for the Kidney chapter.

─── SOURCE ────────────────────────────────────────────────────────────────────
SRTR 2023 ADR supporting figures:
  File: data/raw/Kidney_Figures_Supporting_Information.xlsx
  Obtained from: https://srtr.transplant.hrsa.gov/annual_reports/Default.aspx

Transplant cohort for survival figures: 2016–2018 recipients.

Key figures used (Kidney chapter):
  Figure KI 22  — 3-year waitlist outcomes (2018–2020 listing cohort)
  Figure KI 24  — pretransplant mortality rate overall by year
  Figure KI 25  — pretransplant mortality rate by age
  Figure KI 26  — pretransplant mortality rate by race
  Figure KI 30  — pretransplant mortality rate by DSA (2023)
  Figure KI 53  — DDKT graft survival by recipient age (KM, 2016–2018 cohort)
  Figure KI 61  — LDKT graft survival by recipient age (KM, 2016–2018 cohort)
  Figure KI 70  — DD patient survival by recipient age (KM, 2016–2018 cohort)
  Figure KI 71  — DD patient survival by race (KM, 2016–2018 cohort)
  Figure KI 76  — LD patient survival by recipient age (KM, 2016–2018 cohort)

Outputs:
  data/processed/srtr_graft_survival.csv
  data/processed/srtr_patient_survival.csv
  data/processed/srtr_params.json
"""

import sys
import json
import warnings
import pandas as pd
import numpy as np
from pathlib import Path

# Allow running from repo root or src/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import DATA_RAW, DATA_PROC, median_to_annual_tx_prob

warnings.filterwarnings("ignore", category=UserWarning)

SRTR_EXCEL = DATA_RAW / "Kidney_Figures_Supporting_Information.xlsx"

# ── CONFIRMED SRTR 2023 ADR FALLBACK VALUES ───────────────────────────────────
# Parsed directly from Kidney_Figures_Supporting_Information.xlsx.
# Transplant cohort for survival figures: 2016–2018 recipients.

SRTR_FALLBACK = {

    # ── WAITLIST MORTALITY ────────────────────────────────────────────────
    # SRTR 2023 ADR Figure KI 24 (2023 data year)
    "pretx_mort_per_100py_overall_2023": 5.0,
    # DSA range (Figure KI 30, 2023): 2.01–7.30/100 PY
    "pretx_mort_per_100py_dsa_min": 2.01,
    "pretx_mort_per_100py_dsa_max": 7.30,
    # By age (Figure KI 25, 2023)
    "pretx_mort_per_100py_age1834": 2.04,
    "pretx_mort_per_100py_age3549": 3.01,
    "pretx_mort_per_100py_age5064": 5.20,
    "pretx_mort_per_100py_age65p":  7.48,
    # By race (Figure KI 26, 2023)
    "pretx_mort_per_100py_black":    4.62,
    "pretx_mort_per_100py_white":    5.71,
    "pretx_mort_per_100py_hispanic": 4.59,

    # ── WAITLIST 3-YEAR OUTCOMES (FIGURE KI 22, 2018–2020 LISTING COHORT) ─
    "wl_3yr_still_waiting":  0.299,
    "wl_3yr_ddkt":           0.308,
    "wl_3yr_ldkt":           0.135,
    "wl_3yr_died":           0.068,
    "wl_3yr_removed_other":  0.191,
    # Cause-specific removal rate back-calculated from the KI 22 3-year CIF
    # accounting for competing transplantation and death (see _solve_removal_rate).
    # Old incorrect formula 1-(1-0.191)^(1/3) = 6.8% ignored competing risks.
    # Corrected value: ~12.6%/yr reproduces 19.1% at 3 yr in the full model.
    "wl_annual_removal_competing": 0.1260,

    # ── MEDIAN WAIT TIMES ─────────────────────────────────────────────────
    # National mean waiting time at transplant, Punjala 2024 (Transplant Proc
    # 56:1740-1751), Table 3: post-KAS250 (5/2021-4/2022) 58 months = 1765 d;
    # pre-KAS250 (8/2018-7/2019) 61 months = 1857 d.
    "wl_std_median_days": 1765,
    "wl_std_median_days_prekas250": 1857,
    # Prior living donors (PLD) — Wainright 2017 AJT, UNOS abstract 2015
    "wl_pld_median_days_overall": 100,
    "wl_pld_median_days_from_activation": 23,
    # By blood type (Schold AJT 2023)
    "wl_median_days_ab":          int(4.48 * 365.25),
    "wl_median_days_overall_km":  int(4.05 * 365.25),

    # ── DDKT GRAFT SURVIVAL BY RECIPIENT AGE (FIGURE KI 53, 5-YEAR KM) ───
    # Unadjusted Kaplan-Meier, 2016–2018 transplant cohort
    "ddkt_graft_5yr_age1834": 0.822,
    "ddkt_graft_5yr_age3549": 0.835,
    "ddkt_graft_5yr_age5064": 0.768,
    "ddkt_graft_5yr_age65p":  0.661,
    # Convenience mid-range for combined 35–64 group
    "ddkt_graft_5yr_age3564": 0.802,

    # ── LDKT GRAFT SURVIVAL BY RECIPIENT AGE (FIGURE KI 61, 5-YEAR KM) ───
    "ldkt_graft_5yr_age1834": 0.901,
    "ldkt_graft_5yr_age3549": 0.917,
    "ldkt_graft_5yr_age5064": 0.893,
    "ldkt_graft_5yr_age65p":  0.802,
    "ldkt_graft_5yr_age3564": 0.905,

    # ── PATIENT SURVIVAL POST-DDKT BY AGE (FIGURE KI 70, 5-YEAR KM) ──────
    "posttx_dd_5yr_patient_surv_age1834": 0.957,
    "posttx_dd_5yr_patient_surv_age3549": 0.914,
    "posttx_dd_5yr_patient_surv_age5064": 0.820,
    "posttx_dd_5yr_patient_surv_age65p":  0.701,
    # Annual mortality derived: 1 - surv^(1/5)
    "posttx_dd_annual_mort_age1834": 0.009,
    "posttx_dd_annual_mort_age3549": 0.018,
    "posttx_dd_annual_mort_age5064": 0.039,
    "posttx_dd_annual_mort_age65p":  0.069,

    # ── PATIENT SURVIVAL POST-DDKT BY RACE (FIGURE KI 71, 5-YEAR KM) ─────
    "posttx_dd_5yr_patient_surv_black": 0.837,
    "posttx_dd_5yr_patient_surv_white": 0.825,
    "posttx_dd_annual_mort_black": 0.036,
    "posttx_dd_annual_mort_white": 0.038,

    # ── PATIENT SURVIVAL POST-LDKT BY AGE (FIGURE KI 76, 5-YEAR KM) ──────
    "posttx_ld_5yr_patient_surv_age1834": 0.979,
    "posttx_ld_5yr_patient_surv_age3549": 0.961,
    "posttx_ld_5yr_patient_surv_age5064": 0.917,
    "posttx_ld_5yr_patient_surv_age65p":  0.819,

    # Assumed national DDKT 1-year graft survival, paired with the
    # age-specific 5-year values above to derive post-year-1 annual graft
    # failure below.  Overridden by the true age-specific 1-yr KM value when
    # parsed directly from the ADR Excel (see main()).
    "ddkt_graft_1yr_assumed": 0.955,
    # Long-term median graft survival — Schold JD et al., AJT 2021;21(5):1729–1738
    # (SRTR data 1995–2017, 2014-era cohort half-life 11.7 yr for DDKT)
    "ddkt_median_graft_surv_yr_2014era": 11.7,
}

# Age-stratified post-year-1 annual graft failure, derived as
# 1 - (5yr_survival / 1yr_survival)^(1/4) from the Figure KI 53 values above.
# A single flat rate cannot represent this: the four age bands span roughly
# 3.3%-8.8%/yr once actually computed. (An earlier version hardcoded a flat
# 2.5%/yr "conservative" value that undershot every age band; see paper
# design.tex for the correction.) The flat "graft_annual_fail_postyear1" key
# is kept only as a fallback for an unmatched/unspecified age.
for _age_key in ("age1834", "age3549", "age5064", "age65p"):
    SRTR_FALLBACK[f"graft_annual_fail_postyear1_{_age_key}"] = round(
        1 - (SRTR_FALLBACK[f"ddkt_graft_5yr_{_age_key}"]
             / SRTR_FALLBACK["ddkt_graft_1yr_assumed"]) ** 0.25, 4)
SRTR_FALLBACK["graft_annual_fail_postyear1"] = round(sum(
    SRTR_FALLBACK[f"graft_annual_fail_postyear1_{_k}"]
    for _k in ("age1834", "age3549", "age5064", "age65p")) / 4, 4)


# ── HELPERS ───────────────────────────────────────────────────────────────────

def _solve_removal_rate(
    target_removed_3yr: float,
    wl_tx_annual: float,
    wl_mort_annual: float,
    n_iter: int = 80,
) -> float:
    """
    Back-calculate the cause-specific annual removal probability from the
    3-year cumulative fraction removed, accounting for competing depletion
    by transplantation and death.

    The naive formula 1-(1-CIF)^(1/3) is wrong because it treats removal
    as the only competing event.  Here we solve for r such that propagating
    a cohort through  die → transplant → remove  for 3 annual cycles
    yields cumulative_removed == target_removed_3yr.
    """
    def cum_removed(r: float) -> float:
        WL = 1.0; total = 0.0
        for _ in range(3):
            WL -= WL * wl_mort_annual
            WL -= WL * wl_tx_annual
            rem = WL * r;  WL -= rem;  total += rem
        return total

    lo, hi = 0.0, 1.0
    for _ in range(n_iter):
        mid = (lo + hi) / 2.0
        if cum_removed(mid) < target_removed_3yr:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 4)


def _find_km_header(df):
    """Return (hdr_row_idx, yr_col_idx) for a KM figure sheet."""
    for i, row in df.iterrows():
        for j, v in enumerate(row):
            if str(v) == "Years after transplant" or str(v) == "Years after listing":
                return i, j
    return None, None


def km_at_year(xl, sheet, target_year):
    """
    Extract all column values from a KM figure sheet at the time point
    closest to target_year. Returns {col_name: float_percent} dict.
    """
    df = xl.parse(sheet, header=None)
    hdr_row, yr_col = _find_km_header(df)
    if hdr_row is None:
        return {}

    cols = df.iloc[hdr_row].tolist()
    data = df.iloc[hdr_row + 1:].copy()
    data.iloc[:, yr_col] = pd.to_numeric(data.iloc[:, yr_col], errors="coerce")
    data = data.dropna(subset=[data.columns[yr_col]])

    idx = (data.iloc[:, yr_col] - target_year).abs().idxmin()
    row = data.loc[idx]

    result = {}
    for ci, col in enumerate(cols):
        label = str(col)
        if label in ("nan", "Years after transplant", "Years after listing"):
            continue
        v = pd.to_numeric(row.iloc[ci], errors="coerce")
        if not pd.isna(v):
            result[label] = round(v, 4)
    return result


def km_to_rows(xl, sheet, donor_type, timepoints=(1, 3, 5)):
    """Build a list of {age_group, timepoint, survival, donor_type} dicts."""
    rows = []
    for tp in timepoints:
        vals = km_at_year(xl, sheet, tp)
        for col, pct in vals.items():
            rows.append({
                "age_group": col,
                "timepoint": f"{tp}yr",
                "survival": round(pct / 100, 5),
                "donor_type": donor_type,
            })
    return rows


def parse_time_series_last(xl, sheet, col_name):
    """Return the last (most recent) value in a year × metric time-series sheet."""
    df = xl.parse(sheet, header=None)
    # Header row has 'Year' in column 0
    hdr_row = None
    for i, row in df.iterrows():
        if str(row.iloc[0]).strip() == "Year":
            hdr_row = i
            break
    if hdr_row is None:
        return {}
    cols = df.iloc[hdr_row].tolist()
    data = df.iloc[hdr_row + 1:].copy()
    data.iloc[:, 0] = pd.to_numeric(data.iloc[:, 0], errors="coerce")
    data = data.dropna(subset=[data.columns[0]]).sort_values(data.columns[0])
    last_row = data.iloc[-1]
    result = {}
    for ci, col in enumerate(cols):
        label = str(col)
        if label in ("nan", "Year"):
            continue
        if col_name and label != col_name:
            continue
        v = pd.to_numeric(last_row.iloc[ci], errors="coerce")
        if not pd.isna(v):
            result[label] = round(v, 4)
    return result


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("=== 03_download_srtr.py ===\n")

    params = dict(SRTR_FALLBACK)
    graft_rows = []
    patient_rows = []

    if not SRTR_EXCEL.exists():
        print(f"  {SRTR_EXCEL.name} not found in data/raw/ — using hardcoded fallback values.\n")
    else:
        print(f"  Parsing {SRTR_EXCEL.name} ...\n")
        xl = pd.ExcelFile(SRTR_EXCEL)

        # ── Pretransplant mortality ───────────────────────────────────────
        mort_overall = parse_time_series_last(xl, "KI-F24-mort-adult-waiting-all", "overall")
        if mort_overall:
            params["pretx_mort_per_100py_overall_2023"] = round(mort_overall["overall"], 2)

        mort_age = parse_time_series_last(xl, "KI-F25-mort-adult-waiting-age", None)
        age_map = {
            "18-34 years": "pretx_mort_per_100py_age1834",
            "35-49":        "pretx_mort_per_100py_age3549",
            "50-64":        "pretx_mort_per_100py_age5064",
            "65+":          "pretx_mort_per_100py_age65p",
        }
        for col, key in age_map.items():
            if col in mort_age:
                params[key] = round(mort_age[col], 2)

        mort_race = parse_time_series_last(xl, "KI-F26-mort-adult-waiting-race", None)
        race_map = {
            "Black":    "pretx_mort_per_100py_black",
            "White":    "pretx_mort_per_100py_white",
            "Hispanic": "pretx_mort_per_100py_hispanic",
        }
        for col, key in race_map.items():
            if col in mort_race:
                params[key] = round(mort_race[col], 2)

        # ── 3-year waitlist outcomes ──────────────────────────────────────
        wl3 = km_at_year(xl, "KI-F22-3yr-outcomes-adult-waiti", 3)
        outcome_map = {
            "Still on waiting list": "wl_3yr_still_waiting",
            "DD transplant":         "wl_3yr_ddkt",
            "LD transplant":         "wl_3yr_ldkt",
            "Died":                  "wl_3yr_died",
            "Removed from list":     "wl_3yr_removed_other",
        }
        for col, key in outcome_map.items():
            if col in wl3:
                params[key] = round(wl3[col] / 100, 4)
        if "wl_3yr_removed_other" in params:
            import math
            wl_mort_a = 1.0 - math.exp(
                -params.get("pretx_mort_per_100py_overall_2023", 5.0) / 100)
            wl_tx_a   = float(median_to_annual_tx_prob(
                float(params.get("wl_std_median_days", 1765))))
            params["wl_annual_removal_competing"] = _solve_removal_rate(
                params["wl_3yr_removed_other"], wl_tx_a, wl_mort_a)

        # ── DDKT graft survival (Figure KI 53) ───────────────────────────
        ddkt_rows = km_to_rows(xl, "KI-F53-tx-adult-GF-DD-5yr-age", "DDKT")
        if ddkt_rows:
            graft_rows.extend(ddkt_rows)
            age_suffix_map = {
                "18-34 years": "age1834",
                "35-49":       "age3549",
                "50-64":       "age5064",
                "65+":         "age65p",
            }
            # Collect 1-yr and 5-yr survival per age band from the same sheet
            surv_by_age = {}
            for r in ddkt_rows:
                ag = r["age_group"]
                if ag in age_suffix_map and r["timepoint"] in ("1yr", "5yr"):
                    surv_by_age.setdefault(ag, {})[r["timepoint"]] = r["survival"]

            for ag, suffix in age_suffix_map.items():
                vals = surv_by_age.get(ag, {})
                if "5yr" in vals:
                    params[f"ddkt_graft_5yr_{suffix}"] = vals["5yr"]
                # Post-year-1 annual graft failure: use the *actual*
                # age-specific 1-yr KM value from this sheet rather than the
                # assumed national constant used in the SRTR_FALLBACK default.
                if "1yr" in vals and "5yr" in vals and vals["1yr"] > 0:
                    params[f"graft_annual_fail_postyear1_{suffix}"] = round(
                        1 - (vals["5yr"] / vals["1yr"]) ** 0.25, 4)

            _fail_keys = [f"graft_annual_fail_postyear1_{s}"
                          for s in age_suffix_map.values()]
            if all(k in params for k in _fail_keys):
                params["graft_annual_fail_postyear1"] = round(
                    sum(params[k] for k in _fail_keys) / len(_fail_keys), 4)

        # ── LDKT graft survival (Figure KI 61) ───────────────────────────
        ldkt_rows = km_to_rows(xl, "KI-F61-tx-adult-GF-LD-5yr-age", "LDKT")
        if ldkt_rows:
            graft_rows.extend(ldkt_rows)
            for r in ldkt_rows:
                if r["timepoint"] == "5yr":
                    ag = r["age_group"]
                    surv = r["survival"]
                    key_map = {
                        "18-34 years": "ldkt_graft_5yr_age1834",
                        "35-49":       "ldkt_graft_5yr_age3549",
                        "50-64":       "ldkt_graft_5yr_age5064",
                        "65+":         "ldkt_graft_5yr_age65p",
                    }
                    if ag in key_map:
                        params[key_map[ag]] = surv

        # Convenience combined 35–64 graft survival
        for prefix in ("ddkt", "ldkt"):
            s3549 = params.get(f"{prefix}_graft_5yr_age3549")
            s5064 = params.get(f"{prefix}_graft_5yr_age5064")
            if s3549 and s5064:
                params[f"{prefix}_graft_5yr_age3564"] = round((s3549 + s5064) / 2, 3)

        # ── DD patient survival by age (Figure KI 70) ────────────────────
        dd_pt_rows = km_to_rows(xl, "KI-F70-tx-adult-pat-surv-DD-5y-", "DDKT_patient")
        if dd_pt_rows:
            patient_rows.extend(dd_pt_rows)
            ps_key_map = {
                "18-34 years": ("posttx_dd_5yr_patient_surv_age1834", "posttx_dd_annual_mort_age1834"),
                "35-49":       ("posttx_dd_5yr_patient_surv_age3549", "posttx_dd_annual_mort_age3549"),
                "50-64":       ("posttx_dd_5yr_patient_surv_age5064", "posttx_dd_annual_mort_age5064"),
                "65+":         ("posttx_dd_5yr_patient_surv_age65p",  "posttx_dd_annual_mort_age65p"),
            }
            for r in dd_pt_rows:
                if r["timepoint"] == "5yr" and r["age_group"] in ps_key_map:
                    surv = r["survival"]
                    surv_key, mort_key = ps_key_map[r["age_group"]]
                    params[surv_key] = surv
                    params[mort_key] = round(1 - surv ** (1 / 5), 4)

        # ── DD patient survival by race (Figure KI 71) ───────────────────
        dd_race_rows = km_to_rows(xl, "KI-F71-tx-adult-pat-surv-DD-5y-", "DDKT_patient_race")
        if dd_race_rows:
            patient_rows.extend(dd_race_rows)
            for r in dd_race_rows:
                if r["timepoint"] == "5yr":
                    ag = r["age_group"]
                    surv = r["survival"]
                    if ag == "Black":
                        params["posttx_dd_5yr_patient_surv_black"] = surv
                        params["posttx_dd_annual_mort_black"] = round(1 - surv ** (1 / 5), 4)
                    elif ag == "White":
                        params["posttx_dd_5yr_patient_surv_white"] = surv
                        params["posttx_dd_annual_mort_white"] = round(1 - surv ** (1 / 5), 4)

        # ── LD patient survival by age (Figure KI 76) ────────────────────
        ld_pt_rows = km_to_rows(xl, "KI-F76-tx-adult-pat-surv-LD-5y-", "LDKT_patient")
        if ld_pt_rows:
            patient_rows.extend(ld_pt_rows)
            ld_ps_key_map = {
                "18-34 years": "posttx_ld_5yr_patient_surv_age1834",
                "35-49":       "posttx_ld_5yr_patient_surv_age3549",
                "50-64":       "posttx_ld_5yr_patient_surv_age5064",
                "65+":         "posttx_ld_5yr_patient_surv_age65p",
            }
            for r in ld_pt_rows:
                if r["timepoint"] == "5yr" and r["age_group"] in ld_ps_key_map:
                    params[ld_ps_key_map[r["age_group"]]] = r["survival"]

        print(f"  Parsed {len(graft_rows)} graft survival rows")
        print(f"  Parsed {len(patient_rows)} patient survival rows")

    # ── SAVE GRAFT SURVIVAL CSV ───────────────────────────────────────────
    if graft_rows:
        gdf = pd.DataFrame(graft_rows)
        out = DATA_PROC / "srtr_graft_survival.csv"
        gdf.to_csv(out, index=False)
        print(f"  Saved → {out}")

    # ── SAVE PATIENT SURVIVAL CSV ─────────────────────────────────────────
    if patient_rows:
        pdf = pd.DataFrame(patient_rows)
        out_pt = DATA_PROC / "srtr_patient_survival.csv"
        pdf.to_csv(out_pt, index=False)
        print(f"  Saved → {out_pt}")

    # ── SAVE SCALAR PARAMS JSON ───────────────────────────────────────────
    params_out = DATA_PROC / "srtr_params.json"
    with open(params_out, "w") as f:
        json.dump(params, f, indent=2)
    print(f"\n  Saved scalar parameters → {params_out}")

    # ── KEY VALUES SUMMARY ────────────────────────────────────────────────
    p = params
    print("\n  Key values written to srtr_params.json:")
    print(f"    Pretx mortality (2023):          {p['pretx_mort_per_100py_overall_2023']:.1f}/100 PY")
    print(f"    Std wait median (post-KAS250):   {p['wl_std_median_days']} days")
    print(f"    PLD wait median:                 {p['wl_pld_median_days_overall']} days")
    print(f"    DDKT 5-yr graft surv (18–34):    {p['ddkt_graft_5yr_age1834']:.1%}")
    print(f"    DDKT 5-yr graft surv (35–49):    {p['ddkt_graft_5yr_age3549']:.1%}")
    print(f"    DDKT 5-yr graft surv (50–64):    {p['ddkt_graft_5yr_age5064']:.1%}")
    print(f"    DDKT 5-yr graft surv (65+):      {p['ddkt_graft_5yr_age65p']:.1%}")
    print(f"    DD 5-yr patient surv (50–64):    {p['posttx_dd_5yr_patient_surv_age5064']:.1%}")
    print(f"    DD 5-yr patient surv (65+):      {p['posttx_dd_5yr_patient_surv_age65p']:.1%}")
    print(f"    Graft failure post-yr1 (18–34):  {p['graft_annual_fail_postyear1_age1834']:.1%}/yr")
    print(f"    Graft failure post-yr1 (35–49):  {p['graft_annual_fail_postyear1_age3549']:.1%}/yr")
    print(f"    Graft failure post-yr1 (50–64):  {p['graft_annual_fail_postyear1_age5064']:.1%}/yr")
    print(f"    Graft failure post-yr1 (65+):    {p['graft_annual_fail_postyear1_age65p']:.1%}/yr")
    print("\nDone.")


if __name__ == "__main__":
    main()
