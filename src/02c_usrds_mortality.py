"""
02c_usrds_mortality.py
──────────────────────
Parses USRDS 2023 ADR, ESRD Volume, Chapter 6 (Mortality) figure/table data
exports and derives age-stratified dialysis mortality for the Markov model.

─── SOURCE ────────────────────────────────────────────────────────────────────
USRDS 2023 ADR, ESRD Volume, Chapter 6: Mortality
  https://usrds-adr.niddk.nih.gov/2023/end-stage-renal-disease/6-mortality

Manual download (CSV export from each figure/table page) into
data/raw/usrds_ch6/:
  esrd-Figure-6.5b.csv               — mortality per 1,000 PY by age/race and
                                       modality (All ESRD/Dialysis/Transplant),
                                       2019–2021, adjusted and unadjusted
  esrd-Figure-6.4.csv                — mortality per 1,000 PY by year since
                                       ESRD onset, incident cohorts
  esrd-table-6.1-esrd.csv            — expected remaining years of life on
                                       dialysis / with transplant, by age & sex
  esrd-table-6.1-us-population.csv   — same, general US population
  esrd-Figure-6.7.csv                — 5-yr survival of incident ESRD patients
                                       by modality (not age-stratified)

─── DERIVATIONS ───────────────────────────────────────────────────────────────
Age-band dialysis mortality (Figure 6.5b, Dialysis, Adjusted, 2019):
  annual probability = 1 − exp(−rate / 1000). 2019 is used rather than
  2020–2021 because COVID-19 raised dialysis mortality ~20% in those years;
  2021 values are kept for sensitivity.

First-year multiplier (Figure 6.4, Hemodialysis, adjusted, censored at
transplant, 2019 incident cohort): rate in year 1 ÷ rate in year 2. Year 2 of
the same incident cohort is one year older, so the ratio isolates the
early-dialysis excess without the age confounding of a comparison against
the prevalent population (whose age is already handled by the age bands).

Output:
  data/processed/usrds_mortality_params.json
"""

import sys
import csv
import json
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import DATA_RAW, DATA_PROC

RAW = DATA_RAW / "usrds_ch6"

AGE_BANDS = {"18-44": "age1844", "45-64": "age4564", "65-74": "age6574", "75+": "age75p"}


def _read_rows(name: str) -> list:
    with open(RAW / name, newline="", encoding="utf-8-sig") as f:
        return list(csv.reader(f))


def rate_to_prob(rate_per_1000py: float) -> float:
    """Constant-hazard conversion: per-1,000-PY rate → annual probability."""
    return 1.0 - math.exp(-rate_per_1000py / 1000.0)


def parse_fig_6_5b() -> dict:
    """Age-band mortality per 1,000 PY: {(tab, adjustment, band, year): rate}."""
    rows = _read_rows("esrd-Figure-6.5b.csv")
    header, out = rows[0], {}
    years = header[3:]
    for r in rows[1:]:
        tab, adj, cat = r[0], r[1], r[2]
        if cat not in AGE_BANDS:
            continue
        for yr, val in zip(years, r[3:]):
            out[(tab, adj, AGE_BANDS[cat], yr)] = float(val)
    return out


def parse_fig_6_4() -> dict:
    """Mortality per 1,000 PY by year since ESRD onset, keyed by
    (censor_at_tx, modality, adjusted, cohort_year) → {year_since_onset: rate}."""
    rows = _read_rows("esrd-Figure-6.4.csv")
    censor, modality, adjusted, cohort = (r[1:] for r in rows[:4])
    out = {}
    for j in range(len(cohort)):
        key = (censor[j], modality[j], adjusted[j], cohort[j])
        out[key] = {int(r[0]): float(r[j + 1]) for r in rows[4:] if r[j + 1]}
    return out


def parse_table_6_1() -> dict:
    """Expected remaining years of life, 2019, by 5-yr age group and sex, for
    dialysis, transplant, and the general US population."""
    esrd = _read_rows("esrd-table-6.1-esrd.csv")
    pop  = _read_rows("esrd-table-6.1-us-population.csv")
    out = {}
    # ESRD file: columns Age | Dial F | Dial M | Tx F | Tx M (2019 block first)
    for r in esrd[3:]:
        age = r[0]
        for col, (mod, sex) in enumerate([("dialysis", "female"), ("dialysis", "male"),
                                          ("transplant", "female"), ("transplant", "male")], 1):
            val = r[col].rstrip("a")   # 'a' footnote: estimate from small cell
            if val:
                out[f"le_{mod}_{sex}_{age}"] = float(val)
    # Population file: columns Age | F | M (2019 block first)
    for r in pop[2:]:
        out[f"le_uspop_female_{r[0]}"] = float(r[1])
        out[f"le_uspop_male_{r[0]}"]   = float(r[2])
    return out


def parse_fig_6_7() -> dict:
    """Incident-cohort survival (%) at 12/36/60 months, 2017 cohort, by modality."""
    rows = _read_rows("esrd-Figure-6.7.csv")
    modality, cohort = rows[0][1:], rows[1][1:]
    by_month = {int(r[0]): r[1:] for r in rows[2:] if r[0]}
    out = {}
    for j, (mod, yr) in enumerate(zip(modality, cohort)):
        if yr != "2017":
            continue
        tag = mod.lower().replace(" ", "_").replace("-", "_")
        for m in (12, 36, 60):
            if m in by_month and by_month[m][j]:
                out[f"surv_{m}mo_{tag}_2017"] = float(by_month[m][j]) / 100.0
    return out


def main():
    print("=== 02c_usrds_mortality.py ===\n")
    f65 = parse_fig_6_5b()
    f64 = parse_fig_6_4()

    params = {}
    # Age-band dialysis mortality: base case 2019 adjusted; 2021 for sensitivity
    for band in AGE_BANDS.values():
        for yr in ("2019", "2021"):
            rate = f65[("Dialysis", "Adjusted", band, yr)]
            sfx = "" if yr == "2019" else "_2021"
            params[f"dialysis_rate_per1000py_{band}{sfx}"] = rate
            params[f"dialysis_mort_{band}{sfx}"] = round(rate_to_prob(rate), 5)
        params[f"transplant_rate_per1000py_{band}"] = f65[("Transplant", "Adjusted", band, "2019")]

    # First-year excess: 2019 incident HD cohort, adjusted, censored at transplant
    hd19 = f64[("Yes", "Hemodialysis", "Yes", "2019")]
    params["dialysis_rate_yr1_hd_2019"] = hd19[1]
    params["dialysis_rate_yr2_hd_2019"] = hd19[2]
    params["dialysis_yr1_mult"] = round(hd19[1] / hd19[2], 4)

    params.update(parse_table_6_1())
    params.update(parse_fig_6_7())
    params["_source"] = ("USRDS 2023 ADR, ESRD Vol. Ch.6: Figures 6.4, 6.5b, 6.7; "
                         "Table 6.1. CSV exports in data/raw/usrds_ch6/.")

    out = DATA_PROC / "usrds_mortality_params.json"
    with open(out, "w") as f:
        json.dump(params, f, indent=2)
    print(f"  Saved → {out}\n")

    print("  Dialysis mortality, 2019 adjusted (Fig 6.5b):")
    for label, band in AGE_BANDS.items():
        print(f"    {label:>6}: {params[f'dialysis_rate_per1000py_{band}']:6.1f}/1000 PY"
              f"  →  {params[f'dialysis_mort_{band}']:.1%}/yr")
    print(f"  First-year multiplier (Fig 6.4, 2019 HD): "
          f"{hd19[1]:.1f} / {hd19[2]:.1f} = {params['dialysis_yr1_mult']:.3f}")
    print("\nDone.")


if __name__ == "__main__":
    main()
