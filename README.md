#  Post-Transplant Survival, Not Waitlist Time, Drives the Life-Expectancy Cost of Living Kidney Donation

Quantitative model assessing whether living donor priority policies offset the
increased ESRD risk associated with kidney donation

## Research question

Does the US kidney allocation priority afforded to prior living donors
(~100-day median wait vs ~985 days for standard candidates) offset the
elevated ESRD risk from donation (~3.7–4.5× matched controls), making
donation life-expectancy neutral or net positive?

**Short answer from the model:** No. Priority offsets only a fraction of the ESRD cost.
The dominant cost is the permanent post-transplant survival gap, which no allocation
policy can address.

## Repository structure

```
kidney-donor-le/
├── README.md
├── requirements.txt
├── data/
│   ├── raw/               # Downloaded source files (not committed — see below)
│   └── processed/         # Cleaned parameter tables produced by src/01–05 scripts
│       ├── lifetable_combined_2021.csv
│       ├── lifetable_female_2021.csv
│       ├── lifetable_male_2021.csv
│       ├── literature_params.json
│       ├── params.json
│       ├── srtr_graft_survival.csv
│       ├── srtr_params.json
│       ├── srtr_patient_survival.csv
│       ├── usrds_ckd_by_insurance.csv
│       ├── usrds_ckd_kdigo_grid.csv
│       ├── usrds_ckd_kdigo_trends.csv
│       ├── usrds_ckd_params.json
│       ├── usrds_ckd_risk_behaviors.csv
│       ├── usrds_esrd7_params.json
│       └── usrds_params.json
├── paper/                           # LaTeX files
├── src/
│   ├── 01_download_lifetables.py    # CDC NVSR 2021 US life tables
│   ├── 02_download_usrds.py         # USRDS ADR 2025 CKD Chapter 1 (general CKD)
│   ├── 02b_download_usrds_esrd7.py  # USRDS ADR 2025 ESRD Chapter 7 (transplant);
│   │                                #   derives calibrated waitlist listing probability
│   ├── 03_download_srtr.py          # SRTR 2023 ADR waitlist & post-transplant outcomes
│   │                                #   outputs: srtr_params.json,
│   │                                #            srtr_graft_survival.csv,
│   │                                #            srtr_patient_survival.csv
│   ├── 04_literature_params.py      # Parameters from Muzaale 2014, Grams 2016, etc.
│   ├── 05_assemble_parameters.py    # Combines all sources into params.json
│   ├── 06_markov_simulation.py      # Main Monte Carlo simulation (PSA + OWSA)
│   ├── 07_cohort_markov.py          # Analytic (deterministic) cohort Markov model
│   ├── 08_calibration.py            # Face-validity checks vs SRTR KI 22 / KI 70
│   ├── 09_voucher_cohort_markov.py  # Voucher-holder cohort: LE benefit of priority
│   │                                #   access for a healthy non-donor family member
│   ├── 10_esrd_conditional_cohort_markov.py  # ESRD-conditional cohort: LE benefit
│   │                                          #   of priority vs standard waitlist,
│   │                                          #   conditioned on ESRD onset
│   └── utils.py                     # Shared helpers
├── results/                         # Output figures and tables
│   ├── kidney_model_results.csv          # Summary table of all analyses
│   ├── calibration_report.txt            # Face-validity calibration summary
│   ├── kidney_model_A_distributions.png  # Input parameter distributions
│   ├── kidney_model_B_psa.png            # Probabilistic sensitivity analysis
│   ├── kidney_model_C_tornado.png        # One-way sensitivity (tornado) diagram
│   ├── kidney_model_D_age_subgroup.png   # ΔLE by age at donation
│   ├── kidney_model_E_race_subgroup.png  # ΔLE by race
│   ├── kidney_model_age_race_matrix.png  # ΔLE heat map: age × race
│   ├── kidney_model_age_race_sex_matrix.png
│   ├── kidney_model_age_sex_matrix.png   # ΔLE heat map: age × sex
│   ├── kidney_model_sex_race_matrix.png  # ΔLE heat map: sex × race
│   ├── kidney_model_race.png             # Race subgroup summary
│   ├── kidney_model_sex.png              # Sex subgroup summary
│   ├── cohort_markov_state_occupancy.png      # State occupancy over time (cohort Markov)
│   ├── cohort_markov_survival.png             # Survival curves (cohort Markov)
│   ├── voucher_cohort_state_occupancy.png     # State occupancy — voucher vs control
│   ├── voucher_cohort_survival.png            # Survival curves — voucher vs control
│   ├── voucher_cohort_owsa_tornado.png        # OWSA tornado — voucher LE benefit
│   ├── voucher_cohort_age_sweep.png           # Voucher LE benefit by age at designation
│   ├── esrd_conditional_state_occupancy.png   # State occupancy — priority vs standard | ESRD
│   ├── esrd_conditional_survival.png          # Survival from ESRD onset
│   ├── esrd_conditional_owsa_tornado.png      # OWSA tornado — priority waitlist LE benefit | ESRD
│   └── esrd_conditional_age_sweep.png         # Priority waitlist benefit by age at ESRD onset
└── tests/
    └── test_model.py               # Validation checks
```

## Data sources

| Source | Access | What we extract |
|--------|--------|-----------------|
| CDC NVSR 2021 US Life Tables | Free download (FTP) | Age/sex qx mortality rates |
| USRDS 2025 ADR CKD Chapter 1 | Free (usrds-adr.niddk.nih.gov) | General CKD prevalence (contextual) |
| USRDS 2025 ADR ESRD Chapter 7 | Free (usrds-adr.niddk.nih.gov) | Waitlisting rates, preemptive listing, transplant outcomes |
| SRTR 2023 ADR Kidney chapter | Free (srtr.transplant.hrsa.gov) | Waitlist mortality, wait times, graft & patient survival |
| Muzaale 2014 (JAMA 311:579) | PubMed/PMC free | Donor ESRD risk by race |
| Grams 2016 (NEJM 374:411) | Published table | Pre-donation ESRD baseline by race/sex |
| Massie 2017 (JASN 28:2749) | PubMed | Post-donation ESRD hazard ratios |
| Wainright 2017 (AJT 17:1103) | Published | PLD median wait post-KAS |
| USRDS ADR (dialysis survival) | Hardcoded from published tables | HD 1-yr and subsequent-year mortality |

Raw data files are **not committed** to the repo. Run `src/01_download_lifetables.py`
through `src/04_literature_params.py` in order to populate `data/raw/` and
`data/processed/`. SRTR data requires free registration; see notes in each script.
USRDS CKD Chapter 1 data requires registration; dialysis survival parameters are
hardcoded directly from published USRDS ADR tables.

## Quickstart

```bash
pip install -r requirements.txt

# Step 1: Download and process all parameters
python src/01_download_lifetables.py
python src/02_download_usrds.py         # requires manual USRDS download — see script
python src/02b_download_usrds_esrd7.py  # encodes USRDS ESRD Ch.7; no download needed
python src/03_download_srtr.py          # requires manual SRTR download — see script
                                        # outputs srtr_params.json, srtr_graft_survival.csv,
                                        #         srtr_patient_survival.csv
python src/04_literature_params.py
python src/05_assemble_parameters.py

# Step 2: Run simulation
python src/06_markov_simulation.py      # Monte Carlo (PSA + OWSA + subgroups)
python src/07_cohort_markov.py          # Analytic cohort Markov (optional, --age N)
python src/08_calibration.py            # Face-validity checks vs SRTR KI 22 / KI 70
python src/09_voucher_cohort_markov.py  # Voucher-holder LE benefit (optional, --age N)
python src/10_esrd_conditional_cohort_markov.py  # Priority benefit | ESRD (optional, --age N)

# Results appear in results/
```

## Key findings

See [results/kidney_model_results.csv](results/kidney_model_results.csv) for current figures (re-run the simulation
after any parameter changes). Base-case results:

| Analysis | LE Donor (yr) | LE Non-donor (yr) | ΔLE (days) |
|----------|--------------|-------------------|------------|
| Base case (age 40, overall) | 35.63 | 38.14 | −914.3 |
| Age 25 | 48.65 | 51.75 | −1131.6 |
| Age 35 | 39.90 | 42.62 | −991.5 |
| Age 45 | 31.51 | 33.74 | −814.0 |
| Age 55 | 23.69 | 25.31 | −591.2 |
| Race: White | 35.68 | 38.13 | −892.5 |
| Race: Black | 35.39 | 38.01 | −950.7 |
| Sex: Female | 38.08 | 40.46 | −888.3 |
| Sex: Male | 33.55 | 36.07 | −915.9 |
| Sensitivity: Donor mort HR=1.30 flat (Mjøen) | — | — | −1084.5 |
| Sensitivity: Post-Tx LDKT quality | — | — | −889.8 |

Key takeaways:
- **Priority offsets** only ~0.2 days of the ESRD risk cost (base case −909.3 vs no-priority −909.5 days)
- **Dominant cost:** Time-varying donor all-cause mortality HR (flat at 1.0 through year 10, ramping to Mjøen's 1.30 by year 15) — unaffected by allocation policy
- **Age pattern:** Net harm is greatest for young donors and falls substantially by age 55, since a younger donor spends more of their remaining life under the elevated HR
- **Race:** Black donors bear somewhat more cost than white donors (−951 vs −893 days at age 40), though this ESRD-driven gap is now a minor contributor next to the mortality-HR effect
- **Critical uncertainty:** The donor mortality HR assumption dominates every other source of uncertainty — holding it flat at the literature's individual point estimates ranges from a multi-year *benefit* (HR=0.60, O'Keeffe pooled estimate) to a cost exceeding 4.6 years (HR=1.52, Mjøen upper 95% CI)
- **LDKT vs DDKT:** If post-Tx survival at LDKT quality, harm falls to −889.8 days (vs −909.3 base)
- **PSA:** P(donation beneficial) = 0%; mean ΔLE −910.8 days [−936.1, −889.3], reflecting narrow residual uncertainty once the mortality HR is held fixed at its base-case time-varying profile

## Testing

```bash
python3 -m pytest tests/ -v
```

Tests require `data/processed/` to be populated (run scripts 01–05 first).

## Citation
This model cites the primary data sources:
- Muzaale AD et al. JAMA 2014;311:579–586
- Grams ME et al. NEJM 2016;374:411–421
- Massie AB et al. JASN 2017;28:2749–2755
- SRTR 2023 Annual Data Report
- USRDS 2025 Annual Data Report (ESRD Volume, Chapters 1 & 7)
- CDC National Vital Statistics Reports Vol 72 No 12 (2023)
