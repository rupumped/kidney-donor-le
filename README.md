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
| Base case (age 40, overall) | 38.05 | 38.15 | −38.3 |
| Age 25 | 51.64 | 51.77 | −47.4 |
| Age 35 | 42.52 | 42.63 | −41.8 |
| Age 45 | 33.67 | 33.75 | −31.0 |
| Age 55 | 25.27 | 25.31 | −17.3 |
| Race: White | 38.08 | 38.15 | −25.1 |
| Race: Black | 37.87 | 38.09 | −73.4 |
| Sex: Female | 40.45 | 40.48 | −29.9 |
| Sex: Male | 35.97 | 36.09 | −43.6 |
| Sensitivity: Donor mort HR=1.30 (Mjøen) | — | — | −1067.9 |
| Sensitivity: Post-Tx LDKT quality | — | — | −22.7 |

Key takeaways:
- **Priority offsets** only ~0.6 days of the ESRD risk cost (base case −38.3 vs no-priority −38.9 days)
- **Dominant cost:** Post-transplant survival gap (~3–7%/yr mortality vs healthy peers depending on age) — unaffected by allocation policy
- **Age pattern:** Net harm is greatest for young donors and approaches zero by age 55
- **Race:** Black donors bear substantially more absolute ESRD cost than white donors (−73 vs −25 days)
- **Critical uncertainty:** If Mjøen 2014 all-cause mortality HR (1.30) is correct rather than Muzaale/US data (HR ≈ 1.0), net harm increases to ~3 years
- **LDKT vs DDKT:** If post-Tx survival at LDKT quality, harm falls to −23 days (vs −38 base)
- **PSA:** P(donation beneficial) ≈ 48%; mean ΔLE +29.5 days [−1086, +1068] reflecting this uncertainty

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
