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
├── src/
│   ├── 01_download_lifetables.py   # CDC NVSR 2021 US life tables
│   ├── 02_download_usrds.py        # USRDS ADR 2025 CKD Chapter 1 (general CKD prevalence)
│   ├── 03_download_srtr.py         # SRTR 2023 ADR waitlist & post-transplant outcomes
│   │                               # Also outputs: srtr_graft_survival.csv,
│   │                               #               srtr_patient_survival.csv
│   ├── 04_literature_params.py     # Parameters extracted from Muzaale 2014, Grams 2016, etc.
│   ├── 05_assemble_parameters.py   # Combines all sources into params.json
│   ├── 06_markov_simulation.py     # Main Monte Carlo simulation
│   └── utils.py                    # Shared helpers
├── results/                        # Output figures and tables
└── tests/
    └── test_model.py               # Validation checks
```

## Data sources

| Source | Access | What we extract |
|--------|--------|-----------------|
| CDC NVSR 2021 US Life Tables | Free download (FTP) | Age/sex qx mortality rates |
| USRDS 2025 ADR CKD Chapter 1 | Free (usrds-adr.niddk.nih.gov) | General CKD prevalence (contextual) |
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
python src/02_download_usrds.py        # requires manual USRDS download — see script
python src/03_download_srtr.py         # requires manual SRTR download — see script
                                       # outputs srtr_params.json, srtr_graft_survival.csv,
                                       #         srtr_patient_survival.csv
python src/04_literature_params.py
python src/05_assemble_parameters.py

# Step 2: Run simulation
python src/06_markov_simulation.py

# Results appear in results/
```

## Key findings

See `results/kidney_model_results.csv` for current figures (re-run the simulation
after any parameter changes). Representative base-case results at age 40:

- **Base case (age 40, overall):** ΔLE slightly negative (donor vs matched non-donor)
- **Priority offsets:** a minority fraction of the ESRD risk cost
- **Dominant cost:** Post-transplant survival gap (~3–7%/yr mortality vs healthy peers
  depending on age) — unaffected by allocation policy
- **Age pattern:** Net harm is greatest for young donors and approaches zero by age 55
- **Race:** Black donors bear substantially more absolute ESRD cost than white donors
- **Critical uncertainty:** If Mjøen 2014 all-cause mortality HR (1.30) is correct
  rather than Muzaale/US data (HR ≈ 1.0), net harm increases substantially

## Citation

If using this model, please cite the primary data sources:
- Muzaale AD et al. JAMA 2014;311:579–586
- Grams ME et al. NEJM 2016;374:411–421
- Massie AB et al. JASN 2017;28:2749–2755
- SRTR 2023 Annual Data Report
- USRDS 2023 Annual Data Report (AJKD 2024)
- CDC National Vital Statistics Reports Vol 72 No 12 (2023)
