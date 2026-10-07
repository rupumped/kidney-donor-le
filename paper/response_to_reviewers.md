# Response to Reviewers

**Manuscript:** Life Expectancy Changes Due to Kidney Donation Across Race, Sex, and Age

Thank you for your careful reading. Your comment identified the most important weakness in the model's treatment of outcomes after ESRD, and addressing it changed the analysis substantially. Changes in the revised manuscript appear in red. Below I reproduce the comment, summarize the revisions, and then respond to each part in turn.

---

## Comment 1

> I think my biggest concern is that the inputs are from the wrong population. Dialysis mortality, listing rates, waitlist mortality, and post-transplant survival come from all transplant candidates, not prior donors. Prior donors with ESKD are closely followed and more likely to be listed preemptively or find a living donor, they are also much healthier (generally). The dialysis mortality figures also aren't age-stratified and are driven by older patients, so they likely overstate harm for a donor-age cohort. Muzaale's paper is a good starting point for donor-specific inputs.

### Summary of response

I agree on every point and have revised the model accordingly:

1. **Donor-specific inputs after ESRD.** The donor arm now uses inputs from Muzaale et al. (2016), *Transplantation* 100:1306, which followed 99 US living donors who developed ESRD. Two inputs are taken directly from the study, and two are calibrated to its observed listing and mortality.
2. **Age-stratified dialysis mortality.** The flat 22% (first year) and 17% (later years) rates have been replaced with USRDS age-specific rates, which are much lower at donor ages (6.8%/year at ages 18–44).
3. **Age-stratified waitlist mortality.** The all-ages SRTR rate has been replaced with SRTR's age-specific rates.

Making these changes exposed related problems, which I have also corrected (point 5 below). The headline result is largely unchanged: about 890 days (2.4 years) versus 912 days (2.5 years). It remains dominated by the donor all-cause mortality hazard ratio. However, the cost that runs through ESRD fell by about a fifth, the racial disparity narrowed considerably, and the conditional value of priority access changed.

---

### 1. Inputs drawn from all candidates rather than prior donors

The donor arm now uses donor-specific inputs at every step after ESRD onset (Methods, *Outcomes after ESRD in prior donors*).

| Input | Original (all candidates) | Revised (prior donors) | Basis |
|---|---|---|---|
| Listed before dialysis at ESRD onset | 9.4% (USRDS, ages 18–44) | **21.2%** | 21 of 99 donors (Muzaale 2016) |
| Annual listing probability on dialysis | 0.15 | **0.31** | Calibrated to 43.4% listed within 12 months |
| Dialysis mortality | General population | **× 0.78** | Calibrated to observed 1/3/5/10-year mortality after ESRD |
| Post-transplant mortality | SRTR DDKT recipients | **× 0.7** | Adjusted HR, prior donors vs matched non-donor recipients (95% CI 0.2–2.4) |

The listing probability and dialysis multiplier were fitted jointly to the study's donors, starting from ESRD onset at their median age of 50. For each candidate multiplier, I solved for the listing probability that reproduces 43.4% listed within 12 months. I then kept the multiplier that best matched observed mortality at 1, 3, 5 and 10 years (9.1, 20.2, 25.5 and 32.4%). The fitted model predicts 7.6, 17.2, 23.8 and 37.2%, all within the published 95% confidence intervals. The fitted multiplier (0.78) is consistent with the study's overall mortality HR of 0.7 (95% CI 0.4–1.0).

Because the source cohort is small, all four inputs are now sampled in the probabilistic sensitivity analysis, with uncertainty taken from the published confidence intervals. The one-way sensitivity analysis includes a scenario that replaces all four with general-population inputs.

These inputs apply only to the donor arm. The voucher and ESRD-conditional analyses are meant to isolate the value of waitlist priority, and voucher holders are not donors, so neither arm of those analyses uses them.

### 2. Close follow-up, preemptive listing, and living-donor transplantation

Close follow-up and preemptive listing are captured by the revised preemptive fraction (21.2%) and listing probability (0.31/year).

On finding a living donor, Muzaale et al. report the opposite of what one might expect. Prior donors received living-donor kidneys *less* often than matched non-donors (13% vs 39%), and standard-criteria deceased-donor kidneys far more often (87%). This is consistent with priority allocation making a deceased-donor kidney available within months (median wait 2.8 vs 21.5 months). The model represents this through the priority waitlist and the donor post-transplant mortality multiplier rather than through a separate living-donor pathway. The one-way sensitivity analysis also substitutes LDKT-quality post-transplant survival, which changes the result by 8 days.

### 3. Prior donors are healthier

This is captured by the calibrated dialysis-mortality multiplier (0.78) and the post-transplant mortality multiplier (0.7).

As an additional check, I ran the general-population pathway against Muzaale et al.'s matched non-donors with ESRD, a comparison that was not used in fitting. The model reproduces their 3- and 5-year mortality (25.9 vs 26.1%; 37.4 vs 34.9%). It under-predicts first-year mortality (10.6 vs 15.6%) and over-predicts 10-year mortality (57.5 vs 43.0%). That pattern fits front-loaded deaths in a heterogeneous ESRD population, which a cohort model that treats all ESRD patients alike cannot capture. Differences in transplant access explain little of the gap: giving the general pathway the matched non-donors' 21.5-month median wait reduces 10-year mortality only from 57% to 53%. I note this in the Methods and Discussion. It has little effect on the main result, because so few non-donors develop ESRD, but it may overstate the conditional value of priority access.

### 4. Dialysis mortality is not age-stratified and is driven by older patients

This was correct, and the effect was large.

- **Dialysis.** Mortality now comes from the USRDS 2023 Annual Data Report, Figure 6.5b: 2019 adjusted rates for prevalent dialysis patients. Per year, these are 6.8% at ages 18–44, 11.4% at 45–64, 17.8% at 65–74 and 25.1% at 75+. I used 2019 because COVID-19 inflated 2020–2021 mortality by roughly 20%. First-year mortality is 1.06 times the later-year rate, from Figure 6.4. As a check, these rates reproduce USRDS remaining life expectancy on dialysis (Table 6.1) within 15% for every 5-year age group from 40–44 to 80–84.
- **Waitlist.** Mortality now comes from SRTR 2023 Figure KI 25, also 2019: 1.98, 2.95, 4.97 and 7.43 deaths per 100 patient-years at ages 18–34, 35–49, 50–64 and 65+, in place of a single all-ages rate of 5.0. The race subgroups multiply these by the race-specific ratio from Figure KI 26.
- **Smoothing.** All age-banded inputs are now interpolated smoothly between band midpoints. Step functions at band edges had created artifacts. For example, post-transplant mortality jumping from 3.9% to 6.9% at exactly age 65 made the conditional value of priority access rise between ages 65 and 73.

### 5. Related problems found while addressing this comment

- **ESRD was protective at very old ages.** The original model's fixed post-ESRD mortality rates fell below general-population mortality at advanced ages; post-transplant mortality of 6.9%/year is exceeded by the life table from age 83. In the original model, 877 per million people in ESRD states were still alive at age 110. Mortality in every post-ESRD state is now floored at the life-table rate for the same age. This was worth about 17 days in the original base case and would have grown once the lower donor-specific rates were introduced.
- **The non-donor ESRD risk gradient was borrowed.** The non-donor baseline was scaled with age at 1.40 per decade, a figure taken from non-Black donors in Massie et al. (2017). It now follows Grams et al. (2016) projections at ages 20, 40 and 60 by race and sex. Grams shows Black non-donor risk flattening after age 40 (about 1.13-fold per decade, against 1.45 for White). Correcting this removed a small, spurious reversal of the racial disparity at age 55. The voucher analysis had also given every age the age-40 risk; it now scales risk to the entry age.

---

### Effect on results

| Result | Original | Revised |
|---|---|---|
| Base-case ΔLE, age 40 (analytic) | −912 days (2.5 yr) | **−891 days (2.4 yr)** |
| ESRD-attributable cost (HR = 1.0 flat) | −68 days | **−55 days** |
| Same, with general-population inputs in the donor arm | — | −92 days |
| Probabilistic analysis, mean (95% CrI) | −911 (−936 to −889) | −896 (−941 to −862) |
| Age 25 / age 55 | −1,132 / −591 days | −1,092 / −592 days |
| Black vs White donors, age 25 | 28% larger | **16% larger** |
| Black vs White donors, age 40 | 7% larger | **2% larger** |
| Black vs White donors, age 55 | Slightly reversed | **Equal (583 vs 587 days)** |
| Male vs female donors, age 40 | 3% | 2% |
| Losing PLD priority (change in ΔLE) | −0.2 days | −6.7 days |
| Priority value given ESRD at age 40 / 60 | +539 / +109 days | **+434 / +154 days** |
| Voucher value at age 25 / 40 | +1.2 / <1 day | +0.6 / +0.5 days |

In the original model, general-population inputs overstated the cost of donation that runs through ESRD by about 70% (92 vs 55 days). The headline changes little because most of the cost comes from the donor all-cause mortality hazard ratio, which is unaffected by this comment. The race and sex findings changed more: once donors with ESRD fare as well as the evidence suggests, race and sex are small modifiers except at young donation ages. The Abstract, Highlights, Results and Discussion have been rewritten to reflect this.

### Changes to the manuscript

- **Introduction:** new sentence on why general ESRD outcomes do not transfer to prior donors.
- **Methods:**
  - *Dialysis survival and waitlist listing* rewritten;
  - new subsection *Outcomes after ESRD in prior donors*;
  - waitlist mortality now age-stratified;
  - new subsection *Age dependence and background floor for post-ESRD states*;
  - *Subgroup analyses* updated for the Grams age gradient and the race-specific waitlist adjustment;
  - donor inputs added to the probabilistic and one-way sensitivity analyses;
  - voucher and ESRD-conditional methods clarified.
- **Results:** all values updated. The racial-disparity and voucher sections are rewritten, and a new paragraph quantifies the effect of donor-specific inputs.
- **Discussion:** updated values and a new limitations paragraph on outcomes after ESRD.
- **Tables 1–2 and all figures:** regenerated. Table 1 adds rows for the prior-donor inputs.

### Remaining limitations

- The donor-specific inputs rest on a single cohort of 99 donors with ESRD between 1994 and 2011. I have propagated its uncertainty through the probabilistic analysis, but it cannot address era effects.
- The non-donor arm consists of people screened to be as healthy as donors, but it still uses general-population outcomes after ESRD. No source of outcomes for this group is available to me.
- Race-specific waitlist mortality assumes race and age act multiplicatively, because SRTR does not publish rates by both.
