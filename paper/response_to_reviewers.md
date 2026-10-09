# Response to Comments

Thank you for your careful reading. Your comments identified important weaknesses in the model's treatment of outcomes after ESRD. Changes in the revised manuscript appear in red. Below I reproduce each comment, summarize the revisions, and respond in turn.

---

## Comment 1

> The headline rests on one study. The 2.5-year figure comes almost entirely from the Mjøen HR, while the other three cohorts in the O'Keeffe meta-analysis imply very little change in lifespan. Your time-varying ramp is a reasonable way to reconcile them, but it's still an assumption. At HR = 1.0 the cost is 68 days. I'd lead with the range or explicit scenarios rather than a single number, as it is likely to arouse the ire of reviewers.

Fair point. I have reframed the paper accordingly. Specifically:

- **Highlights:** The first bullet now reads: "The life-expectancy cost of kidney donation depends primarily on which published donor mortality hazard ratio (HR) is assumed: under two months at HR~$=1.0$ (no excess mortality) to over four years at the upper confidence bound of the longest-follow-up cohort. Under a time-varying base case that reconciles the existing cohort data (excess mortality emerging only after roughly a decade), the estimate is 2.4~years." The third bullet, which previously stated the range, now makes the research-priority point instead to avoid repetition.

- **Abstract (Results):** The paragraph now opens by explaining the HR disagreement and stating the range (under two months to over four years), then presents the base-case 890-day estimate as the product of a specific, literature-motivated reconciliation.

- **Discussion:** The opening sentence now reads "Under the base-case time-varying hazard assumption, the mean life-expectancy cost…" rather than stating the number without qualification.

---

## Comment 2

> The paper converts HRs from papers >=a decade old into life-years. It's useful for counseling, but it isn't new data. If Mjøen's late excess mortality were real, you'd expect it to have changed donor practice by now, and the paper should address why it hasn't. I'd frame it as "what do these estimates imply?" rather than "what does donation cost?" We are not aware of any large, longer-follow-up donor studies since 2014, so your call for updated cohort data is well taken. Unfortunately, the political reality of linked SRTR data to CMS here in the US is such that it's not something we can do right now (believe, me would have done it already).

I have reframed the contribution as translation, not new estimation: The Introduction's concluding sentence has been revised to make explicit that this paper asks what published risk estimates *imply* for individual donors, rather than claiming to estimate the causal effect of donation de novo. I also now explicitly state "This paper does not adjudicate between the competing HR estimates; it quantifies what is at stake in that adjudication" in the Discussion.

---

## Comment 3

> I think my biggest concern is that the inputs are from the wrong population. Dialysis mortality, listing rates, waitlist mortality, and post-transplant survival come from all transplant candidates, not prior donors. Prior donors with ESKD are closely followed and more likely to be listed preemptively or find a living donor, they are also much healthier (generally). The dialysis mortality figures also aren't age-stratified and are driven by older patients, so they likely overstate harm for a donor-age cohort. Muzaale's paper is a good starting point for donor-specific inputs.

I agree on every point and have revised the model accordingly:
1. **Donor-specific inputs after ESRD.** The donor arm now uses inputs from Muzaale et al. (2016), *Transplantation* 100:1306, which followed 99 US living donors who developed ESRD. Two inputs are taken directly from the study, and two are calibrated to its observed listing and mortality.
2. **Age-stratified dialysis mortality.** I replaced the flat 22% (first year) and 17% (later years) rates with USRDS age-specific rates, which are much lower at donor ages (6.8%/year at ages 18–44).
3. **Age-stratified waitlist mortality.** I replaced the all-ages SRTR rate with SRTR's age-specific rates.

The headline result is largely unchanged: about 890 days (2.4 years) versus 912 days (2.5 years). It remains dominated by the donor all-cause mortality hazard ratio. However, the cost that runs through ESRD fell by about a fifth, the racial disparity narrowed considerably, and the conditional value of priority access changed.

---

## Comment 4

> NKR vouchers give priority for a living-donor kidney, but you model them with the deceased-donor priority wait and DDKT survival. The 1,765-day standard wait also comes from a study of long waiting times, which may skew it upward.

Both points are correct, and I have revised the model accordingly.

1. LDKT survival for the voucher arm.

The NKR voucher mechanism is explicitly designed to facilitate living-donor transplantation via a kidney-exchange chain. Modeling the voucher arm with DDKT post-transplant survival understated its benefit. The revised voucher arm uses LDKT patient survival (SRTR 2023 ADR Figure KI-76: 0.4%, 0.8%, 1.7%, and 3.9%/year at ages 18–34, 35–49, 50–64, and 65+) and LDKT graft survival (SRTR 2023 ADR Figure KI-61), replacing the prior DDKT figures.

The donor arm's base-case retains DDKT survival. Muzaale et al. show that prior donors who develop ESRD receive living-donor kidneys *less* often than matched non-donors (13% vs 39%), predominantly receiving standard-criteria deceased-donor kidneys (87%), consistent with priority allocation making a DDKT available quickly.

This correction substantially increased the conditional value of priority access. LDKT outcomes are meaningfully better than DDKT outcomes at every age, so the life-years added by faster access to a living-donor kidney are larger than those from faster access to a deceased-donor kidney.

**2. Standard waitlist mean wait time.**

I replaced the Punjala 2024 mean of 1,765 days with a value derived directly from SRTR 2023 ADR Figure KI-22 (2018–2020 listing cohort, 3-year outcomes) using the same competing exponential model already applied to the waitlist removal rate.