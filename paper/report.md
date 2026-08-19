# Proofreading Report

**Paper:** Life Expectancy Changes Due to Kidney Donation Across Race, Sex, and Age  
**Date:** 2026-08-19  
**Source:** /home/nick/kidney-donor-le/paper/main.tex  

---

## Scorecard

| Check | Errors | Warnings | Info | Total |
|-------|--------|----------|------|-------|
| 1. Paper Structure | 2 | 3 | 2 | 7 |
| 2. Math Notation | 1 | 2 | 0 | 3 |
| 3. Statistical Relevance | 0 | 2 | 1 | 3 |
| 4. Figures and Tables | 1 | 6 | 1 | 8 |
| 5. Grammar and Style | 4 | 4 | 0 | 8 |
| 6. Abbreviations | 10 | 2 | 0 | 12 |
| **Total** | **18** | **19** | **4** | **41** |

---

## Overall Assessment

This is a technically strong paper with a clear contribution, well-constructed analysis, and internally consistent quantitative results. The writing is precise and the argument flows well. The most urgent issues are in Check 6 (abbreviations): a cluster of major data-source abbreviations — USRDS, SRTR, NKR, DDKT, HR, PSA, CKD, CIF, QALY — are used throughout without formal introduction, which will trigger immediate rejection at many journals. Check 5 contains two outright grammar errors (a double-verb sentence and a pronoun switch from "I" to "we") and a serious citation error ("Grams 2018" is cited for the PSA hazard-ratio prior but no such entry exists in the bibliography). The table uses `\hline` throughout instead of `booktabs`, and four of the five result figures are PNG raster rather than vector format.

### Top Issues to Address

1. **[ERROR] Check 6 — USRDS and SRTR used in abstract without introduction** — quickest rejection trigger
2. **[ERROR] Check 5 — "Grams 2018" cited for PSA HR prior** (design.tex:49) with no \cite{} and no matching bibliography entry — likely should be O'Keeffe et al. 2018
3. **[ERROR] Check 5 — Grammar: "I derived sex-specific donor ESRD rates were derived"** (design.tex:43) — double-verb sentence; also switches "I" → "we" in same paragraph
4. **[ERROR] Check 6 — HR, PSA, DDKT, CKD, NKR, QALY, CIF, ADR used without introduction** — at least 8 abbreviations need "(ABBREV)" parenthetical on first use
5. **[ERROR] Check 4 — All table rules use `\hline` instead of booktabs** — replace with \toprule, \midrule, \bottomrule
6. **[ERROR] Check 5 — Mixed British/American English** — "haemodialysis"/"modelled"/"parameterised"/"favourable" (British) vs. "parameterized"/"counseling" (American) — choose one variant consistently

---

## Check 1: Paper Structure

**Summary:** The paper is a single-author medical journal submission (SAGE) following IMRAD format, which does not conventionally include contribution bullet-lists or section roadmaps. These structural elements are flagged per the checklist, but many findings here reflect genre conventions rather than errors. The introduction has an explicit gap statement and a clear statement of what the paper does. The Discussion is well-structured, with quantitative results, contextualization, and limitations. The main avoidable issues are the absence of any novelty framing and the word "contribut..." never appearing.

### Findings

- **[ERROR] background.tex — The word "contribut" (contribute/contribution) does not appear anywhere in the introduction.**  
  > "This paper addresses that gap using a Markov microsimulation model..."  
  Suggestion: Add a sentence such as "The specific contributions of this paper are: (1) we quantify the individual life-expectancy cost of donation stratified by age, race, and sex; and (2) we estimate the life-expectancy value of priority-waitlist mechanisms both unconditionally and conditional on ESRD onset." *Note: standard medical papers do not require this, but the journal may.*
Response: This is not a problem for this journal. Skip.

- **[ERROR] background.tex — No formal contribution list in the introduction.**  
  The introduction ends with a description of what the paper does ("This paper addresses that gap...") but does not enumerate contributions as a list.  
  Suggestion: Convert the final paragraph into a bulleted or numbered list of specific, falsifiable contributions.
Response: Skip

- **[WARN] background.tex — No explicit novelty statement ("To the best of our knowledge, this is the first...").**  
  The introduction describes the gap but does not claim priority of this analysis.  
  Suggestion: Add "To the best of our knowledge, this is the first individual-level life-expectancy analysis of kidney donation stratified simultaneously by age, race, and sex" in the introduction's final paragraph.
Response: Skip.

- **[WARN] results.tex — No explicit hypotheses stated before results.**  
  The Results section proceeds directly to subgroup findings without stating expected directions. *Note: formal hypotheses are not required in medical modeling papers, but their absence may draw reviewer comment.*  
  Suggestion: Add brief a priori hypothesis statements (e.g., "We hypothesized that younger and Black donors would bear a disproportionate life-expectancy cost") before the results.
Resposne: Skip. The instruction from the journal is "Always end the Introduction section with a clear statement of the study’s objectives or hypotheses."

- **[WARN] conc.tex — Discussion does not explicitly frame results as improvements over prior work.**  
  The Discussion quantifies results and discusses implications, but does not note that no prior study directly quantified LE cost of donation across demographics. A direct statement would strengthen the novelty framing.  
  Suggestion: Add a sentence in the first Discussion paragraph: "Prior to this analysis, no published study had translated post-donation ESRD risk into individual-level life-expectancy units stratified by age, race, and sex simultaneously."
Resposne: Skip.

- **[INFO] background.tex — No section roadmap at the end of the introduction.**  
  Many SAGE medical journals do not require a roadmap, but consider adding one if the target journal does.
Resposne: Skip.

- **[INFO] abstract.tex:5 — "recipient" is ambiguous in the abstract.**  
  > "...conditional on the recipient's actually developing ESRD..."  
  "Recipient" in kidney-transplant literature means the person receiving the kidney; here it refers to the designated voucher beneficiary (a family member who may never donate). Readers may be confused.  
  Suggestion: Replace "the recipient's actually developing ESRD" with "the designated beneficiary's actually developing ESRD" or "the voucher holder's actually developing ESRD."
Resposne: Good point. Fix.

---

## Check 2: Math Symbols and Notation

**Detected notation convention: Minimal (epidemiological modeling paper with scalar quantities only; no vectors or matrices).** The paper uses a small number of mathematical symbols, and the notation is generally adequate. The main issues are one undefined formula symbol and an inconsistent typesetting of the primary outcome variable.

### Findings

- **[ERROR] design.tex:66 — $\mathrm{CIF}$ used in a formula without definition.**  
  > `$1-(1-\mathrm{CIF})^{1/3}$`  
  The symbol CIF (cumulative incidence function) appears inside a formula without prior prose definition of the form "where CIF denotes..." The abbreviation CIF also lacks a formal "(CIF)" introduction elsewhere in the text (see Check 6).  
  Suggestion: Precede with "...where $\mathrm{CIF}$ denotes the three-year cumulative incidence..." or replace with a descriptive expression.
Response: I have now fixed this.

- **[WARN] design.tex:29 — $S_5$ and $S_1$ used without explicit definition.**  
  > `derived as $1-(S_5/S_1)^{1/4}$ from age-specific 1- and 5-year DDKT graft survival`  
  The subscripts are recoverable from context ("1- and 5-year...graft survival"), but $S_n$ is never formally defined as the $n$-year survival probability.  
  Suggestion: Add "where $S_n$ denotes the $n$-year graft survival probability" immediately after the formula.
Response: Good point. Fix.

- **[WARN] Multiple locations — Inconsistent typesetting of the primary outcome symbol.**  
  `$\Delta\text{LE}$` (LE inside `\text{}`, rendered as upright letters) appears in results.tex:5,17,21 and tables.tex:7, but `$\Delta$LE` (LE outside math mode, rendered as ordinary text but adjacent to math) appears in results.tex:21 ("shifted $\Delta$LE by up to 15~days") and figures.tex:28 ("$\Delta$LE resulting from testing"). These produce visually different output.  
  Suggestion: Pick one form and use it throughout. Prefer `$\Delta\text{LE}$` for consistency with the first use in results.tex.
Response: Agreed. Fix.

---

## Check 3: Statistical Relevance

**Summary:** This is primarily a deterministic Markov model paper. A PSA (200 draws) is presented, and its distributional output is shown in fig:psa as a histogram — appropriate for communicating uncertainty. Deterministic OWSA results in the tornado chart and ESRD-conditional sweep are correctly presented as single-valued outputs. The age-race-sex matrix similarly shows deterministic scenario outputs. The main finding in this check is a comparative claim in the introduction that lacks in-sentence quantitative support.

### Findings

- **[WARN] background.tex:1 — Comparative claim without quantitative support in the sentence.**  
  > "Living donor kidney transplantation (LDKT) remains the treatment of choice for end-stage renal disease (ESRD), offering better graft function and longer allograft survival than deceased-donor transplantation."  
  The claim is supported by a citation but no numbers appear in the sentence itself. The comparison motivates the paper's framing.  
  Suggestion: Add concrete numbers, e.g., "...offering better graft function (10-year graft survival ~70% vs. ~55% for deceased-donor recipients~\cite{srtr_2023}) and longer allograft survival."
Response: Skip.

- **[WARN] figures/kidney_model_age_race_sex_matrix.png — Results matrix shows no uncertainty.**  
  The matrix presents 21 deterministic ΔLE values with no confidence intervals or PSA-derived ranges. Because the PSA (fig:psa) shows the credible interval is very wide (±~1,000 days) under HR uncertainty, readers seeing the matrix alone may over-interpret the precision of the point estimates.  
  Suggestion: Add a note to the caption or nearby text that PSA uncertainty for these subgroup estimates is also dominated by the HR prior; or add a brief supplementary matrix showing PSA credible intervals.
Response: Skip.

- **[INFO] results.tex:17 — PSA mean ΔLE (+25.8 days) is positive while base case is −40.5 days.**  
  The paper explains this as a consequence of HR prior draws below 1.0, but some readers may be surprised that the PSA mean and base case have opposite signs. The explanation is adequate; consider making it more prominent in the text.
Response: I re-rean the PSA, and the mean Delta LE is approximately 0. Report that.

---

## Check 4: Figures and Tables

**Summary:** All references to figures and tables are valid and no dangling labels were found. The main structural issue is pervasive use of `\hline` in the table. Four of five result figures are in raster (PNG) format. The Markov model figure caption is a single minimal sentence. Float placement uses `[ht]` throughout.

### Findings

- **[ERROR] tables.tex:9,11,17,21,25,29,32,36 — `\hline` used throughout; booktabs not used.**  
  > `\hline` (appears 8 times in the table environment)  
  Suggestion: Replace with `\toprule` (top), `\midrule` (between header and body rows), and `\bottomrule` (bottom). Also add `\usepackage{booktabs}` to main.tex preamble. The current table has two `\hline` rows for every section separator — replace each pair with a single `\midrule`.
Response: Agreed. Implement booktabs.

- **[WARN] figures/kidney_model_age_race_sex_matrix.png — Raster format (PNG).**  
  Suggestion: Export as PDF or SVG. The heatmap with many text labels will pixelate when scaled in print.
Response: Skip

- **[WARN] figures/kidney_model_B_psa.png — Raster format (PNG).**  
  Suggestion: Export as PDF or SVG.
Response: Skip

- **[WARN] figures/kidney_model_C_tornado.png — Raster format (PNG).**  
  Suggestion: Export as PDF or SVG.
Response: Skip

- **[WARN] figures/esrd_conditional_age_sweep.png — Raster format (PNG).**  
  Suggestion: Export as PDF or SVG.
Response: Skip

- **[WARN] figures.tex:4,11,18,25,32; tables.tex:4 — All figure and table environments use `[ht]` placement, which includes the discouraged `h` (here) option.**  
  Suggestion: Replace `[ht]` with `[t]` (or `[!t]`) throughout to let LaTeX choose the best top-of-page slot without the "here" fallback that can interrupt text flow.
Response: Skip. The figures are in their own section at the end.

- **[WARN] figures.tex:7 — fig:markov caption is too brief and not self-contained.**  
  > `\caption{The Markov model for kidney donation.}`  
  A reader encountering the figure without the surrounding text would not know what health states are shown, what the arrows represent, or how the donor and non-donor arms differ.  
  Suggestion: Expand to something like: "Markov state-transition model for kidney donation. Boxes represent the seven mutually exclusive health states; arrows show permitted transitions. The donor arm (upper pathway) routes ESRD-onset patients to the Priority waitlist; the non-donor arm routes to the Standard waitlist. The Dead state is absorbing."
Response: Good. Implement.

- **[WARN] tables.tex:8 — Numeric columns use `r` (right-aligned) rather than `c` (center-aligned).**  
  > `\begin{tabular}{lrr}`  
  Suggestion: Change to `{lcc}` and verify that decimal alignment is adequate; alternatively use `S` from the `siunitx` package to align on the decimal point.
Response: Agreed. Implement right alignment.

- **[WARN] figures/esrd_conditional_age_sweep.png — x-axis label "Age at ESRD onset" has no unit.**  
  The axis shows integer tick values (40, 45, ..., 75) without stating years.  
  Suggestion: Change x-axis label to "Age at ESRD onset (years)".
Response: Skip. I have fixed.

- **[INFO] figures/kidney_model_C_tornado.png; figures/esrd_conditional_age_sweep.png — No background grid lines.**  
  Light grid lines would aid reading values off the axes.
Response: Skip.

---

## Check 5: Grammar and Style

**Detected English variant: Mixed — predominantly British ("haemodialysis," "modelled," "parameterised," "favourable") with American intrusions ("parameterized," "counseling"). This is an [ERROR] — see finding 5.3 below.**

### Findings

- **[ERROR] design.tex:43 — Double-verb sentence structure.**  
  > "I derived sex-specific donor ESRD rates were derived by applying the Massie~2017 male-sex HR..."  
  This sentence fuses two constructions ("I derived X" and "X were derived by...") and is ungrammatical.  
  Suggestion: "Sex-specific donor ESRD rates were derived by applying the Massie~2017 male-sex HR (1.88)~\cite{massie_quantifying_2017} to the overall donor rate, weighted by the observed donor sex mix (63.5\% female; SRTR 2023 ADR Figure~KI~104).\cite{srtr_2023}"
Response: Make it first-person singular and active voice. "I derived X"

- **[ERROR] design.tex:43 — Pronoun switch from "I" to "we" within the same paragraph.**  
  > "For non-donor baselines, **we** used the Grams~2016 White sex-specific rates"  
  The rest of the paper (and this section) uses first-person singular "I."  
  Suggestion: Change "we used" to "I used."
Response: Agreed. Fix.

- **[ERROR] design.tex:49 — "Grams~2018 meta-analysis CI" cited without a `\cite{}` and with no matching bibliography entry.**  
  > "Donor all-cause mortality HR: log-normal, $\mu = {-0.016}$, $\sigma = 0.143$, derived from the Grams~2018 meta-analysis CI;"  
  The bibliography contains `grams_kidney-failure_2016` (Grams et al. 2016, a kidney-failure risk paper) but no Grams 2018 entry. The donor all-cause mortality HR uncertainty is discussed elsewhere in the paper using the O'Keeffe et al.\ 2018 meta-analysis (`okeeffe_mid-_2018`). The log-normal parameters (μ = −0.016, σ = 0.143) imply a 95% CI of approximately (0.74, 1.30), which does not match the O'Keeffe pooled RR 0.60 (CI 0.31–1.10) but could reflect a different analysis. This appears to be a citation error.  
  Suggestion: Verify which study provided these log-normal parameters and add the correct `\cite{}` reference. If the source is O'Keeffe 2018 with a restricted prior, explain the derivation explicitly. If a separate Grams 2018 meta-analysis exists, add it to references.bib.
Response: It should be O'Keeffe. Please fix with the correct values.

- **[ERROR] Multiple locations — Mixed British/American English.**  
  British forms: "haemodialysis" (design.tex:10,18; conc.tex:3), "modelled" (abstract.tex:5, design.tex:26, conc.tex:11, acks\_declarations.tex:23), "parameterised" (design.tex:48), "parameterisation" (design.tex:26), "favourable" (conc.tex:3).  
  American forms: "parameterized" (abstract.tex:3), "counseling" (conc.tex:5).  
  Most notably, "parameterized" and "parameterised" both appear in the same document.  
  Suggestion: Choose British English throughout (consistent with "haemodialysis" and "modelled"). Change "parameterized" → "parameterised" (abstract.tex:3) and "counseling" → "counselling" (conc.tex:5).
Response: Fix to use American English.

- **[WARN] design.tex:66 — Unclosed parenthesis.**  
  > "The model predicted 38.3\% transplanted at three years versus 44.3\% observed (DDKT + LDKT combined. No independent cohort..."  
  The opening parenthesis "(" is never closed before the full stop.  
  Suggestion: Change to "...versus 44.3\% observed (DDKT + LDKT combined). No independent cohort..."
Response: Agreed. Fix.

- **[WARN] figures.tex:28 — Em-dash in figure caption.**  
  > "\textbf{One-way sensitivity analysis — tornado chart.}"  
  Suggestion: Replace the em-dash with a colon: "\textbf{One-way sensitivity analysis: tornado chart.}" — or move the subtitle into the caption body.
Response: Agreed. Fix.

- **[WARN] background.tex:1 — Comparative claim without quantitative support in the sentence** (see also Check 3).  
  > "offering better graft function and longer allograft survival than deceased-donor transplantation"  
  Suggestion: Add specific numbers or cite values inline.
Response: Skip.

- **[WARN] abstract.tex:1 — First-person singular "I" is used throughout the abstract and methods.**  
  > "I constructed a Markov microsimulation model..."  
  *Note: For a single-author medical paper, "I" is standard practice and many SAGE journals accept it. Flag is included per checklist rule but no change may be needed — confirm with the target journal's style guide.*
Response: Should be "I" throughout.

---

## Check 6: Abbreviations

**Summary:** This check identified the most impactful issues in the paper. At least ten abbreviations are used without formal introduction in the body text, including the two primary data sources (USRDS, SRTR) that appear un-introduced in the very first methods sentence. Several are also used in the abstract without introduction. Fixing these is critical before submission.

### 6.2 — Abstract

- **[ERROR] abstract.tex:3 — "USRDS" used without introduction in the abstract.**  
  > "parameterized from USRDS, SRTR, and published cohort studies..."  
  Suggestion: Change to "parameterized from the United States Renal Data System (USRDS), the Scientific Registry of Transplant Recipients (SRTR), and published cohort studies..."

- **[ERROR] abstract.tex:3 — "SRTR" used without introduction in the abstract.**  
  See above. Fix jointly with USRDS.

### 6.3 — First use in body

- **[ERROR] background.tex:7 — "CKD" used without introduction.**  
  > "progression to CKD stage~3 or higher"  
  Suggestion: Change to "progression to chronic kidney disease (CKD) stage~3 or higher."

- **[ERROR] background.tex:11 — "NKR" used without introduction.**  
  > "the NKR voucher"  
  "National Kidney Registry" appears earlier (background.tex:9) without "(NKR)."  
  Suggestion: Add "(NKR)" the first time "National Kidney Registry" appears: "The National Kidney Registry (NKR) has extended this priority concept..."

- **[ERROR] design.tex:18 — "USRDS" used without introduction in body text.**  
  > "consistent with USRDS estimates"  
  Full form "United States Renal Data System" never appears in the body.  
  Suggestion: Introduce on first use: "United States Renal Data System (USRDS)" (likely already fixed if abstract introduction is added per 6.2 above; but body still needs its own introduction per convention).

- **[ERROR] design.tex:20 — "ADR" (Annual Data Report) used without introduction.**  
  > "USRDS 2025 ADR Figure~7.13"  
  Suggestion: First use: "USRDS 2025 Annual Data Report (ADR) Figure~7.13" and use "ADR" subsequently.

- **[ERROR] design.tex:29 — "DDKT" (deceased-donor kidney transplantation) used without introduction.**  
  > "Kaplan--Meier estimates for DDKT recipients"  
  "Deceased-donor transplantation" appears in background.tex without "(DDKT)."  
  Suggestion: In background.tex where deceased-donor transplantation first appears: "...than deceased-donor kidney transplantation (DDKT)." Then use "DDKT" consistently in the methods.

- **[ERROR] design.tex:32 — "HR" (hazard ratio) used without formal introduction.**  
  > "I also examined an HR of 1.30"  
  "Hazard ratio" is written out earlier (background.tex:5, design.tex:15) but never parenthetically introduces "(HR)."  
  Suggestion: At the first prose use of "hazard ratio" in background.tex: "...a hazard ratio (HR) of 1.30" — then use "HR" throughout.

- **[ERROR] design.tex:46 — "PSA" used without formal introduction in body text.**  
  > "I conducted a PSA with 200 parameter draws."  
  The subsection heading "Probabilistic sensitivity analysis" does not constitute a parenthetical introduction.  
  Suggestion: Change the opening sentence to: "I conducted a probabilistic sensitivity analysis (PSA) with 200 parameter draws."

- **[ERROR] design.tex:66 — "CIF" used without introduction.**  
  > "$1-(1-\mathrm{CIF})^{1/3}$" and "Grams race-specific non-donor CIFs"  
  Suggestion: Introduce on first use: "...the naïve formula $1-(1-\mathrm{CIF})^{1/3}$, where CIF denotes the three-year cumulative incidence fraction..." and use "CIF" thereafter.

- **[ERROR] conc.tex:12 — "QALY" used without introduction.**  
  > "a QALY-adjusted analysis would likely widen..."  
  Suggestion: Change to "a quality-adjusted life-year (QALY)-adjusted analysis" on first use.

- **[ERROR] acks_declarations.tex:23 — "NHANES" used without introduction.**  
  > "drawn from NHANES~III"  
  Suggestion: Change to "the Third National Health and Nutrition Examination Survey (NHANES~III)."

- **[ERROR] design.tex:15 — "CI" (confidence interval) used without formal introduction.**  
  > "2.96 (95\%~CI 2.25--3.89)"  
  *Note: "CI" for confidence interval is extremely common in medical literature and many journals treat it as a standard notation requiring no introduction. Verify whether the target journal requires explicit introduction.*  
  Suggestion: Add "(CI)" on first use of "confidence interval" if the journal requires it: "95% confidence interval (CI)."

### 6.4 — Subsequent uses

- **[WARN] design.tex:29 — "SRTR" used in body without introduction, then used repeatedly.**  
  "SRTR" appears throughout the methods (design.tex:26,29,32,37,40,43) without ever being introduced as "Scientific Registry of Transplant Recipients (SRTR)." This is likely the most-cited data source in the paper.  
  Suggestion: Introduce on first use in background.tex or at design.tex:18: "Scientific Registry of Transplant Recipients (SRTR)."

- **[WARN] figures.tex:28 (fig:tornado caption) — "HR" used in figure caption as if previously defined.**  
  > "The donor all-cause mortality HR (top bar, truncated at left)..."  
  This is acceptable once HR is introduced in the body per the fix for 6.3 above. No separate action needed beyond fixing the body introduction.
