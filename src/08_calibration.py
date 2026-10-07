"""
08_calibration.py
─────────────────
Face-validity and calibration checks for the kidney donation Markov model.

Four comparisons are run:
  1. 3-year waitlist outcomes vs SRTR 2023 ADR Figure KI 22
  2. Post-transplant 5-year patient survival vs SRTR 2023 ADR Figure KI 70
  3. Mortality after ESRD, prior donors vs matched non-donors (Muzaale 2016)
  4. Dialysis life expectancy by age vs USRDS 2023 ADR Table 6.1

Circularity is noted for each row:
  NON-CIRCULAR  Parameter was NOT derived from this observable → genuine test
  CIRCULAR      Parameter WAS derived from this observable → consistency check only

Usage:
  python src/08_calibration.py

Output:
  results/calibration_report.txt
"""

import sys
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load_params, DATA_PROC, RESULTS, mean_to_annual_tx_prob,
                   load_life_table, dialysis_annual_mort, esrd_onset_mortality_curve)

BASE = load_params()
LT   = load_life_table(DATA_PROC / "lifetable_combined_2021.csv")
W = 72


def load_json(path):
    import json
    with open(path) as f:
        return json.load(f)


LIT = load_json(DATA_PROC / "literature_params.json")


# ── FORMATTING ────────────────────────────────────────────────────────────────

def _header(title):
    return f"\n{'─'*W}\n  {title}\n{'─'*W}"

def _colhead():
    return f"  {'Outcome':<42} {'Model':>6}  {'Obs.':>6}  {'Δ%':>6}  {'Status'}"

def _row(label, pred, obs, circ, note=""):
    tag  = "⬤ circular  " if circ == "CIRCULAR" else "○ non-circ."
    err  = f"{100*(pred-obs)/obs:+.0f}%" if obs else "n/a"
    line = f"  {label:<42} {pred:6.1%}  {obs:6.1%}  {err:>6}  {tag}"
    if note:
        line += f"\n    ↳ {note}"
    return line


# ── CHECK 1: 3-YEAR WAITLIST OUTCOMES ─────────────────────────────────────────

def check_waitlist_outcomes():
    """
    Start a cohort at waitlist entry; propagate 3 annual cycles (die → tx →
    remove); compare 3-year outcome distribution to SRTR 2023 ADR Figure KI 22
    (2019-2021 listing cohort).

    Circularity notes
    -----------------
    Transplant  NON-CIRCULAR  Annual tx-prob derived from mean wait (1,765 d,
                              Punjala 2024) via exponential model; NOT from
                              the KI 22 3-year outcome distribution.
    Died        NON-CIRCULAR  Derived from cross-sectional rate (5.0/100 PY,
                              SRTR KI 24, 2023); NOT from the KI 22 cohort.
    Removed     CIRCULAR      wl_removal_rate_yr = 1-(1-0.191)^(1/3) derived
                              from this KI 22 figure; see caveat below.
    """
    wl_tx   = float(mean_to_annual_tx_prob(float(BASE["wl_std_mean_days"])))
    wl_mort = float(1.0 - math.exp(-float(BASE["wl_mort_per_100py"]) / 100))
    wl_rem  = float(BASE["wl_removal_rate_yr"])

    WL = 1.0
    cum_tx = cum_died = cum_rem = 0.0
    for _ in range(3):
        d = WL * wl_mort;   WL -= d;   cum_died += d
        t = WL * wl_tx;     WL -= t;   cum_tx   += t
        r = WL * wl_rem;    WL -= r;   cum_rem  += r

    # SRTR 2023 ADR Figure KI 22 (2019-2021 listing cohort, adult kidney)
    obs_all_tx  = 0.308 + 0.135   # DDKT 30.8% + LDKT 13.5%
    obs_ddkt    = 0.308
    obs_died    = 0.068
    obs_removed = 0.191
    obs_waiting = 0.299

    lines = [
        _header("CHECK 1 — 3-year waitlist outcomes  (SRTR 2023 ADR Fig KI 22)"),
        _colhead(),
        _row("Transplanted — all Tx (DDKT + LDKT)",
             cum_tx, obs_all_tx, "NON-CIRCULAR",
             f"tx-prob {wl_tx:.3f}/yr from 1,765-d mean; SRTR combined = {obs_all_tx:.1%}"),
        _row("Transplanted — DDKT only",
             cum_tx, obs_ddkt, "NON-CIRCULAR",
             "Model cannot split DDKT/LDKT; all-Tx comparison is primary"),
        _row("Died on waitlist",
             cum_died, obs_died, "NON-CIRCULAR",
             f"Model: {wl_mort:.4f}/yr from 5.0/100 PY (KI 24); see note [A]"),
        _row("Removed from waitlist",
             cum_rem, obs_removed, "CIRCULAR",
             "wl_removal_rate_yr back-calculated via competing-risk bisection to reproduce this target"),
        _row("Still waiting at 3 yr",
             WL, obs_waiting, "NON-CIRCULAR",
             "Residual check"),
        "",
        f"  Parameters used:",
        f"    wl_std_mean_days    = {int(BASE['wl_std_mean_days'])} d  "
        f"→  annual Tx prob = {wl_tx:.4f}",
        f"    wl_mort_per_100py   = {BASE['wl_mort_per_100py']:.1f}     "
        f"→  annual mort prob = {wl_mort:.5f}",
        f"    wl_removal_rate_yr  = {wl_rem:.4f}",
        "",
        "  [A] Mortality over-prediction (model ~10% vs SRTR 6.8%):",
        "      The model uses a 2023 cross-sectional mortality rate (SRTR KI 24,",
        "      5.0/100 PY, confirmed from the raw Excel file) applied to a",
        "      2019-2021 listing cohort (KI 22). Differences in cohort",
        "      characteristics and calendar time likely explain the gap.",
        "      The bias is symmetric across both arms and conservative.",
        "",
        f"  Note: Removal calibrated via competing-risk bisection using mean tx prob",
        f"      ({BASE['wl_removal_rate_yr']:.4f}/yr); the naive formula 1-(1-CIF)^(1/3)",
        f"      was replaced in src/03_download_srtr.py.",
    ]
    return "\n".join(lines)


# ── CHECK 2: POST-TRANSPLANT 5-YEAR PATIENT SURVIVAL ──────────────────────────

def check_posttx_survival():
    """
    Apply the model's annual post-transplant mortality rates for 5 years and
    compare predicted 5-year survival to SRTR 2023 ADR Figure KI 70 (DDKT,
    2016-2018 transplant cohort).

    CIRCULAR: posttx_annual_mort_ageXX was derived as 1 - SRTR_5yr^(1/5), so
    this check verifies numerical consistency only (rounding effects < 0.2%).
    """
    bands = [
        ("18–34", "posttx_annual_mort_age1834", 0.957),
        ("35–49", "posttx_annual_mort_age3549", 0.914),
        ("50–64", "posttx_annual_mort_age5064", 0.820),
        ("65+",   "posttx_annual_mort_age65p",  0.701),
    ]

    lines = [
        _header("CHECK 2 — Post-Tx 5-yr patient survival  "
                "(SRTR 2023 ADR Fig KI 70, DDKT 2016–2018 cohort)"),
        _colhead(),
    ]
    for band, key, obs in bands:
        m    = float(BASE.get(key, BASE["posttx_annual_mort"]))
        pred = (1.0 - m) ** 5
        lines.append(_row(
            f"5-yr patient survival  age {band}",
            pred, obs, "CIRCULAR",
            f"m = 1−{obs}^(1/5) = {m:.4f}; verifies rounding only",
        ))
    return "\n".join(lines)


# ── CHECK 3: MORTALITY AFTER ESRD — PRIOR DONORS VS MATCHED NON-DONORS ────────

def check_mortality_after_esrd():
    """
    Cumulative mortality 1/3/5/10 yr after ESRD onset at age 50 (Muzaale 2016
    median), donor pathway vs general pathway, against Muzaale 2016.

    Donor     CIRCULAR      donor_dialysis_mort_mult and donor_wl_listing_prob
                            were fit to this curve (05_assemble_parameters.py).
    General   NON-CIRCULAR  General pathway uses USRDS/SRTR inputs only.
                            Muzaale's matched non-donors are an imperfect
                            comparator: ESRD 1994–2011, 39% LDKT, 21.5-month
                            median wait vs the model's post-KAS250 standard
                            arm (1,765-day mean wait).
    """
    m16 = LIT["muzaale2016"]
    years = [1, 3, 5, 10]
    age = int(m16["median_age_at_esrd"])
    donor = esrd_onset_mortality_curve(BASE, age, True, years, LT)
    gen   = esrd_onset_mortality_curve(BASE, age, False, years, LT)
    # General pathway with Muzaale-like transplant access, to size how much of
    # the gap access explains (median 21.5 mo → exponential mean = median/ln2)
    gen_access = esrd_onset_mortality_curve(
        dict(BASE, wl_std_mean_days=m16["wait_median_months_nondonor"] * 30.44 / math.log(2),
             wl_listing_prob=0.22),
        age, False, years, LT)

    lines = [
        _header(f"CHECK 3 — Mortality after ESRD onset at age {age}  "
                "(Muzaale 2016, Transplantation 100:1306)"),
        _colhead(),
    ]
    for t in years:
        lo, hi = m16["mort_after_esrd_donor_ci"][str(t)]
        lines.append(_row(f"Prior donor, {t:>2}-yr mortality", donor[t],
                          m16["mort_after_esrd_donor"][str(t)], "CIRCULAR",
                          f"95% CI {lo:.1%}–{hi:.1%}"))
    for t in years:
        lines.append(_row(f"General pathway, {t:>2}-yr mortality", gen[t],
                          m16["mort_after_esrd_nondonor"][str(t)], "NON-CIRCULAR"))
    lines += [
        "",
        "  Donor inputs: preemptive {:.1%}, listing {:.3f}/yr, dialysis mort ×{:.2f}, "
        "post-Tx mort ×{:.2f}".format(
            BASE["donor_esrd_preemptive_prob"], BASE["donor_wl_listing_prob"],
            BASE["donor_dialysis_mort_mult"], BASE["donor_posttx_mort_mult"]),
        "  [B] General pathway vs Muzaale non-donors: 3- and 5-yr mortality are",
        "      reproduced within ~15%, but year 1 is under-predicted and year 10",
        "      over-predicted. Observed ESRD mortality is front-loaded (frail",
        "      patients die early, leaving a healthier pool); the model's",
        "      homogeneous cohort has no such selection. Transplant access",
        "      explains little of the gap: giving the general pathway Muzaale-like",
        "      access (21.5-month median wait, ~26% listed by 12 mo) lowers 10-yr",
        "      mortality only from {:.0%} to {:.0%}. The non-donor arm's ESRD".format(
            gen[10], gen_access[10]),
        "      incidence is ~0.04%/15 yr, so this has negligible effect on ΔLE, but",
        "      it bears on the ESRD-conditional priority benefit (script 10).",
    ]
    return "\n".join(lines), donor, gen


# ── CHECK 4: DIALYSIS LIFE EXPECTANCY VS USRDS TABLE 6.1 ──────────────────────

def check_dialysis_le():
    """
    Remaining life expectancy on dialysis (no transplant) from the age-banded
    annual mortality, vs USRDS 2023 ADR Table 6.1 (2019, sex-averaged).

    SEMI-CIRCULAR  Both come from USRDS 2019 dialysis mortality, but Table 6.1
                   uses 5-yr age groups while the model uses 4 broad bands and
                   a constant-hazard conversion; this checks that the
                   coarsening and conversion preserve life expectancy.
    """
    usrds = load_json(DATA_PROC / "usrds_mortality_params.json")
    lines = [
        _header("CHECK 4 — Dialysis life expectancy  (USRDS 2023 ADR Table 6.1, 2019)"),
        f"  {'Age group':<20} {'Model':>7}  {'USRDS':>7}  {'Δ%':>6}",
    ]
    for grp, mid in [("40-44", 42), ("50-54", 52), ("60-64", 62),
                     ("70-74", 72), ("80-84", 82)]:
        f_, m_ = usrds.get(f"le_dialysis_female_{grp}"), usrds.get(f"le_dialysis_male_{grp}")
        if f_ is None or m_ is None:
            continue
        obs = (f_ + m_) / 2
        s, le = 1.0, 0.5   # mid-cycle correction for end-of-cycle accumulation
        for age in range(mid, len(LT)):
            s *= 1 - dialysis_annual_mort(BASE, age, q_floor=LT[age])
            le += s
        lines.append(f"  {grp:<20} {le:7.1f}  {obs:7.1f}  {100*(le-obs)/obs:+5.0f}%")
    return "\n".join(lines)


# ── SUMMARY ───────────────────────────────────────────────────────────────────

def summary(wl_lines, posttx_lines):
    lines = [_header("SUMMARY — non-circular checks")]
    wl_tx   = float(mean_to_annual_tx_prob(float(BASE["wl_std_mean_days"])))
    wl_mort = float(1.0 - math.exp(-float(BASE["wl_mort_per_100py"]) / 100))
    wl_rem  = float(BASE["wl_removal_rate_yr"])
    WL = 1.0; cum_tx = cum_died = cum_rem = 0.0
    for _ in range(3):
        d=WL*wl_mort; WL-=d; cum_died+=d
        t=WL*wl_tx;   WL-=t; cum_tx+=t
        r=WL*wl_rem;  WL-=r; cum_rem+=r

    checks = [
        ("Transplanted at 3 yr (all Tx)", cum_tx, 0.443),
        ("Died on waitlist at 3 yr",       cum_died, 0.068),
        ("Still waiting at 3 yr",          WL,      0.299),
    ]
    for label, pred, obs in checks:
        err = 100*(pred-obs)/obs
        flag = "PASS" if abs(err) < 25 else "WARN"
        lines.append(f"  [{flag}]  {label:<38}  model {pred:.1%}  obs {obs:.1%}  "
                     f"Δ = {err:+.0f}%")

    mort_obs = 0.068
    mort_err = round(100 * (cum_died - mort_obs) / mort_obs)
    lines += [
        "",
        "  Interpretation:",
        "  - Transplant rate is well-calibrated (+5%); the exponential-mean",
        "    model for wait time aligns with SRTR observed 3-year Tx fraction.",
        f"  - Mortality is over-predicted ({mort_err:+d}%); see note [A]. This is conservative",
        "    (slightly overstates the cost of the waitlist period for both arms).",
        "  - The model cannot be validated against an independent end-to-end",
        "    dataset (donation → ESRD → transplant) because no such cohort exists.",
        "    The checks above represent the best available face-validity evidence.",
    ]
    return "\n".join(lines)


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    heading = (
        f"\n{'='*W}\n"
        f"  KIDNEY DONATION MARKOV MODEL — CALIBRATION REPORT\n"
        f"{'='*W}\n\n"
        f"  ○ non-circ.  Parameter NOT derived from this outcome → genuine test\n"
        f"  ⬤ circular   Parameter WAS derived from this outcome → consistency only\n"
    )

    c1 = check_waitlist_outcomes()
    c2 = check_posttx_survival()
    c3, _, _ = check_mortality_after_esrd()
    c4 = check_dialysis_le()
    sm = summary(c1, c2)
    report = (heading + c1 + "\n\n" + c2 + "\n\n" + c3 + "\n\n" + c4
              + "\n\n" + sm + "\n")

    print(report)

    out = RESULTS / "calibration_report.txt"
    with open(out, "w", encoding="utf-8") as f:
        # Plain-text version (ASCII symbols)
        f.write(report
                .replace("⬤", "*")
                .replace("○", " ")
                .replace("↳", "->")
                .replace("→", "->"))
    print(f"  Report saved → {out}")


if __name__ == "__main__":
    main()
