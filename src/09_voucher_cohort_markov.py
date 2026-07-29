"""
09_voucher_cohort_markov.py
============================
Analytic cohort Markov model quantifying the life-expectancy benefit of a
kidney-donation "voucher" held by a healthy, non-donor family member.

Background: Living kidney donors can designate non-donor family members to
receive a voucher granting them the same priority waitlist access as prior
living donors (PLD) — approximately 100-day median wait vs 985 days on the
standard deceased-donor waitlist — should those family members ever develop
ESRD and need a transplant.

Research question: How much does holding such a voucher improve life expectancy
for a healthy non-donor, compared to an otherwise identical person without one?

Model design:
  Arm A  Voucher holder — non-donor ESRD risk, priority waitlist if ESRD develops
  Arm B  Control       — non-donor ESRD risk, standard waitlist if ESRD develops

Both arms carry baseline (non-donor) ESRD risk and background mortality HR = 1.0.
The sole difference between arms is waitlist priority upon ESRD onset.
Post-transplant outcomes use DDKT survival (base case) since the UNOS PLD
priority is for deceased-donor allocation; LDKT outcomes tested in sensitivity.

Non-Markovian dialysis mortality is handled by state-splitting (same approach
as 07_cohort_markov.py):
  D1  ESRD / dialysis year 1  (22% annual mortality)
  D2  ESRD / dialysis year 2+ (17% annual mortality)

States:
  H    Healthy
  D1   ESRD year 1
  D2   ESRD year 2+
  WL   Waitlist (priority for voucher arm, standard for control arm)
  PT   Post-transplant
  Dead Absorbing (tracked implicitly)

Sensitivity analyses performed:
  1. Voucher wait time:             50 d / 102.6 d (base) / 200 d
  2. Preemptive listing probability: 5.8% (conservative) / 9.4% (base, informed)
  3. Per-cycle listing rate:         0.05 / 0.15 (base) / 0.30
  4. Post-Tx quality:               DDKT (base) / LDKT (living-donor chain)
  5. Non-donor ESRD risk:           overall / Black-race / White-race

Usage:
  python src/09_voucher_cohort_markov.py
  python src/09_voucher_cohort_markov.py --age 25
  python src/09_voucher_cohort_markov.py --age 40 --n 500000
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load_life_table, load_params, DATA_PROC, RESULTS,
                   weibull_annual_prob, weibull_scale_from_cumrisk_competing,
                   median_to_annual_tx_prob)

# ── CONSTANTS ─────────────────────────────────────────────────────────────────
N_PER_ARM = 1_000_000
MAX_AGE   = 100
CYCLE_YRS = 1.0

_lt_path = DATA_PROC / "lifetable_combined_2021.csv"
LIFE_TABLE_QX = load_life_table(_lt_path if _lt_path.exists() else None)
BASE = load_params()

# Color palette
_BLUE  = "#2B6CB0"   # voucher holder
_GREY  = "#718096"   # control (no voucher)
_LIGHT = "#F1EFE8"
_STATE_COLORS = {
    "H":  "#4DB8A0",
    "D1": "#E87A5D",
    "D2": "#C4503A",
    "WL": "#534AB7",
    "PT": "#BA7517",
}


# ── HELPERS ───────────────────────────────────────────────────────────────────
def _ptx_mort(p: dict, age: int, ldkt: bool = False) -> float:
    """Age-stratified annual post-Tx mortality. DDKT base; LDKT for sensitivity."""
    prefix = "posttx_ld_" if ldkt else "posttx_"
    fallback = float(p.get("posttx_ld_annual_mort" if ldkt else "posttx_annual_mort",
                           p["posttx_annual_mort"]))
    if age < 35:   return float(p.get(f"{prefix}annual_mort_age1834", fallback))
    elif age < 50: return float(p.get(f"{prefix}annual_mort_age3549", fallback))
    elif age < 65: return float(p.get(f"{prefix}annual_mort_age5064", fallback))
    else:          return float(p.get(f"{prefix}annual_mort_age65p",  fallback))


def _wl_mort(p: dict) -> float:
    return float(1.0 - np.exp(-p["wl_mort_per_100py"] / 100))


# ── COHORT SIMULATION ─────────────────────────────────────────────────────────
def run_arm(p: dict, n: float, age_at_entry: int, voucher: bool,
            esrd_cum_risk_15: float = None, ldkt: bool = False,
            life_table: np.ndarray = None):
    """
    Propagate a cohort of healthy non-donors forward through annual Markov cycles.

    Parameters
    ----------
    p : dict
        Model parameters (from params.json).
    n : float
        Initial cohort size (result is independent of this value).
    age_at_entry : int
        Age at cohort entry (years).
    voucher : bool
        True  → priority waitlist (~100-day median) upon ESRD onset.
        False → standard waitlist (~985-day median) upon ESRD onset.
    esrd_cum_risk_15 : float, optional
        15-year competing-risk ESRD CIF. Defaults to non-donor overall.
    ldkt : bool
        True → use LDKT post-Tx mortality (sensitivity). Default: DDKT.
    life_table : np.ndarray, optional
        Overrides the module-level life table.

    Returns
    -------
    mean_le : float
        Mean remaining life-years from age_at_entry.
    state_trace : list of dict
        Per-cycle state occupancy counts.
    """
    lt = life_table if life_table is not None else LIFE_TABLE_QX

    median_days = float(p["wl_pld_median_days"] if voucher else p["wl_std_median_days"])
    wl_tx      = float(median_to_annual_tx_prob(median_days))
    wl_mort    = _wl_mort(p)
    wl_remove  = float(p["wl_removal_rate_yr"])
    wl_listing = float(p.get("wl_listing_prob", 0.15))
    dial_mort1 = float(p["dialysis_1yr_mort"])
    dial_mort  = float(p["dialysis_annual_mort"])
    graft_fail = float(p.get("graft_annual_fail_postyear1", 0.025))

    # Voucher holders are part of the living-donation ecosystem; they're informed
    # and monitored, making preemptive listing more likely (base: 9.4%, same as
    # donor-like patients). Control uses the general ESRD population rate (5.8%).
    preemptive_p = float(p.get(
        "esrd_preemptive_prob_pld" if voucher else "esrd_preemptive_prob_std", 0.0
    ))

    # Both arms carry non-donor (baseline) ESRD risk
    cum_risk_15 = esrd_cum_risk_15 if esrd_cum_risk_15 is not None \
                  else float(p["esrd_15yr_nondonor"])
    wbl_k   = float(p["weibull_shape"])
    wbl_lam = weibull_scale_from_cumrisk_competing(
        cum_risk_15, wbl_k, lt, age_at_entry, bg_hr=1.0
    )

    H  = float(n)
    D1 = 0.0
    D2 = 0.0
    WL = 0.0
    PT = 0.0

    total_ly   = 0.0
    state_trace = []
    yr = 0

    while H + D1 + D2 + WL + PT > 0.5:
        age   = age_at_entry + yr
        q_bg  = lt[min(age, MAX_AGE)]          # background mortality (HR=1.0 both arms)
        p_esrd = weibull_annual_prob(float(yr), wbl_lam, wbl_k)
        ptx_m  = _ptx_mort(p, age, ldkt=ldkt)

        # H (Healthy)
        H_die        = H * q_bg
        H_surv       = H - H_die
        H_esrd       = H_surv * p_esrd
        H_preemptive = H_esrd * preemptive_p   # → WL directly (no dialysis)
        H_to_D1      = H_esrd - H_preemptive   # → D1 (start dialysis)
        H_stay       = H_surv - H_esrd

        # D1 (ESRD year 1) — all survivors exit D1 after exactly one cycle
        D1_die    = D1 * dial_mort1
        D1_surv   = D1 - D1_die
        D1_listed = D1_surv * wl_listing       # → WL
        D1_to_D2  = D1_surv - D1_listed        # → D2 (survived year 1, not yet listed)

        # D2 (ESRD year 2+)
        D2_die    = D2 * dial_mort
        D2_surv   = D2 - D2_die
        D2_listed = D2_surv * wl_listing       # → WL
        D2_stay   = D2_surv - D2_listed        # → D2

        # WL (Waitlist)
        WL_die    = WL * wl_mort
        WL_surv   = WL - WL_die
        WL_tx     = WL_surv * wl_tx            # → PT
        WL_after  = WL_surv - WL_tx
        WL_remove = WL_after * wl_remove       # → D2 (removed, returns to dialysis)
        WL_stay   = WL_after - WL_remove

        # PT (Post-transplant)
        PT_die  = PT * ptx_m
        PT_surv = PT - PT_die
        PT_fail = PT_surv * graft_fail         # → D2 (graft failure)
        PT_stay = PT_surv - PT_fail

        # Update state counts
        H  = H_stay
        D1 = H_to_D1
        D2 = D1_to_D2 + D2_stay + WL_remove + PT_fail
        WL = D1_listed + D2_listed + WL_stay + H_preemptive
        PT = WL_tx + PT_stay

        alive = H + D1 + D2 + WL + PT
        total_ly += alive * CYCLE_YRS

        state_trace.append({
            "yr": yr, "age": age,
            "H": H, "D1": D1, "D2": D2, "WL": WL, "PT": PT,
            "alive": alive,
        })

        yr += 1
        if yr > MAX_AGE - age_at_entry + 10:
            break

    return total_ly / n, state_trace


# ── PLOTTING ──────────────────────────────────────────────────────────────────
def make_fig_state_occupancy(trace_v, trace_c, n, age_at_entry):
    """Fraction of cohort in each state over time, voucher vs control."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    fig.patch.set_facecolor("#FAFAF8")

    for ax, trace, title in zip(
        axes, [trace_v, trace_c], ["Voucher holder (priority)", "Control (standard)"]
    ):
        ages = [t["age"] for t in trace]
        for key, color in _STATE_COLORS.items():
            vals = [t[key] / n for t in trace]
            ax.plot(ages, vals, color=color, lw=1.8, label=key)
        ax.set_facecolor(_LIGHT)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#B4B2A9")
        ax.tick_params(colors="#5F5E5A", labelsize=9)
        ax.set_xlabel("Age", fontsize=9)
        ax.set_ylabel("Fraction of cohort", fontsize=9)
        ax.set_title(title, fontsize=11, fontweight="bold", color="#2C2C2A")
        ax.legend(fontsize=9, frameon=False)

    fig.suptitle(
        "State occupancy — voucher-holder cohort Markov (non-donor ESRD risk)",
        fontsize=12, fontweight="bold", color="#2C2C2A"
    )
    fig.tight_layout()
    return fig


def make_fig_survival(trace_v, trace_c, n, age_at_entry):
    """Survival curves: voucher holder vs control."""
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor("#FAFAF8")
    ax.set_facecolor(_LIGHT)

    ages_v = [age_at_entry] + [t["age"] for t in trace_v]
    surv_v = [1.0] + [t["alive"] / n for t in trace_v]
    ages_c = [age_at_entry] + [t["age"] for t in trace_c]
    surv_c = [1.0] + [t["alive"] / n for t in trace_c]

    ax.plot(ages_v, surv_v, color=_BLUE, lw=2.0, label="Voucher holder (priority WL)")
    ax.plot(ages_c, surv_c, color=_GREY, lw=2.0, label="Control (standard WL)",
            linestyle="--")
    ax.legend(fontsize=10, frameon=False)
    ax.set_xlabel("Age", fontsize=10)
    ax.set_ylabel("Survival fraction", fontsize=10)
    ax.set_ylim(0, 1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#B4B2A9")
    ax.tick_params(colors="#5F5E5A", labelsize=9)
    ax.set_title("Survival curves — voucher holder vs control",
                 fontsize=12, fontweight="bold", color="#2C2C2A")
    fig.tight_layout()
    return fig


def make_fig_tornado(owsa_results, base_diff_days):
    """One-way sensitivity analysis tornado plot."""
    labels = [r["label"] for r in owsa_results]
    lows   = [r["low_days"]  - base_diff_days for r in owsa_results]
    highs  = [r["high_days"] - base_diff_days for r in owsa_results]

    # Sort by total swing (descending)
    swings = [abs(h - l) for h, l in zip(highs, lows)]
    order  = sorted(range(len(swings)), key=lambda i: swings[i])
    labels = [labels[i] for i in order]
    lows   = [lows[i]   for i in order]
    highs  = [highs[i]  for i in order]

    fig, ax = plt.subplots(figsize=(10, max(4, len(labels) * 0.55 + 1.5)))
    fig.patch.set_facecolor("#FAFAF8")
    ax.set_facecolor(_LIGHT)

    y_pos = np.arange(len(labels))
    for i, (lo, hi) in enumerate(zip(lows, highs)):
        # Draw two segments both anchored at x=0 so the baseline always bisects
        # or touches every bar, even when both deviations are on the same side.
        # lo segment (lighter) extends from 0 to lo_delta (leftward when negative)
        # hi segment (darker)  extends from 0 to hi_delta (rightward when positive)
        # matplotlib accepts negative widths: barh(y, w<0, left=0) → [w, 0]
        ax.barh(i, lo, left=0, height=0.5, color=_BLUE, alpha=0.45)
        ax.barh(i, hi, left=0, height=0.5, color=_BLUE, alpha=0.80)
    ax.axvline(0, color="#2C2C2A", lw=1.2, linestyle="-")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Change in ΔLE from base case (days)", fontsize=10)
    ax.set_title(
        "One-way sensitivity — voucher LE benefit\n"
        f"(base ΔLE = +{base_diff_days:.1f} days)",
        fontsize=11, fontweight="bold", color="#2C2C2A"
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#B4B2A9")
    ax.tick_params(colors="#5F5E5A", labelsize=9)
    fig.tight_layout()
    return fig


def make_fig_age_sweep(age_results):
    """ΔLE (voucher benefit) vs age at voucher designation."""
    ages  = [r["age"] for r in age_results]
    diffs = [r["diff_days"] for r in age_results]

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor("#FAFAF8")
    ax.set_facecolor(_LIGHT)

    ax.plot(ages, diffs, color=_BLUE, lw=2.2, marker="o", markersize=6)

    for age, d in zip(ages, diffs):
        ax.annotate(f"+{d:.1f}d", xy=(age, d), xytext=(2, 6),
                    textcoords="offset points", fontsize=8, color="#2C2C2A")

    ax.set_ylim(bottom=0)
    ax.set_xlabel("Age at voucher designation", fontsize=10)
    ax.set_ylabel("LE benefit of voucher (days)", fontsize=10)
    ax.set_title(
        "Voucher life-expectancy benefit by age\n"
        "(non-donor ESRD risk, DDKT outcomes, base-case parameters)",
        fontsize=11, fontweight="bold", color="#2C2C2A"
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#B4B2A9")
    ax.tick_params(colors="#5F5E5A", labelsize=9)
    fig.tight_layout()
    return fig


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main(age_at_entry: int = 40, n: int = N_PER_ARM):
    print("=" * 60)
    print("VOUCHER-HOLDER COHORT MARKOV MODEL (analytic, deterministic)")
    print("=" * 60)
    print(f"Age at entry    : {age_at_entry}")
    print(f"n per arm       : {n:,}  (2× total)")
    print(f"ESRD 15-yr CIF  : {BASE['esrd_15yr_nondonor']*100:.4f}%  (non-donor baseline)")
    print()

    p = BASE.copy()

    # ── BASE CASE ─────────────────────────────────────────────────────────────
    le_v,  trace_v = run_arm(p, n, age_at_entry, voucher=True)
    le_c,  trace_c = run_arm(p, n, age_at_entry, voucher=False)
    diff = le_v - le_c

    print("BASE CASE RESULTS")
    print("-" * 40)
    print(f"  LE voucher holder : {le_v:.6f}  remaining years (from age {age_at_entry})")
    print(f"  LE control        : {le_c:.6f}  remaining years")
    print(f"  ΔLE (voucher ben.): {diff:+.6f} years  ({diff * 365.25:+.2f} days)")
    print()

    base_diff_days = diff * 365.25

    # ── AGE SWEEP ─────────────────────────────────────────────────────────────
    print("AGE SWEEP (base-case parameters)")
    print("-" * 40)
    age_results = []
    for age in [25, 30, 35, 40, 45, 50, 55, 60]:
        le_a, _ = run_arm(p, n, age, voucher=True)
        le_b, _ = run_arm(p, n, age, voucher=False)
        d_days = (le_a - le_b) * 365.25
        age_results.append({"age": age, "diff_days": d_days,
                            "le_voucher": le_a, "le_control": le_b})
        print(f"  Age {age:2d}: LE voucher={le_a:.4f}yr, control={le_b:.4f}yr, "
              f"ΔLE={d_days:+.2f} days")
    print()

    # ── ESRD RISK SUBGROUPS ───────────────────────────────────────────────────
    print("SUBGROUP ANALYSIS BY ESRD RISK")
    print("-" * 40)
    subgroups = [
        ("Overall non-donor",    p["esrd_15yr_nondonor"]),
        ("Black non-donor",      p["esrd_15yr_nondonor_black"]),
        ("White non-donor",      p["esrd_15yr_nondonor_white"]),
    ]
    for label, risk in subgroups:
        le_a, _ = run_arm(p, n, age_at_entry, voucher=True,  esrd_cum_risk_15=risk)
        le_b, _ = run_arm(p, n, age_at_entry, voucher=False, esrd_cum_risk_15=risk)
        d_days = (le_a - le_b) * 365.25
        print(f"  {label:<25s}  15-yr CIF={risk*100:.3f}%  ΔLE={d_days:+.2f} days")
    print()

    # ── ONE-WAY SENSITIVITY ANALYSIS (OWSA) ───────────────────────────────────
    print("ONE-WAY SENSITIVITY ANALYSIS")
    print("-" * 40)

    # Pessimistic DDKT baseline: post-Tx mortality rates ×1.25.
    # Used as the "low" end of the Post-Tx quality row so both ends of that
    # bar deviate from the base case (avoiding a one-sided bar).
    ptx_pessimistic = {
        k: p[k] * 1.25
        for k in ("posttx_annual_mort_age1834", "posttx_annual_mort_age3549",
                  "posttx_annual_mort_age5064", "posttx_annual_mort_age65p",
                  "posttx_annual_mort")
    }

    # Each entry: (label, lo_overrides, hi_overrides, ldkt_lo, ldkt_hi)
    owsa_scenarios = [
        (
            "Voucher wait time (days)",
            {"wl_pld_median_days": 200.0}, {"wl_pld_median_days": 50.0},
            False, False,
        ),
        (
            "Standard wait time (days)",
            {"wl_std_median_days": 750.0}, {"wl_std_median_days": 1200.0},
            False, False,
        ),
        (
            "Per-cycle listing prob (from dialysis)",
            {"wl_listing_prob": p.get("wl_listing_prob_sens_low", 0.05)},
            {"wl_listing_prob": p.get("wl_listing_prob_sens_high", 0.30)},
            False, False,
        ),
        (
            "Preemptive listing prob (voucher)",
            {"esrd_preemptive_prob_pld": p["esrd_preemptive_prob_std"]},
            {"esrd_preemptive_prob_pld": 0.15},
            False, False,
        ),
        (
            "Post-Tx mortality (DDKT ×1.25 → LDKT)",
            ptx_pessimistic, {},   # low = pessimistic DDKT, high = LDKT
            False, True,
        ),
    ]

    owsa_results = []
    for label, lo_overrides, hi_overrides, ldkt_lo, ldkt_hi in owsa_scenarios:
        p_lo = {**p, **lo_overrides}
        le_v_lo, _ = run_arm(p_lo, n, age_at_entry, voucher=True,  ldkt=ldkt_lo)
        le_c_lo, _ = run_arm(p_lo, n, age_at_entry, voucher=False, ldkt=ldkt_lo)

        p_hi = {**p, **hi_overrides}
        le_v_hi, _ = run_arm(p_hi, n, age_at_entry, voucher=True,  ldkt=ldkt_hi)
        le_c_hi, _ = run_arm(p_hi, n, age_at_entry, voucher=False, ldkt=ldkt_hi)

        low_days  = (le_v_lo - le_c_lo) * 365.25
        high_days = (le_v_hi - le_c_hi) * 365.25
        owsa_results.append({"label": label, "low_days": low_days, "high_days": high_days})
        print(f"  {label}:  low={low_days:+.2f} d    high={high_days:+.2f} d")
    print()

    # ── FIGURES ───────────────────────────────────────────────────────────────
    fig_occ = make_fig_state_occupancy(trace_v, trace_c, n, age_at_entry)
    occ_path = RESULTS / "voucher_cohort_state_occupancy.png"
    fig_occ.savefig(occ_path, dpi=150, bbox_inches="tight",
                    facecolor=fig_occ.get_facecolor())
    plt.close()
    print(f"Figure saved: {occ_path}")

    fig_surv = make_fig_survival(trace_v, trace_c, n, age_at_entry)
    surv_path = RESULTS / "voucher_cohort_survival.png"
    fig_surv.savefig(surv_path, dpi=150, bbox_inches="tight",
                     facecolor=fig_surv.get_facecolor())
    plt.close()
    print(f"Figure saved: {surv_path}")

    fig_tornado = make_fig_tornado(owsa_results, base_diff_days)
    tornado_path = RESULTS / "voucher_cohort_owsa_tornado.png"
    fig_tornado.savefig(tornado_path, dpi=150, bbox_inches="tight",
                        facecolor=fig_tornado.get_facecolor())
    plt.close()
    print(f"Figure saved: {tornado_path}")

    fig_age = make_fig_age_sweep(age_results)
    age_path = RESULTS / "voucher_cohort_age_sweep.png"
    fig_age.savefig(age_path, dpi=150, bbox_inches="tight",
                    facecolor=fig_age.get_facecolor())
    plt.close()
    print(f"Figure saved: {age_path}")

    # ── SUMMARY ───────────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Base case ΔLE (voucher − control): {base_diff_days:+.2f} days"
          f"  (age {age_at_entry})")
    print()
    print("  Interpretation:")
    print(f"    A voucher-holding non-donor confers approximately")
    print(f"    {base_diff_days:.1f} days additional life expectancy vs")
    print(f"    an otherwise identical non-donor without a voucher.")
    print(f"    This benefit is entirely attributable to faster")
    print(f"    transplantation (priority waitlist) in the ~{BASE['esrd_15yr_nondonor']*100:.2f}%")
    print(f"    of voucher holders who develop ESRD over 15 years,")
    print(f"    with no countervailing elevated ESRD risk (unlike donors).")
    print()
    print("  Comparison to prior living donor analysis (07_cohort_markov.py):")
    print("    Donor ΔLE ≈ −38 days  (net harm: elevated ESRD risk dominates)")
    print(f"    Voucher ΔLE ≈ +{base_diff_days:.1f} days (net benefit: priority with no added risk)")
    print()
    print("  Largest sensitivity drivers (see tornado plot):")
    swings = [(abs(r["high_days"] - r["low_days"]), r["label"]) for r in owsa_results]
    for sw, lab in sorted(swings, reverse=True):
        print(f"    {lab}: ±{sw/2:.1f} days swing")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cohort Markov model of kidney-donation voucher life-expectancy benefit"
    )
    parser.add_argument("--age", type=int, default=40,
                        help="Age at voucher designation / cohort entry (default: 40)")
    parser.add_argument("--n", type=int, default=N_PER_ARM,
                        help="Cohort size per arm (default: 1_000_000)")
    args = parser.parse_args()
    main(age_at_entry=args.age, n=args.n)
