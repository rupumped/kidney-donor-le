"""
10_esrd_conditional_cohort_markov.py
======================================
Analytic cohort Markov model restricted to individuals who have already
developed ESRD. Quantifies the life-expectancy benefit of priority waitlist
access (identical to the prior-living-donor benefit) for this population,
independent of the probability of ever reaching ESRD.

Motivation: scripts 07 and 09 model populations that START healthy and may
or may not develop ESRD over their lifetime. The population-level ΔLE is
therefore diluted by the large fraction who never develop ESRD. This script
conditions entirely on ESRD onset, answering:

  "For the person who has just been diagnosed with ESRD, how much longer
   does priority waitlist access let them live compared to standard access?"

Arms:
  Priority  — priority waitlist (≈100-day median, voucher / prior-donor benefit)
  Standard  — standard waitlist (≈985-day median)

Entry state:
  All n individuals enter with newly diagnosed ESRD.
  Preemptive-listing fraction (those listed before starting dialysis)
  transitions directly to WL at cycle 0; the remainder begin in D1.

No H (Healthy) state exists. State-specific mortalities (dialysis, waitlist,
post-Tx) are all-cause rates measured in those populations and incorporate
background aging; no separate life-table term is added to disease states.
Post-transplant mortality is age-stratified as in scripts 07 and 09.

States:
  D1   ESRD / dialysis year 1   (22 %/yr all-cause mortality)
  D2   ESRD / dialysis year 2+  (17 %/yr all-cause mortality)
  WL   Waitlist — priority or standard
  PT   Post-transplant
  Dead Absorbing (tracked implicitly)

Sensitivity analyses:
  1. Priority wait time:      50 d / 102.6 d (base) / 200 d
  2. Standard wait time:      750 d / 985 d (base) / 1 200 d
  3. Per-cycle listing rate:  0.05 / 0.15 (base) / 0.30
  4. Preemptive listing prob: 5.8 % (standard) / 9.4 % (base, priority arm)
  5. Post-Tx quality:         DDKT (base) / LDKT
  6. Dialysis mortality:      −20 % / base / +20 %  (uniform scale on both yr-1 and yr-2+)

Usage:
  python src/10_esrd_conditional_cohort_markov.py
  python src/10_esrd_conditional_cohort_markov.py --age 55
  python src/10_esrd_conditional_cohort_markov.py --age 65 --n 500000
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load_params, DATA_PROC, RESULTS,
                   median_to_annual_tx_prob, load_life_table)

# ── CONSTANTS ─────────────────────────────────────────────────────────────────
N_PER_ARM = 1_000_000
MAX_AGE   = 100
CYCLE_YRS = 1.0

BASE = load_params()
_lt_path = DATA_PROC / "lifetable_combined_2021.csv"
LIFE_TABLE_QX = load_life_table(_lt_path if _lt_path.exists() else None)

# Color palette (consistent with rest of project)
_BLUE  = "#2B6CB0"   # priority
_GREY  = "#718096"   # standard
_LIGHT = "#F1EFE8"
_STATE_COLORS = {
    "D1": "#E87A5D",
    "D2": "#C4503A",
    "WL": "#534AB7",
    "PT": "#BA7517",
}


# ── HELPERS ───────────────────────────────────────────────────────────────────
def _ptx_mort(p: dict, age: int, ldkt: bool = False) -> float:
    """Age-stratified annual post-Tx all-cause mortality. DDKT base; LDKT sensitivity."""
    prefix   = "posttx_ld_" if ldkt else "posttx_"
    fallback = float(p.get("posttx_ld_annual_mort" if ldkt else "posttx_annual_mort",
                           p["posttx_annual_mort"]))
    if age < 35:   return float(p.get(f"{prefix}annual_mort_age1834", fallback))
    elif age < 50: return float(p.get(f"{prefix}annual_mort_age3549", fallback))
    elif age < 65: return float(p.get(f"{prefix}annual_mort_age5064", fallback))
    else:          return float(p.get(f"{prefix}annual_mort_age65p",  fallback))


# ── COHORT SIMULATION ─────────────────────────────────────────────────────────
def run_arm(p: dict, n: float, age_at_esrd: int, priority: bool,
            dial_mort_scale: float = 1.0, ldkt: bool = False):
    """
    Propagate a cohort that has already developed ESRD through annual cycles.

    Parameters
    ----------
    p : dict
        Model parameters.
    n : float
        Initial cohort size (result is independent of this).
    age_at_esrd : int
        Age at ESRD onset / cohort entry.
    priority : bool
        True  → priority waitlist (~100-day median).
        False → standard waitlist (~985-day median).
    dial_mort_scale : float
        Multiplicative scale on both dialysis mortality rates (sensitivity).
    ldkt : bool
        True → LDKT post-Tx mortality (sensitivity). Default: DDKT.

    Returns
    -------
    mean_le : float
        Mean remaining life-years from ESRD onset.
    state_trace : list of dict
        Per-cycle state occupancy counts.
    tx_count : float
        Total number transplanted (ever reached PT) as fraction of n.
    """
    median_days = float(p["wl_pld_median_days"] if priority else p["wl_std_median_days"])
    wl_tx      = float(median_to_annual_tx_prob(median_days))
    wl_mort    = float(1.0 - np.exp(-p["wl_mort_per_100py"] / 100))
    wl_remove  = float(p["wl_removal_rate_yr"])
    wl_listing = float(p.get("wl_listing_prob", 0.15))
    dial_mort1 = float(p["dialysis_1yr_mort"])  * dial_mort_scale
    dial_mort  = float(p["dialysis_annual_mort"]) * dial_mort_scale
    graft_fail = float(p.get("graft_annual_fail_postyear1", 0.025))

    # Preemptive listing: fraction of ESRD-onset patients listed before dialysis.
    # Priority arm uses the donor-like (informed) rate; standard arm uses the
    # general ESRD population rate.
    preemptive_p = float(p.get(
        "esrd_preemptive_prob_pld" if priority else "esrd_preemptive_prob_std", 0.0
    ))

    # ── Initial state: everyone enters with ESRD (no Healthy state) ──────────
    D1 = float(n) * (1.0 - preemptive_p)   # start dialysis
    WL = float(n) * preemptive_p            # listed before dialysis (preemptive)
    D2 = 0.0
    PT = 0.0

    total_ly   = 0.0
    tx_count   = 0.0
    state_trace = []
    yr = 0

    while D1 + D2 + WL + PT > 0.5:
        age   = age_at_esrd + yr
        ptx_m = _ptx_mort(p, age, ldkt=ldkt)

        # ── D1 (ESRD year 1) — all survivors leave D1 after one cycle ────────
        D1_die    = D1 * dial_mort1
        D1_surv   = D1 - D1_die
        D1_listed = D1_surv * wl_listing   # → WL
        D1_to_D2  = D1_surv - D1_listed    # → D2

        # ── D2 (ESRD year 2+) ─────────────────────────────────────────────────
        D2_die    = D2 * dial_mort
        D2_surv   = D2 - D2_die
        D2_listed = D2_surv * wl_listing   # → WL
        D2_stay   = D2_surv - D2_listed    # → D2

        # ── WL (Waitlist) ─────────────────────────────────────────────────────
        WL_die    = WL * wl_mort
        WL_surv   = WL - WL_die
        WL_tx     = WL_surv * wl_tx        # → PT
        WL_after  = WL_surv - WL_tx
        WL_remove = WL_after * wl_remove   # → D2 (removed, restarts dialysis past yr-1)
        WL_stay   = WL_after - WL_remove

        # ── PT (Post-transplant) ──────────────────────────────────────────────
        PT_die  = PT * ptx_m
        PT_surv = PT - PT_die
        PT_fail = PT_surv * graft_fail     # → D2 (graft failure, past first year)
        PT_stay = PT_surv - PT_fail

        # ── Update state counts ────────────────────────────────────────────────
        D1 = 0.0                                       # no new ESRD entries after cycle 0
        D2 = D1_to_D2 + D2_stay + WL_remove + PT_fail
        WL = D1_listed + D2_listed + WL_stay
        PT = WL_tx + PT_stay

        alive = D1 + D2 + WL + PT
        total_ly += alive * CYCLE_YRS
        tx_count += WL_tx

        state_trace.append({
            "yr": yr + 1, "age": age + 1,   # end-of-cycle age
            "D1": D1, "D2": D2, "WL": WL, "PT": PT,
            "alive": alive,
        })

        yr += 1
        if yr > MAX_AGE - age_at_esrd + 10:
            break

    return total_ly / n, state_trace, tx_count / n


# ── PLOTTING ──────────────────────────────────────────────────────────────────
def make_fig_state_occupancy(trace_p, trace_s, n, age_at_esrd):
    """Fraction of cohort in each state over time."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    fig.patch.set_facecolor("#FAFAF8")

    for ax, trace, title in zip(
        axes, [trace_p, trace_s],
        ["Priority waitlist (voucher / prior donor)", "Standard waitlist (control)"]
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
        "State occupancy — ESRD-conditional cohort Markov",
        fontsize=12, fontweight="bold", color="#2C2C2A"
    )
    fig.tight_layout()
    return fig


def make_fig_survival(trace_p, trace_s, n, age_at_esrd):
    """Survival curves from ESRD onset."""
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor("#FAFAF8")
    ax.set_facecolor(_LIGHT)

    # Cycle 0: full cohort alive (before any transitions)
    ages_p = [age_at_esrd] + [t["age"] for t in trace_p]
    surv_p = [1.0]         + [t["alive"] / n for t in trace_p]
    ages_s = [age_at_esrd] + [t["age"] for t in trace_s]
    surv_s = [1.0]         + [t["alive"] / n for t in trace_s]

    ax.plot(ages_p, surv_p, color=_BLUE, lw=2.0, label="Priority (voucher / prior donor)")
    ax.plot(ages_s, surv_s, color=_GREY, lw=2.0, label="Standard", linestyle="--")
    ax.legend(fontsize=10, frameon=False)
    ax.set_xlabel("Age", fontsize=10)
    ax.set_ylabel("Survival fraction", fontsize=10)
    ax.set_ylim(0, 1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#B4B2A9")
    ax.tick_params(colors="#5F5E5A", labelsize=9)
    ax.set_title("Survival from ESRD onset — priority vs standard waitlist",
                 fontsize=12, fontweight="bold", color="#2C2C2A")
    fig.tight_layout()
    return fig


def make_fig_tornado(owsa_results, base_diff_days):
    """One-way sensitivity analysis tornado plot."""
    labels = [r["label"] for r in owsa_results]
    lows   = [r["low_days"]  - base_diff_days for r in owsa_results]
    highs  = [r["high_days"] - base_diff_days for r in owsa_results]

    swings = [abs(h - l) for h, l in zip(highs, lows)]
    order  = sorted(range(len(swings)), key=lambda i: swings[i])
    labels = [labels[i] for i in order]
    lows   = [lows[i]   for i in order]
    highs  = [highs[i]  for i in order]

    fig, ax = plt.subplots(figsize=(10, max(4, len(labels) * 0.6 + 1.5)))
    fig.patch.set_facecolor("#FAFAF8")
    ax.set_facecolor(_LIGHT)

    y_pos = np.arange(len(labels))
    for i, (lo, hi) in enumerate(zip(lows, highs)):
        ax.barh(i, lo, left=0, height=0.5, color=_BLUE, alpha=0.45)
        ax.barh(i, hi, left=0, height=0.5, color=_BLUE, alpha=0.80)

    ax.axvline(0, color="#2C2C2A", lw=1.2)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Change in ΔLE from base case (days)", fontsize=10)
    ax.set_title(
        f"One-way sensitivity — priority vs standard waitlist LE benefit\n"
        f"(base ΔLE = +{base_diff_days:.1f} days, conditional on ESRD)",
        fontsize=11, fontweight="bold", color="#2C2C2A"
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#B4B2A9")
    ax.tick_params(colors="#5F5E5A", labelsize=9)
    fig.tight_layout()
    return fig


def make_fig_age_sweep(age_results):
    """ΔLE (priority benefit) vs age at ESRD onset."""
    ages  = [r["age"] for r in age_results]
    diffs = [r["diff_days"] for r in age_results]

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor("#FAFAF8")
    ax.set_facecolor(_LIGHT)

    ax.plot(ages, diffs, color=_BLUE, lw=2.2, marker="o", markersize=6)
    ax.set_ylim(bottom=0)

    for age, d in zip(ages, diffs):
        ax.annotate(f"+{d:.0f}d", xy=(age, d), xytext=(2, 6),
                    textcoords="offset points", fontsize=8, color="#2C2C2A")

    ax.set_xlabel("Age at ESRD onset", fontsize=10)
    ax.set_ylabel("LE benefit of priority waitlist (days)", fontsize=10)
    ax.set_title(
        "Priority waitlist LE benefit by age at ESRD onset\n"
        "(DDKT outcomes, base-case parameters, conditional on ESRD)",
        fontsize=11, fontweight="bold", color="#2C2C2A"
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#B4B2A9")
    ax.tick_params(colors="#5F5E5A", labelsize=9)
    fig.tight_layout()
    return fig


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main(age_at_esrd: int = 60, n: int = N_PER_ARM):
    print("=" * 62)
    print("ESRD-CONDITIONAL COHORT MARKOV (analytic, deterministic)")
    print("=" * 62)
    print(f"Age at ESRD onset : {age_at_esrd}")
    print(f"n per arm         : {n:,}  (2× total)")
    print("Entry condition   : all individuals already have ESRD")
    print()

    p = BASE.copy()

    # ── BASE CASE ─────────────────────────────────────────────────────────────
    le_p, trace_p, tx_p = run_arm(p, n, age_at_esrd, priority=True)
    le_s, trace_s, tx_s = run_arm(p, n, age_at_esrd, priority=False)
    diff = le_p - le_s

    print("BASE CASE RESULTS")
    print("-" * 45)
    print(f"  LE priority arm : {le_p:.4f}  remaining years (from ESRD onset at {age_at_esrd})")
    print(f"  LE standard arm : {le_s:.4f}  remaining years")
    print(f"  ΔLE (priority ben.): {diff:+.4f} years  ({diff * 365.25:+.1f} days)")
    print(f"  Ever transplanted — priority: {tx_p*100:.1f}%   standard: {tx_s*100:.1f}%")
    print()

    base_diff_days = diff * 365.25

    # ── AGE SWEEP ─────────────────────────────────────────────────────────────
    print("AGE SWEEP (base-case parameters)")
    print("-" * 45)
    age_results = []
    for age in [40, 45, 50, 55, 60, 65, 70, 75]:
        le_a, _, tx_a = run_arm(p, n, age, priority=True)
        le_b, _, tx_b = run_arm(p, n, age, priority=False)
        d_days = (le_a - le_b) * 365.25
        age_results.append({"age": age, "diff_days": d_days,
                            "le_priority": le_a, "le_standard": le_b,
                            "tx_priority": tx_a, "tx_standard": tx_b})
        print(f"  Age {age:2d}: LE priority={le_a:.3f}yr, standard={le_b:.3f}yr, "
              f"ΔLE={d_days:+.1f}d   Tx: {tx_a*100:.1f}% vs {tx_b*100:.1f}%")
    print()

    # ── ONE-WAY SENSITIVITY ANALYSIS ──────────────────────────────────────────
    print("ONE-WAY SENSITIVITY ANALYSIS")
    print("-" * 45)

    owsa_scenarios: list[tuple[str, dict, dict, bool, bool, float, float]] = [
        # (label, lo_overrides, hi_overrides, ldkt_lo, ldkt_hi, dial_scale_lo, dial_scale_hi)
        (
            "Priority wait time (days)",
            {"wl_pld_median_days": 200.0}, {"wl_pld_median_days": 50.0},
            False, False, 1.0, 1.0,
        ),
        (
            "Standard wait time (days)",
            {"wl_std_median_days": 750.0}, {"wl_std_median_days": 1200.0},
            False, False, 1.0, 1.0,
        ),
        (
            "Per-cycle listing prob (from dialysis)",
            {"wl_listing_prob": p.get("wl_listing_prob_sens_low", 0.05)},
            {"wl_listing_prob": p.get("wl_listing_prob_sens_high", 0.30)},
            False, False, 1.0, 1.0,
        ),
        (
            "Preemptive listing (priority arm)",
            {"esrd_preemptive_prob_pld": p["esrd_preemptive_prob_std"]},  # conservative
            {"esrd_preemptive_prob_pld": 0.15},                            # optimistic
            False, False, 1.0, 1.0,
        ),
        (
            "Post-Tx quality (DDKT vs LDKT)",
            {}, {},
            False, True, 1.0, 1.0,   # ldkt_lo=False (DDKT base), ldkt_hi=True (LDKT)
        ),
        (
            "Dialysis mortality (±20 %)",
            {}, {},
            False, False, 1.20, 0.80,  # lo = +20% (worse), hi = −20% (better)
        ),
    ]

    owsa_results = []
    for label, lo_ov, hi_ov, ldkt_lo, ldkt_hi, dscale_lo, dscale_hi in owsa_scenarios:
        p_lo = {**p, **lo_ov}
        le_p_lo, _, _ = run_arm(p_lo, n, age_at_esrd, priority=True,
                                 dial_mort_scale=dscale_lo, ldkt=ldkt_lo)
        le_s_lo, _, _ = run_arm(p_lo, n, age_at_esrd, priority=False,
                                 dial_mort_scale=dscale_lo, ldkt=ldkt_lo)

        p_hi = {**p, **hi_ov}
        le_p_hi, _, _ = run_arm(p_hi, n, age_at_esrd, priority=True,
                                 dial_mort_scale=dscale_hi, ldkt=ldkt_hi)
        le_s_hi, _, _ = run_arm(p_hi, n, age_at_esrd, priority=False,
                                 dial_mort_scale=dscale_hi, ldkt=ldkt_hi)

        low_days  = (le_p_lo - le_s_lo) * 365.25
        high_days = (le_p_hi - le_s_hi) * 365.25
        owsa_results.append({"label": label, "low_days": low_days, "high_days": high_days})

        print(f"  {label}")
        print(f"    low: {low_days:+.1f} d    high: {high_days:+.1f} d   "
              f"(base: {base_diff_days:+.1f} d)")
    print()

    # ── FIGURES ───────────────────────────────────────────────────────────────
    fig_occ = make_fig_state_occupancy(trace_p, trace_s, n, age_at_esrd)
    occ_path = RESULTS / "esrd_conditional_state_occupancy.png"
    fig_occ.savefig(occ_path, dpi=150, bbox_inches="tight",
                    facecolor=fig_occ.get_facecolor())
    plt.close()
    print(f"Figure saved: {occ_path}")

    fig_surv = make_fig_survival(trace_p, trace_s, n, age_at_esrd)
    surv_path = RESULTS / "esrd_conditional_survival.png"
    fig_surv.savefig(surv_path, dpi=150, bbox_inches="tight",
                     facecolor=fig_surv.get_facecolor())
    plt.close()
    print(f"Figure saved: {surv_path}")

    fig_tornado = make_fig_tornado(owsa_results, base_diff_days)
    tornado_path = RESULTS / "esrd_conditional_owsa_tornado.png"
    fig_tornado.savefig(tornado_path, dpi=150, bbox_inches="tight",
                        facecolor=fig_tornado.get_facecolor())
    plt.close()
    print(f"Figure saved: {tornado_path}")

    fig_age = make_fig_age_sweep(age_results)
    age_path = RESULTS / "esrd_conditional_age_sweep.png"
    fig_age.savefig(age_path, dpi=150, bbox_inches="tight",
                    facecolor=fig_age.get_facecolor())
    plt.close()
    print(f"Figure saved: {age_path}")

    # ── SUMMARY ───────────────────────────────────────────────────────────────
    print()
    print("=" * 62)
    print("SUMMARY")
    print("=" * 62)
    print(f"  Base case ΔLE (priority − standard): {base_diff_days:+.1f} days"
          f"  (ESRD onset age {age_at_esrd})")
    print()
    print("  Comparison across scripts:")
    print("    07 donor vs non-donor (pop. level):    ΔLE ≈ −38 days")
    print(f"    09 voucher vs control (pop. level):    ΔLE ≈ +0.3 days")
    print(f"    10 priority vs standard | ESRD:        ΔLE ≈ {base_diff_days:+.1f} days  ←")
    print()
    print("  Interpretation:")
    print(f"    Among people who have already developed ESRD, priority")
    print(f"    waitlist access adds approximately {base_diff_days:.0f} days of life")
    print(f"    expectancy. This is the 'per-patient-who-needs-it' benefit;")
    print(f"    it dilutes to ~+0.3 days at the population level (script 09)")
    print(f"    because only {BASE['esrd_15yr_nondonor']*100:.3f}% of non-donors develop")
    print(f"    ESRD within 15 years and most never reach the waitlist.")
    print()
    print("  Largest sensitivity drivers (see tornado plot):")
    swings = [(abs(r["high_days"] - r["low_days"]), r["label"]) for r in owsa_results]
    for sw, lab in sorted(swings, reverse=True):
        print(f"    {lab}: ±{sw/2:.1f} days swing")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="ESRD-conditional cohort Markov: priority vs standard waitlist LE"
    )
    parser.add_argument("--age", type=int, default=60,
                        help="Age at ESRD onset / cohort entry (default: 60)")
    parser.add_argument("--n", type=int, default=N_PER_ARM,
                        help="Cohort size per arm (default: 1_000_000)")
    args = parser.parse_args()
    main(age_at_esrd=args.age, n=args.n)
