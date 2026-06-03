"""
07_cohort_markov.py
===================
Analytic (deterministic) cohort Markov counterpart to 06_markov_simulation.py.

Instead of random draws for each individual, propagates state *counts* forward
one year at a time and computes mean life-expectancy exactly.

Non-Markovian elements handled by state-splitting:
  D1  ESRD / dialysis, year 1 only  (elevated first-year mortality)
  D2  ESRD / dialysis, year 2+      (standard dialysis mortality)
  Return from waitlist removal and graft failure always enters D2.

States per arm:
  H    Healthy
  D1   ESRD year 1
  D2   ESRD year 2+
  WL   Waitlist  (priority for donors, standard for non-donors)
  PT   Post-transplant
  Dead Absorbing (tracked implicitly as n - alive)

Usage:
  python src/07_cohort_markov.py
  python src/07_cohort_markov.py --age 25
"""

import sys
import argparse
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load_life_table, load_params, DATA_PROC,
                   weibull_annual_prob, weibull_scale_from_cumrisk_competing,
                   median_to_annual_tx_prob)

# ── CONSTANTS ─────────────────────────────────────────────────────────────────
N_PER_ARM = 1_000_000   # cohort size per arm (result is independent of this)
MAX_AGE   = 100
CYCLE_YRS = 1.0

_lt_path = DATA_PROC / "lifetable_combined_2021.csv"
LIFE_TABLE_QX = load_life_table(_lt_path if _lt_path.exists() else None)
BASE = load_params()


# ── HELPERS ───────────────────────────────────────────────────────────────────
def _wl_tx_prob(p, priority: bool) -> float:
    median_days = p["wl_pld_median_days"] if priority else p["wl_std_median_days"]
    return float(median_to_annual_tx_prob(float(median_days)))


def _wl_mort(p) -> float:
    return float(1.0 - np.exp(-p["wl_mort_per_100py"] / 100))


def _ptx_mort(p, age: int) -> float:
    if age < 35:   return float(p.get("posttx_annual_mort_age1834", p["posttx_annual_mort"]))
    elif age < 50: return float(p.get("posttx_annual_mort_age3549", p["posttx_annual_mort"]))
    elif age < 65: return float(p.get("posttx_annual_mort_age5064", p["posttx_annual_mort"]))
    else:          return float(p.get("posttx_annual_mort_age65p",  p["posttx_annual_mort"]))


# ── COHORT SIMULATION ─────────────────────────────────────────────────────────
def run_arm(p, n: float, age_at_entry: int, donor: bool):
    """
    Propagate a cohort of n people forward through annual Markov cycles.

    Returns
    -------
    mean_le : float
        Mean remaining life-years from age_at_entry.
    state_trace : list of dicts
        Per-cycle state counts, for inspection.
    """
    # Pre-compute time-invariant transition probabilities
    wl_tx      = _wl_tx_prob(p, priority=donor)
    wl_mort    = _wl_mort(p)
    wl_remove  = float(p["wl_removal_rate_yr"])
    wl_listing = float(p.get("wl_listing_prob", 1.0))
    dial_mort1 = float(p["dialysis_1yr_mort"])
    dial_mort  = float(p["dialysis_annual_mort"])
    bg_hr      = float(p.get("donor_mort_hr", 1.0)) if donor else 1.0
    graft_fail = float(p.get("graft_annual_fail_postyear1", 0.025))

    # Weibull ESRD hazard calibrated to competing-risk 15-yr cumulative incidence
    cum_risk_15 = float(p["esrd_15yr_donor_overall"] if donor else p["esrd_15yr_nondonor"])
    wbl_k   = float(p["weibull_shape"])
    wbl_lam = weibull_scale_from_cumrisk_competing(
        cum_risk_15, wbl_k, LIFE_TABLE_QX, age_at_entry, bg_hr
    )

    # Initial state vector — everyone healthy
    H  = float(n)
    D1 = 0.0   # ESRD year 1
    D2 = 0.0   # ESRD year 2+
    WL = 0.0   # waitlist
    PT = 0.0   # post-transplant

    total_ly   = 0.0
    state_trace = []
    yr = 0

    while H + D1 + D2 + WL + PT > 0.5:
        age = age_at_entry + yr

        # ── Year-varying transition probabilities ──────────────────────────
        q_bg   = LIFE_TABLE_QX[min(age, MAX_AGE)] * bg_hr
        p_esrd = weibull_annual_prob(float(yr), wbl_lam, wbl_k)
        ptx_m  = _ptx_mort(p, age)

        # ── H (Healthy) ───────────────────────────────────────────────────
        H_die  = H * q_bg
        H_surv = H - H_die
        H_esrd = H_surv * p_esrd    # → D1 (new ESRD onset)
        H_stay = H_surv - H_esrd    # → H

        # ── D1 (ESRD year 1) ──────────────────────────────────────────────
        # After one full year, survivors leave D1 entirely: listed → WL, rest → D2
        D1_die    = D1 * dial_mort1
        D1_surv   = D1 - D1_die
        D1_listed = D1_surv * wl_listing   # → WL
        D1_to_D2  = D1_surv - D1_listed    # → D2 (survived year 1, not yet listed)

        # ── D2 (ESRD year 2+) ─────────────────────────────────────────────
        D2_die    = D2 * dial_mort
        D2_surv   = D2 - D2_die
        D2_listed = D2_surv * wl_listing   # → WL
        D2_stay   = D2_surv - D2_listed    # → D2

        # ── WL (Waitlist) ─────────────────────────────────────────────────
        WL_die    = WL * wl_mort
        WL_surv   = WL - WL_die
        WL_tx     = WL_surv * wl_tx          # → PT
        WL_after  = WL_surv - WL_tx
        WL_remove = WL_after * wl_remove     # → D2 (past first year, skips D1)
        WL_stay   = WL_after - WL_remove     # → WL

        # ── PT (Post-transplant) ──────────────────────────────────────────
        PT_die  = PT * ptx_m
        PT_surv = PT - PT_die
        PT_fail = PT_surv * graft_fail       # → D2 (graft failure, past first year)
        PT_stay = PT_surv - PT_fail          # → PT

        # ── Update state counts ───────────────────────────────────────────
        H  = H_stay
        D1 = H_esrd                                      # only fresh ESRD
        D2 = D1_to_D2 + D2_stay + WL_remove + PT_fail   # all returning flows
        WL = D1_listed + D2_listed + WL_stay
        PT = WL_tx + PT_stay

        # ── Life-years: end-of-cycle convention (matches simulation) ──────
        alive = H + D1 + D2 + WL + PT
        total_ly += alive * CYCLE_YRS

        state_trace.append({
            "yr": yr, "age": age,
            "H": H, "D1": D1, "D2": D2, "WL": WL, "PT": PT,
            "alive": alive,
        })

        yr += 1
        if yr > MAX_AGE - age_at_entry + 10:   # safety cap
            break

    return total_ly / n, state_trace


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main(age_at_entry: int = 40, n: int = N_PER_ARM):
    print("=" * 55)
    print("COHORT MARKOV MODEL (analytic, deterministic)")
    print("=" * 55)
    print(f"Age at entry : {age_at_entry}")
    print(f"n per arm    : {n:,}  (2× total)")
    print()

    p = BASE.copy()

    le_d,  trace_d  = run_arm(p, n, age_at_entry, donor=True)
    le_nd, trace_nd = run_arm(p, n, age_at_entry, donor=False)
    diff = le_d - le_nd

    print(f"  LE donor:      {le_d:.4f}  remaining years (from age {age_at_entry})")
    print(f"  LE non-donor:  {le_nd:.4f}  remaining years")
    print(f"  Difference:    {diff:+.4f} years  ({diff * 365.25:+.1f} days)")
    print()

    # Show last year both arms have meaningful occupancy, as a sanity check
    last = trace_d[-1]
    print(f"  Donor arm last cycle: yr={last['yr']}, age={last['age']}, "
          f"alive={last['alive']:.1f}")
    last = trace_nd[-1]
    print(f"  Non-donor arm last:   yr={last['yr']}, age={last['age']}, "
          f"alive={last['alive']:.1f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--age", type=int, default=40,
                        help="Age at donation / entry (default: 40)")
    parser.add_argument("--n", type=int, default=N_PER_ARM,
                        help="Cohort size per arm (default: 1_000_000)")
    args = parser.parse_args()
    main(age_at_entry=args.age, n=args.n)
