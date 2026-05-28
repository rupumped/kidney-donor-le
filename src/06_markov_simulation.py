"""
06_markov_simulation.py
═══════════════════════
Living Kidney Donation & Priority Policy: Markov Monte Carlo Life-Expectancy Model

Compares lifetime outcomes for:
  ARM A: Donate a kidney (elevated ESRD risk, priority waitlist access if ESRD develops)
  ARM B: Do not donate (baseline ESRD risk, standard waitlist if ESRD develops)

Health states:
  0 = Healthy (post-donation or matched non-donor)
  1 = ESRD (on dialysis, awaiting listing)
  2 = Standard waitlist
  3 = Priority waitlist  (ARM A only)
  4 = Post-transplant
  5 = Dead (absorbing)

Parameters loaded from data/processed/params.json (produced by scripts 01–05).
Falls back to hard-coded confirmed values if params.json is absent.

Usage:
  python src/06_markov_simulation.py

Outputs in results/:
  kidney_model_results.png      — six-panel figure
  kidney_model_results.csv      — summary table
  kidney_model_race.png         — race-stratified six-panel figure
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load_life_table, load_params, RESULTS, DATA_PROC,
                   beta_params_from_mean_se, weibull_annual_hazard,
                   weibull_scale_from_cumrisk, median_to_annual_tx_prob)

# ── Global constants ──────────────────────────────────────────────────────────
N_SIM     = 100_000   # individuals per arm (base case)
N_DRAWS   = 500       # PSA parameter draws
MAX_AGE   = 100
CYCLE_YRS = 1.0       # annual Markov cycles

# ── Life table ────────────────────────────────────────────────────────────────
_lt_path = DATA_PROC / "lifetable_combined_2021.csv"
LIFE_TABLE_QX = load_life_table(_lt_path if _lt_path.exists() else None)

# ── Base-case parameters ──────────────────────────────────────────────────────
# Loaded from params.json if available; otherwise hard-coded confirmed values.
BASE = load_params()

# ── PARAMETER SAMPLING FOR PSA ──────────────────────────────────────────────
def sample_params(rng):
    """Draw one set of parameters from uncertainty distributions (PSA)."""
    p = dict(BASE)

    def beta_sample(key, se_frac=0.20):
        mean = float(BASE[key])
        a, b = beta_params_from_mean_se(mean, se_frac)
        return float(rng.beta(a, b))

    p["esrd_15yr_donor_overall"] = beta_sample("esrd_15yr_donor_overall")
    p["esrd_15yr_nondonor"]      = beta_sample("esrd_15yr_nondonor")
    p["esrd_15yr_donor_black"]   = beta_sample("esrd_15yr_donor_black")
    p["esrd_15yr_donor_white"]   = beta_sample("esrd_15yr_donor_white")

    # Weibull shape: log-normal around 1.5 (SE 0.13 on log scale ≈ 13%)
    p["weibull_shape"] = float(np.exp(rng.normal(np.log(1.5), 0.13)))

    # Wait time medians: log-normal (SE 15% standard, 20% PLD)
    p["wl_std_median_days"] = float(np.exp(
        rng.normal(np.log(BASE["wl_std_median_days"]), 0.15)))
    p["wl_pld_median_days"] = float(np.exp(
        rng.normal(np.log(BASE["wl_pld_median_days"]), 0.20)))

    # Waitlist mortality: normal (SE 0.5/100 PY)
    p["wl_mort_per_100py"] = float(np.clip(
        rng.normal(BASE["wl_mort_per_100py"], 0.5), 2.0, 12.0))

    # Post-transplant and dialysis mortality: beta distributions
    p["posttx_annual_mort"]    = beta_sample("posttx_annual_mort")
    p["dialysis_annual_mort"]  = beta_sample("dialysis_annual_mort")

    # Donor mortality HR: log-normal (mean 1.0, SD 0.10 on log scale)
    p["donor_mort_hr"] = float(np.exp(rng.normal(0.0, 0.10)))

    return p


# ── ANNUAL TRANSITION RATES ──────────────────────────────────────────────────
def esrd_annual_hazard(p, donor: bool, age_at_donation: int, age_now: int):
    """
    Weibull annual ESRD hazard at age_now for someone who donated at age_at_donation.
    Calibrated so that integrating over 15 yrs reproduces Muzaale 2014 CIF.
    Returns 0 at time of donation (t=0); accelerates thereafter (k=1.5).
    """
    cum_risk_15 = p["esrd_15yr_donor_overall"] if donor else p["esrd_15yr_nondonor"]
    k   = p["weibull_shape"]
    lam = weibull_scale_from_cumrisk(cum_risk_15, k)
    t   = max(float(age_now - age_at_donation), 0.0)
    return float(np.clip(weibull_annual_hazard(t, lam, k), 0, 0.30))


def waitlist_annual_tx_prob(p, priority: bool):
    """Annual probability of receiving a transplant, given median wait time."""
    median_days = p["wl_pld_median_days"] if priority else p["wl_std_median_days"]
    return median_to_annual_tx_prob(float(median_days))


def waitlist_annual_mort(p):
    """Annual mortality probability while on waitlist."""
    return float(p["wl_mort_per_100py"] / 100)


# ── COHORT SIMULATION (vectorised) ──────────────────────────────────────────
def simulate_cohort(p, n, age_at_entry, donor: bool, rng):
    """
    Simulate n individuals from age_at_entry to MAX_AGE.
    Returns array of life-years lived per person.
    """
    # States: 0=Healthy, 1=ESRD/dialysis, 2=WL-standard, 3=WL-priority,
    #         4=PostTx, 5=Dead
    state   = np.zeros(n, dtype=np.int8)          # start: Healthy
    age     = np.full(n, float(age_at_entry))
    ly      = np.zeros(n)                         # life-years accumulated
    alive   = np.ones(n, dtype=bool)

    # Pre-compute annual transition probs (age-independent simplification
    # for non-ESRD states; ESRD hazard is age-dependent)
    wl_tx_std  = waitlist_annual_tx_prob(p, priority=False)
    wl_tx_pld  = waitlist_annual_tx_prob(p, priority=True)
    wl_mort    = waitlist_annual_mort(p)
    wl_remove  = p["wl_removal_rate_yr"]

    dial_mort1 = p["dialysis_1yr_mort"]   # first year on dialysis
    dial_mort  = p["dialysis_annual_mort"]
    ptx_mort   = p["posttx_annual_mort"]
    bg_hr      = p["donor_mort_hr"] if donor else 1.0

    # Track time in ESRD state for first-year mortality
    esrd_time  = np.zeros(n)

    for yr in range(int(MAX_AGE - age_at_entry)):
        a = int(age_at_entry) + yr
        if a >= MAX_AGE:
            break

        # Background mortality this cycle (all states)
        q_bg = LIFE_TABLE_QX[min(a, MAX_AGE)] * bg_hr

        u = rng.random((n, 6))  # draws for each possible event

        new_state = state.copy()

        # ── Healthy (0) ──────────────────────────────────────────────────
        mask0 = (state == 0) & alive
        if mask0.any():
            p_esrd = esrd_annual_hazard(p, donor, age_at_entry, a)
            # Die from background causes
            die_bg = u[mask0, 0] < q_bg
            # Develop ESRD (conditional on not dying)
            get_esrd = (~die_bg) & (u[mask0, 1] < p_esrd)

            idx = np.where(mask0)[0]
            new_state[idx[die_bg]] = 5
            new_state[idx[~die_bg & get_esrd]] = 1
            esrd_time[idx[~die_bg & get_esrd]] = 0

        # ── ESRD / Dialysis (1) ──────────────────────────────────────────
        mask1 = (state == 1) & alive
        if mask1.any():
            idx = np.where(mask1)[0]
            # Higher mortality in first year
            dm = np.where(esrd_time[idx] < 1, dial_mort1, dial_mort)
            die_dial = u[mask1, 0] < dm
            # Move to waitlist if survive
            go_wl = ~die_dial
            new_state[idx[die_dial]] = 5
            wl_state = 3 if donor else 2
            new_state[idx[go_wl]] = wl_state
            esrd_time[mask1] += 1

        # ── Standard Waitlist (2) ────────────────────────────────────────
        mask2 = (state == 2) & alive
        if mask2.any():
            idx = np.where(mask2)[0]
            die_wl      = u[mask2, 0] < wl_mort
            get_tx      = (~die_wl) & (u[mask2, 1] < wl_tx_std)
            get_removed = (~die_wl) & (~get_tx) & (u[mask2, 2] < wl_remove)
            # Removed → back to dialysis (ESRD state)
            new_state[idx[die_wl]]      = 5
            new_state[idx[get_tx]]      = 4
            new_state[idx[get_removed]] = 1
            esrd_time[idx[get_removed]] = 1  # already past first year

        # ── Priority Waitlist (3) ─────────────────────────────────────────
        mask3 = (state == 3) & alive
        if mask3.any():
            idx = np.where(mask3)[0]
            die_wl      = u[mask3, 0] < wl_mort
            get_tx      = (~die_wl) & (u[mask3, 1] < wl_tx_pld)
            get_removed = (~die_wl) & (~get_tx) & (u[mask3, 2] < wl_remove)
            new_state[idx[die_wl]]      = 5
            new_state[idx[get_tx]]      = 4
            new_state[idx[get_removed]] = 1
            esrd_time[idx[get_removed]] = 1

        # ── Post-transplant (4) ───────────────────────────────────────────
        mask4 = (state == 4) & alive
        if mask4.any():
            idx = np.where(mask4)[0]
            die_ptx = u[mask4, 0] < ptx_mort
            # Graft failure (approximate): 1-yr graft failure ~2-3% after year 1
            graft_fail = (~die_ptx) & (u[mask4, 1] < 0.025)
            new_state[idx[die_ptx]]    = 5
            new_state[idx[graft_fail]] = 1  # back to dialysis
            esrd_time[idx[graft_fail]] = 1

        # ── Accumulate life-years for surviving individuals ───────────────
        survived_cycle = (new_state != 5)
        ly[survived_cycle & alive] += CYCLE_YRS

        alive = survived_cycle
        state = new_state
        age  += CYCLE_YRS

    return ly


# ── RUN BASE CASE ────────────────────────────────────────────────────────────
def run_base_case(age_at_donation=40, n=N_SIM):
    print(f"\nRunning base case — age at donation {age_at_donation}, n={n:,}")
    rng = np.random.default_rng(99)
    p = BASE.copy()

    ly_donor    = simulate_cohort(p, n, age_at_donation, donor=True,  rng=rng)
    ly_nondonor = simulate_cohort(p, n, age_at_donation, donor=False, rng=rng)

    le_donor    = ly_donor.mean()
    le_nondonor = ly_nondonor.mean()
    diff        = le_donor - le_nondonor

    print(f"  LE donor:      {le_donor:.2f} remaining years (from age {age_at_donation})")
    print(f"  LE non-donor:  {le_nondonor:.2f} remaining years")
    print(f"  Difference:    {diff:+.3f} years ({diff*365.25:+.1f} days)")

    return {
        "le_donor": le_donor,
        "le_nondonor": le_nondonor,
        "diff": diff,
        "ly_donor": ly_donor,
        "ly_nondonor": ly_nondonor,
        "age": age_at_donation,
    }


# ── PROBABILISTIC SENSITIVITY ANALYSIS ──────────────────────────────────────
def run_psa(age_at_donation=40, n_draws=N_DRAWS, n_per_draw=5000):
    print(f"\nRunning PSA — {n_draws} parameter draws, {n_per_draw:,}/draw...")
    rng = np.random.default_rng(77)
    diffs = []

    for i in range(n_draws):
        if (i+1) % 100 == 0:
            print(f"  PSA draw {i+1}/{n_draws}")
        p    = sample_params(rng)
        rng_d  = np.random.default_rng(rng.integers(1e9))
        rng_nd = np.random.default_rng(rng.integers(1e9))
        ld  = simulate_cohort(p, n_per_draw, age_at_donation, donor=True,  rng=rng_d)
        lnd = simulate_cohort(p, n_per_draw, age_at_donation, donor=False, rng=rng_nd)
        diffs.append(ld.mean() - lnd.mean())

    diffs = np.array(diffs)
    print(f"  PSA ΔLE: mean={diffs.mean():+.3f}, "
          f"95% CrI [{np.percentile(diffs,2.5):+.3f}, {np.percentile(diffs,97.5):+.3f}]")
    print(f"  P(donation beneficial): {(diffs > 0).mean():.1%}")
    return diffs


# ── ONE-WAY SENSITIVITY ANALYSIS ────────────────────────────────────────────
def run_owsa(age_at_donation=40, n=20_000):
    print("\nRunning one-way sensitivity analysis...")

    # Non-donor 15-yr risk ~3.9/10,000 = 0.00039 (Muzaale matched controls)
    # Donor base 30.8/10,000 = 0.0031; RR ~8x
    # For RR scenarios, scale donor risk proportionally keeping non-donor fixed
    scenarios = {
        "ESRD RR ×2 vs controls":        {"esrd_15yr_donor_overall": 0.00039 * 2},
        "ESRD RR ×4 vs controls":        {"esrd_15yr_donor_overall": 0.00039 * 4},
        "ESRD RR ×8 — base case":        {},
        "ESRD RR ×11 (Mjøen upper)":     {"esrd_15yr_donor_overall": 0.00039 * 11},
        "Std wait 24 months":            {"wl_std_median_days": 730},
        "Std wait 58 months (pre-KAS)":  {"wl_std_median_days": 1760},
        "PLD wait 50 days (optimistic)": {"wl_pld_median_days": 50},
        "PLD wait 200 days":             {"wl_pld_median_days": 200},
        "No priority (PLD=standard)":    {"wl_pld_median_days": 985},
        "Dialysis mort +50%":            {"dialysis_annual_mort": 0.255},
        "Dialysis mort -50%":            {"dialysis_annual_mort": 0.085},
        "Donor mort HR=1.30 (Mjøen)":    {"donor_mort_hr": 1.30},
    }

    results = {}
    seed = 55
    for label, overrides in scenarios.items():
        p = BASE.copy()
        p.update(overrides)
        # Independent RNGs per arm so stochastic noise doesn't cancel
        ld  = simulate_cohort(p, n, age_at_donation, donor=True,
                              rng=np.random.default_rng(seed))
        lnd = simulate_cohort(p, n, age_at_donation, donor=False,
                              rng=np.random.default_rng(seed + 1000))
        diff = ld.mean() - lnd.mean()
        results[label] = diff
        print(f"  {label:<44} ΔLE = {diff:+.3f} yr  ({diff*365.25:+.1f} days)")
        seed += 1

    return results


# ── AGE SUBGROUP ANALYSIS ────────────────────────────────────────────────────
def run_age_subgroups(n=20_000):
    print("\nRunning age-at-donation subgroup analysis...")
    ages = [25, 35, 40, 45, 55]
    results = {}
    for age in ages:
        rng = np.random.default_rng(age * 10)
        p = BASE.copy()
        ld  = simulate_cohort(p, n, age, donor=True,  rng=rng)
        lnd = simulate_cohort(p, n, age, donor=False, rng=rng)
        diff = ld.mean() - lnd.mean()
        results[age] = {"le_donor": ld.mean(), "le_nondonor": lnd.mean(), "diff": diff}
        print(f"  Age {age}: ΔLE = {diff:+.3f} yr  "
              f"(donor {ld.mean():.1f} yr, non-donor {lnd.mean():.1f} yr)")
    return results


# ── RACE SUBGROUP ANALYSIS ───────────────────────────────────────────────────
def run_race_subgroups(age_at_donation=40, n=20_000):
    print("\nRunning race subgroup analysis...")
    scenarios = {
        "White donor":    {"esrd_15yr_donor_overall": 0.00227, "esrd_15yr_nondonor": 0.00010},
        "Overall (base)": {},
        "Black donor":    {"esrd_15yr_donor_overall": 0.00747, "esrd_15yr_nondonor": 0.00039},
    }
    results = {}
    for label, overrides in scenarios.items():
        p = BASE.copy()
        p.update(overrides)
        rng = np.random.default_rng(42)
        ld  = simulate_cohort(p, n, age_at_donation, donor=True,  rng=rng)
        lnd = simulate_cohort(p, n, age_at_donation, donor=False, rng=rng)
        diff = ld.mean() - lnd.mean()
        results[label] = {"le_donor": ld.mean(), "le_nondonor": lnd.mean(), "diff": diff}
        print(f"  {label}: ΔLE = {diff:+.3f} yr")
    return results


# ── PLOTTING ─────────────────────────────────────────────────────────────────
def make_figure(base_res, psa_diffs, owsa_res, age_res, race_res):
    fig = plt.figure(figsize=(16, 14))
    fig.patch.set_facecolor("#FAFAF8")
    gs = gridspec.GridSpec(3, 3, figure=fig,
                           hspace=0.45, wspace=0.38,
                           left=0.07, right=0.97,
                           top=0.93, bottom=0.06)

    TEAL   = "#1D9E75"
    CORAL  = "#D85A30"
    PURPLE = "#534AB7"
    AMBER  = "#BA7517"
    GRAY   = "#888780"
    LIGHT  = "#F1EFE8"

    def style_ax(ax, title):
        ax.set_facecolor(LIGHT)
        ax.spines[["top","right"]].set_visible(False)
        ax.spines[["left","bottom"]].set_color("#B4B2A9")
        ax.tick_params(colors="#5F5E5A", labelsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold",
                     color="#2C2C2A", pad=7)

    age = base_res["age"]

    # ── (A) LE distributions base case ──────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, :2])
    bins = np.arange(0, MAX_AGE - age + 2, 2)
    ax_a.hist(base_res["ly_nondonor"], bins=bins, density=True,
              color=TEAL,  alpha=0.55, label=f"Non-donor  μ={base_res['le_nondonor']:.1f} yr")
    ax_a.hist(base_res["ly_donor"],    bins=bins, density=True,
              color=CORAL, alpha=0.55, label=f"Donor      μ={base_res['le_donor']:.1f} yr")
    ax_a.axvline(base_res["le_nondonor"], color=TEAL,  lw=1.8, ls="--")
    ax_a.axvline(base_res["le_donor"],    color=CORAL, lw=1.8, ls="--")
    ax_a.legend(fontsize=9, frameon=False)
    ax_a.set_xlabel("Remaining life-years from donation age", fontsize=9)
    ax_a.set_ylabel("Density", fontsize=9)
    style_ax(ax_a, f"(A) Lifetime distributions — base case (age {age}, n=100k/arm)")

    # ── (B) PSA scatter ──────────────────────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 2])
    ax_b.hist(psa_diffs * 365.25, bins=40, color=PURPLE, alpha=0.75,
              edgecolor="white", linewidth=0.3)
    ax_b.axvline(0, color=GRAY, lw=1.5, ls="--")
    lo, hi = np.percentile(psa_diffs * 365.25, [2.5, 97.5])
    ax_b.axvline(lo, color=AMBER, lw=1, ls=":")
    ax_b.axvline(hi, color=AMBER, lw=1, ls=":")
    p_ben = (psa_diffs > 0).mean()
    ax_b.set_xlabel("ΔLE (days): donor − non-donor", fontsize=9)
    ax_b.set_ylabel("Count", fontsize=9)
    ax_b.set_title(f"(B) PSA ({N_DRAWS} draws)\nP(beneficial)={p_ben:.0%}, "
                   f"95% CrI [{lo:.0f}, {hi:.0f}] d",
                   fontsize=10, fontweight="bold", color="#2C2C2A", pad=7)
    ax_b.set_facecolor(LIGHT)
    ax_b.spines[["top","right"]].set_visible(False)
    ax_b.spines[["left","bottom"]].set_color("#B4B2A9")
    ax_b.tick_params(colors="#5F5E5A", labelsize=9)

    # ── (C) Tornado / OWSA ───────────────────────────────────────────────────
    ax_c = fig.add_subplot(gs[1, :])
    labels = list(owsa_res.keys())
    vals   = [owsa_res[l] * 365.25 for l in labels]   # convert to days
    colors = [TEAL if v >= 0 else CORAL for v in vals]
    sorted_idx = np.argsort(np.abs(vals))
    labels_s   = [labels[i] for i in sorted_idx]
    vals_s     = [vals[i]   for i in sorted_idx]
    colors_s   = [colors[i] for i in sorted_idx]
    y_pos = np.arange(len(labels_s))
    ax_c.barh(y_pos, vals_s, color=colors_s, alpha=0.8, height=0.65)
    ax_c.axvline(0, color=GRAY, lw=1.2)
    ax_c.set_yticks(y_pos)
    ax_c.set_yticklabels(labels_s, fontsize=9)
    ax_c.set_xlabel("ΔLE (days): donor − non-donor", fontsize=9)
    style_ax(ax_c, "(C) One-way sensitivity analysis — tornado chart")

    # ── (D) Age subgroup ─────────────────────────────────────────────────────
    ax_d = fig.add_subplot(gs[2, 0])
    ages_s = sorted(age_res.keys())
    diffs_d = [age_res[a]["diff"] * 365.25 for a in ages_s]
    bar_cols = [TEAL if d >= 0 else CORAL for d in diffs_d]
    ax_d.bar(range(len(ages_s)), diffs_d, color=bar_cols, alpha=0.8, width=0.6)
    ax_d.axhline(0, color=GRAY, lw=1.2)
    ax_d.set_xticks(range(len(ages_s)))
    ax_d.set_xticklabels([f"Age {a}" for a in ages_s], fontsize=9)
    ax_d.set_ylabel("ΔLE (days)", fontsize=9)
    style_ax(ax_d, "(D) ΔLE by age at donation")

    # ── (E) Race subgroup ─────────────────────────────────────────────────────
    ax_e = fig.add_subplot(gs[2, 1])
    race_labels = list(race_res.keys())
    race_diffs  = [race_res[r]["diff"] * 365.25 for r in race_labels]
    rc = [TEAL if d >= 0 else CORAL for d in race_diffs]
    ax_e.bar(range(len(race_labels)), race_diffs, color=rc, alpha=0.8, width=0.5)
    ax_e.axhline(0, color=GRAY, lw=1.2)
    ax_e.set_xticks(range(len(race_labels)))
    ax_e.set_xticklabels(race_labels, fontsize=8.5)
    ax_e.set_ylabel("ΔLE (days)", fontsize=9)
    style_ax(ax_e, "(E) ΔLE by race")

    # ── (F) Break-even: at what wait-time gap does priority cease to matter? ─
    ax_f = fig.add_subplot(gs[2, 2])
    pld_waits = np.linspace(50, 985, 60)  # PLD wait from 50 to 985 days
    be_diffs  = []
    for pld_w in pld_waits:
        p = BASE.copy()
        p["wl_pld_median_days"] = float(pld_w)
        rng2 = np.random.default_rng(123)
        ld  = simulate_cohort(p, 10_000, 40, donor=True,  rng=rng2)
        lnd = simulate_cohort(p, 10_000, 40, donor=False, rng=rng2)
        be_diffs.append((ld.mean() - lnd.mean()) * 365.25)
    ax_f.plot(pld_waits, be_diffs, color=PURPLE, lw=2)
    ax_f.axhline(0, color=GRAY, lw=1.2, ls="--")
    ax_f.fill_between(pld_waits, be_diffs, 0,
                      where=[d > 0 for d in be_diffs],
                      alpha=0.18, color=TEAL, label="Donor beneficial")
    ax_f.fill_between(pld_waits, be_diffs, 0,
                      where=[d <= 0 for d in be_diffs],
                      alpha=0.18, color=CORAL, label="Donor net harm")
    ax_f.set_xlabel("PLD median wait (days)", fontsize=9)
    ax_f.set_ylabel("ΔLE (days)", fontsize=9)
    ax_f.legend(fontsize=8, frameon=False)
    style_ax(ax_f, "(F) Break-even: PLD wait vs ΔLE (age 40)")

    plt.suptitle(
        "Living Kidney Donation & Priority Policy: Markov Model Results\n"
        "Donor vs Matched Non-Donor Life Expectancy",
        fontsize=13, fontweight="bold", color="#2C2C2A", y=0.98)

    return fig


# ── RESULTS TABLE ─────────────────────────────────────────────────────────────
def make_results_table(base_res, psa_diffs, owsa_res, age_res, race_res):
    rows = []

    # Base case
    rows.append({
        "Analysis": "Base case (age 40, overall)",
        "LE Donor (yr)": f"{base_res['le_donor']:.2f}",
        "LE Non-donor (yr)": f"{base_res['le_nondonor']:.2f}",
        "ΔLE (days)": f"{base_res['diff']*365.25:+.1f}",
        "Note": "n=100k/arm"
    })

    # PSA summary
    rows.append({
        "Analysis": "PSA (probabilistic)",
        "LE Donor (yr)": "—",
        "LE Non-donor (yr)": "—",
        "ΔLE (days)": (f"{psa_diffs.mean()*365.25:+.1f} "
                       f"[{np.percentile(psa_diffs*365.25,2.5):+.1f}, "
                       f"{np.percentile(psa_diffs*365.25,97.5):+.1f}]"),
        "Note": f"P(beneficial)={(psa_diffs>0).mean():.0%}"
    })

    # Age subgroups
    for a, v in sorted(age_res.items()):
        rows.append({
            "Analysis": f"Age at donation: {a}",
            "LE Donor (yr)": f"{v['le_donor']:.2f}",
            "LE Non-donor (yr)": f"{v['le_nondonor']:.2f}",
            "ΔLE (days)": f"{v['diff']*365.25:+.1f}",
            "Note": ""
        })

    # Race subgroups
    for r, v in race_res.items():
        rows.append({
            "Analysis": f"Race: {r}",
            "LE Donor (yr)": f"{v['le_donor']:.2f}",
            "LE Non-donor (yr)": f"{v['le_nondonor']:.2f}",
            "ΔLE (days)": f"{v['diff']*365.25:+.1f}",
            "Note": "age 40"
        })

    # Key OWSA
    for lbl in ["No priority (PLD=standard)", "Donor mort HR=1.30 (Mjøen)",
                "ESRD RR ×11 (Mjøen upper)", "PLD wait 50 days (optimistic)"]:
        rows.append({
            "Analysis": f"Sensitivity: {lbl}",
            "LE Donor (yr)": "—",
            "LE Non-donor (yr)": "—",
            "ΔLE (days)": f"{owsa_res[lbl]*365.25:+.1f}",
            "Note": "OWSA"
        })

    df = pd.DataFrame(rows)
    return df


# ── MAIN ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("="*60)
    print("LIVING KIDNEY DONATION MARKOV MODEL")
    print("="*60)

    base_res  = run_base_case(age_at_donation=40, n=N_SIM)
    psa_diffs = run_psa(age_at_donation=40, n_draws=N_DRAWS, n_per_draw=5000)
    owsa_res  = run_owsa(age_at_donation=40, n=20_000)
    age_res   = run_age_subgroups(n=20_000)
    race_res  = run_race_subgroups(age_at_donation=40, n=20_000)

    # Save figure
    fig = make_figure(base_res, psa_diffs, owsa_res, age_res, race_res)
    fig_path = RESULTS / "kidney_model_results.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"\nFigure saved: {fig_path}")

    # Save results table
    df = make_results_table(base_res, psa_diffs, owsa_res, age_res, race_res)
    csv_path = RESULTS / "kidney_model_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"Results table saved: {csv_path}")
    print("\n" + df.to_string(index=False))
    print("\nDone.")
