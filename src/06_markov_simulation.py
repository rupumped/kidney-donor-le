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
  kidney_model_race.png             — race-stratified LE distributions
  kidney_model_sex.png              — sex-stratified LE distributions
  kidney_model_age_sex_matrix.png   — age × sex ΔLE heatmap
  kidney_model_sex_race_matrix.png  — sex × race ΔLE heatmap
  kidney_model_age_race_sex_matrix.png — age × race faceted by sex (all three variables)
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning)

# Allow running from repo root or src/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (load_life_table, load_params, RESULTS, DATA_PROC,
				   beta_params_from_mean_se, weibull_annual_prob,
				   weibull_scale_from_cumrisk, weibull_scale_from_cumrisk_competing,
				   median_to_annual_tx_prob)

# ── GLOBAL CONSTANTS ──────────────────────────────────────────────────────────
N_SIM     = 5_000_000   # individuals per arm (base case)
N_DRAWS   =       200   # PSA parameter draws
MAX_AGE   =       100
CYCLE_YRS =       1.0   # annual Markov cycles

# ── LIFE TABLE ────────────────────────────────────────────────────────────────
_lt_path = DATA_PROC / "lifetable_combined_2021.csv"
LIFE_TABLE_QX = load_life_table(_lt_path if _lt_path.exists() else None)

# ── BASE-CASE PARAMETERS ──────────────────────────────────────────────────────
# Loaded from params.json if available; otherwise hard-coded confirmed values.
BASE = load_params()

# ── PARAMETER SAMPLING FOR PSA ───────────────────────────────────────────────
def sample_params(rng):
	"""Draw one set of parameters from uncertainty distributions (PSA).

	Two categories of uncertainty are distinguished:
	  - Literature-calibrated: se_frac derived from published 95% CIs via
	    SE = (CI_upper - CI_lower) / (2 * 1.96) / mean.
	  - Scenario-based: large registry sources have negligible sampling SE;
	    the stated spread represents plausible future variation, not estimation
	    error.  These should be interpreted as scenario ranges, not credible
	    intervals.
	"""
	p = dict(BASE)

	def beta_sample(key, se_frac):
		mean = float(BASE[key])
		a, b = beta_params_from_mean_se(mean, se_frac)
		return float(rng.beta(a, b))

	# ESRD 15-yr risks: se_frac from Muzaale 2014 (JAMA 311:579) Table 2 CIs
	# SE = (CI_upper - CI_lower) / (2 * 1.96 * mean)
	p["esrd_15yr_donor_overall"] = beta_sample("esrd_15yr_donor_overall", 0.12)  # CI 24.3–38.5/10k → SE 12%
	p["esrd_15yr_nondonor"]      = beta_sample("esrd_15yr_nondonor",      0.53)  # CI 0.8–8.9/10k → SE 53% (~8 events in 20k controls)
	p["esrd_15yr_donor_black"]   = beta_sample("esrd_15yr_donor_black",   0.20)  # CI 47.8–105.8/10k → SE 20%
	p["esrd_15yr_donor_white"]   = beta_sample("esrd_15yr_donor_white",   0.16)  # CI 15.6–30.1/10k → SE 16%

	# Weibull shape: no published CI; log-normal σ=0.13 is an assumed scenario
	# range (Massie 2017 gives point estimates only; k is calibrated, not fitted)
	p["weibull_shape"] = float(np.exp(rng.normal(np.log(1.5), 0.13)))

	# Registry-derived parameters (SRTR, USRDS) have negligible sampling SE and
	# are left at their base-case values.  Scenario variation is handled by the
	# OWSA, not the PSA.

	# wl_listing_prob: back-calculated from USRDS 2025 ADR; explicitly high
	# uncertainty (OWSA range 0.05–0.30, base 0.15).  se_frac=0.40 gives a PSA
	# 95% CrI of roughly [0.03, 0.27], consistent with the acknowledged range.
	p["wl_listing_prob"] = beta_sample("wl_listing_prob", 0.40)

	# Donor mortality HR: log-normal calibrated to Grams 2018 meta-analysis
	# (PMID 29379948; Ann Intern Med 2018;168:276–284; 52 studies, 118,426 donors).
	# All-cause mortality pooled HR ≈ 0.984 (95% CI 0.743–1.302) vs nondonors.
	# mu  = [ln(0.743) + ln(1.302)] / 2 = -0.016
	# sigma = [ln(1.302) - ln(0.743)] / (2 × 1.96) = 0.143
	p["donor_mort_hr"] = float(np.exp(rng.normal(-0.016, 0.143)))

	return p


# ── ANNUAL TRANSITION RATES ───────────────────────────────────────────────────
def waitlist_annual_tx_prob(p, priority: bool):
	"""Annual probability of receiving a transplant, given median wait time."""
	median_days = p["wl_pld_median_days"] if priority else p["wl_std_median_days"]
	return median_to_annual_tx_prob(float(median_days))


def waitlist_annual_mort(p):
	"""Annual mortality probability while on waitlist (rate per 100 PY → probability)."""
	rate = float(p["wl_mort_per_100py"] / 100)
	return float(1.0 - np.exp(-rate))


# ── COHORT SIMULATION (VECTORISED) ───────────────────────────────────────────
def simulate_cohort(p, n, age_at_entry, donor: bool, rng, life_table=None):
	"""
	Simulate a cohort of n individuals through annual Markov cycles from
	age_at_entry to MAX_AGE and return remaining life-years per person.

	Parameters
	----------
	p : dict
		Parameter set mapping string keys to scalar values.  Consumed keys
		include ESRD 15-year cumulative risks, Weibull shape, waitlist median
		wait times, waitlist/dialysis/post-transplant mortality rates, graft
		failure rate, listing probability, and waitlist removal rate.  Produced
		by BASE.copy() (base case) or sample_params(rng) (PSA draw).
	n : int
		Number of individuals in the cohort.
	age_at_entry : int or float
		Age in years at which every individual enters the simulation in the
		Healthy state (state 0).
	donor : bool
		True  → donor arm: uses donor ESRD risk, applies donor_mort_hr to
				background mortality, and routes ESRD survivors to the
				priority waitlist (state 3).
		False → non-donor arm: uses non-donor ESRD risk, no HR adjustment,
				and routes ESRD survivors to the standard waitlist (state 2).
	rng : numpy.random.Generator
		Random number generator (e.g. numpy.random.default_rng(seed)).
		Passed in so callers control reproducibility across arms and PSA draws.

	Returns
	-------
	ly : numpy.ndarray, shape (n,)
		Remaining life-years lived by each individual from age_at_entry until
		death or MAX_AGE, accumulating one year per survived annual cycle.
	"""
	lt = life_table if life_table is not None else LIFE_TABLE_QX

	# States: 0=Healthy, 1=ESRD/dialysis, 2=WL-standard, 3=WL-priority,
	#         4=PostTx, 5=Dead
	state   = np.zeros(n, dtype=np.int8)          # start: Healthy
	age     = np.full(n, float(age_at_entry))
	ly      = np.zeros(n)                         # life-years accumulated
	alive   = np.ones(n, dtype=bool)

	# Pre-compute competing-risk-calibrated Weibull scale for this cohort
	cum_risk_15 = p["esrd_15yr_donor_overall"] if donor else p["esrd_15yr_nondonor"]
	wbl_k   = p["weibull_shape"]
	wbl_lam = weibull_scale_from_cumrisk_competing(
		cum_risk_15, wbl_k, lt, int(age_at_entry),
		p.get("donor_mort_hr", 1.0) if donor else 1.0
	)

	# Pre-compute age-independent transition probabilities
	wl_tx_std  = waitlist_annual_tx_prob(p, priority=False)
	wl_tx_pld  = waitlist_annual_tx_prob(p, priority=True)
	wl_mort    = waitlist_annual_mort(p)
	wl_remove  = p["wl_removal_rate_yr"]
	wl_listing = p.get("wl_listing_prob", 1.0)
	wl_state   = 3 if donor else 2

	dial_mort1 = p["dialysis_1yr_mort"]   # first year on dialysis
	dial_mort  = p["dialysis_annual_mort"]
	bg_hr      = p.get("donor_mort_hr", 1.0) if donor else 1.0
	graft_fail_rate = p.get("graft_annual_fail_postyear1", 0.025)

	# One-time preemptive listing probability at ESRD onset (USRDS 2025 Fig 7.13)
	preemptive_p = float(p.get(
		"esrd_preemptive_prob_pld" if donor else "esrd_preemptive_prob_std", 0.0))

	# Track time in ESRD state for first-year mortality
	esrd_time  = np.zeros(n)

	# Simulate one annual cycle per iteration; +1 ensures the age-100 cycle runs
	for yr in range(int(MAX_AGE - age_at_entry) + 1):
		a = int(age_at_entry) + yr

		# Background mortality this cycle (all states)
		q_bg = lt[min(a, MAX_AGE)] * bg_hr

		u = rng.random((n, 6))  # draws for each possible event

		new_state = state.copy()

		# ── (0) HEALTHY ───────────────────────────────────────────────────
		mask0 = (state == 0) & alive
		if mask0.any():
			t = max(float(a - age_at_entry), 0.0)
			p_esrd = weibull_annual_prob(t, wbl_lam, wbl_k)
			die_bg   = u[mask0, 0] < q_bg
			get_esrd = (~die_bg) & (u[mask0, 1] < p_esrd)

			idx = np.where(mask0)[0]
			new_state[idx[die_bg]] = 5
			esrd_idx = idx[get_esrd]
			if esrd_idx.size:
				# Branch at ESRD onset: preemptive listing (skip dialysis) vs dialysis yr 1
				is_preemptive = u[esrd_idx, 2] < preemptive_p
				new_state[esrd_idx[is_preemptive]]  = wl_state  # → waitlist, bypass dialysis
				new_state[esrd_idx[~is_preemptive]] = 1          # → dialysis year 1
				esrd_time[esrd_idx[~is_preemptive]] = 0

		# ── (1) ESRD / DIALYSIS ───────────────────────────────────────────
		mask1 = (state == 1) & alive
		if mask1.any():
			idx = np.where(mask1)[0]
			# Higher mortality in first year on dialysis
			dm       = np.where(esrd_time[idx] < 1, dial_mort1, dial_mort)
			die_dial = u[mask1, 0] < dm
			# Listing gate: not all survivors are listed each cycle
			get_listed = (~die_dial) & (u[mask1, 1] < wl_listing)
			new_state[idx[die_dial]]   = 5
			new_state[idx[get_listed]] = wl_state
			# Increment counter for everyone in state 1 this cycle;
			# those transitioning to WL won't re-enter mask1 next cycle
			esrd_time[mask1] += 1

		# ── (2) STANDARD WAITLIST ─────────────────────────────────────────
		mask2 = (state == 2) & alive
		if mask2.any():
			idx = np.where(mask2)[0]
			die_wl      = u[mask2, 0] < wl_mort
			get_tx      = (~die_wl) & (u[mask2, 1] < wl_tx_std)
			get_removed = (~die_wl) & (~get_tx) & (u[mask2, 2] < wl_remove)
			new_state[idx[die_wl]]      = 5
			new_state[idx[get_tx]]      = 4
			# Return to dialysis preserving accumulated esrd_time (past first year)
			new_state[idx[get_removed]] = 1

		# ── (3) PRIORITY WAITLIST ─────────────────────────────────────────
		mask3 = (state == 3) & alive
		if mask3.any():
			idx = np.where(mask3)[0]
			die_wl      = u[mask3, 0] < wl_mort
			get_tx      = (~die_wl) & (u[mask3, 1] < wl_tx_pld)
			get_removed = (~die_wl) & (~get_tx) & (u[mask3, 2] < wl_remove)
			new_state[idx[die_wl]]      = 5
			new_state[idx[get_tx]]      = 4
			new_state[idx[get_removed]] = 1

		# ── (4) POST-TRANSPLANT ───────────────────────────────────────────
		mask4 = (state == 4) & alive
		if mask4.any():
			idx = np.where(mask4)[0]
			# Age-stratified post-transplant mortality (SRTR 2023 DDKT)
			if a < 35:   ptx_mort = p.get("posttx_annual_mort_age1834", p["posttx_annual_mort"])
			elif a < 50: ptx_mort = p.get("posttx_annual_mort_age3549", p["posttx_annual_mort"])
			elif a < 65: ptx_mort = p.get("posttx_annual_mort_age5064", p["posttx_annual_mort"])
			else:        ptx_mort = p.get("posttx_annual_mort_age65p",  p["posttx_annual_mort"])
			die_ptx    = u[mask4, 0] < ptx_mort
			graft_fail = (~die_ptx) & (u[mask4, 1] < graft_fail_rate)
			new_state[idx[die_ptx]]    = 5
			new_state[idx[graft_fail]] = 1  # back to dialysis
			esrd_time[idx[graft_fail]] = 1  # already past first year

		# ── ACCUMULATE LIFE-YEARS ─────────────────────────────────────────
		survived_cycle = (new_state != 5)
		ly[survived_cycle & alive] += CYCLE_YRS

		alive = survived_cycle
		state = new_state
		age  += CYCLE_YRS

	return ly


# ── ANALYTIC COHORT MARKOV ────────────────────────────────────────────────────
def run_arm_analytic(p, age_at_entry: int, donor: bool, life_table=None) -> float:
	"""
	Deterministic cohort Markov — returns mean remaining life-years per person.

	Non-Markovian first-year dialysis mortality is handled by splitting ESRD into:
	  D1: year 1 only (dial_mort1) — only new-onset ESRD enters here.
	  D2: year 2+    (dial_mort)  — waitlist removals and graft failures enter here,
	                                skipping D1 (matching esrd_time >= 1 in the sim).
	Runs until natural extinction (no hard age cap).
	"""
	lt = life_table if life_table is not None else LIFE_TABLE_QX

	wl_tx         = waitlist_annual_tx_prob(p, priority=donor)
	wl_mort_p     = waitlist_annual_mort(p)
	wl_remove     = float(p["wl_removal_rate_yr"])
	wl_listing    = float(p.get("wl_listing_prob", 1.0))
	dial_mort1    = float(p["dialysis_1yr_mort"])
	dial_mort     = float(p["dialysis_annual_mort"])
	bg_hr         = float(p.get("donor_mort_hr", 1.0)) if donor else 1.0
	graft_fail    = float(p.get("graft_annual_fail_postyear1", 0.025))
	preemptive_p  = float(p.get(
		"esrd_preemptive_prob_pld" if donor else "esrd_preemptive_prob_std", 0.0))

	cum_risk_15 = float(p["esrd_15yr_donor_overall"] if donor else p["esrd_15yr_nondonor"])
	wbl_k   = float(p["weibull_shape"])
	wbl_lam = weibull_scale_from_cumrisk_competing(
		cum_risk_15, wbl_k, lt, age_at_entry, bg_hr
	)

	H, D1, D2, WL, PT = 1.0, 0.0, 0.0, 0.0, 0.0
	total_ly = 0.0
	yr = 0

	while H + D1 + D2 + WL + PT > 1e-9:
		age    = age_at_entry + yr
		q_bg   = lt[min(age, MAX_AGE)] * bg_hr
		p_esrd = weibull_annual_prob(float(yr), wbl_lam, wbl_k)

		if age < 35:   ptx_mort = float(p.get("posttx_annual_mort_age1834", p["posttx_annual_mort"]))
		elif age < 50: ptx_mort = float(p.get("posttx_annual_mort_age3549", p["posttx_annual_mort"]))
		elif age < 65: ptx_mort = float(p.get("posttx_annual_mort_age5064", p["posttx_annual_mort"]))
		else:          ptx_mort = float(p.get("posttx_annual_mort_age65p",  p["posttx_annual_mort"]))

		H_die  = H * q_bg;     H_surv = H - H_die
		H_esrd = H_surv * p_esrd;      H_stay = H_surv - H_esrd
		# At ESRD onset, branch: preemptive → WL (skip dialysis), rest → D1
		H_preemptive = H_esrd * preemptive_p
		H_to_D1      = H_esrd - H_preemptive

		D1_die    = D1 * dial_mort1;   D1_surv   = D1 - D1_die
		D1_listed = D1_surv * wl_listing;         D1_to_D2  = D1_surv - D1_listed

		D2_die    = D2 * dial_mort;    D2_surv   = D2 - D2_die
		D2_listed = D2_surv * wl_listing;         D2_stay   = D2_surv - D2_listed

		WL_die    = WL * wl_mort_p;    WL_surv   = WL - WL_die
		WL_tx     = WL_surv * wl_tx;   WL_after  = WL_surv - WL_tx
		WL_remove = WL_after * wl_remove;         WL_stay   = WL_after - WL_remove

		PT_die  = PT * ptx_mort;       PT_surv = PT - PT_die
		PT_fail = PT_surv * graft_fail;            PT_stay = PT_surv - PT_fail

		H  = H_stay
		D1 = H_to_D1
		D2 = D1_to_D2 + D2_stay + WL_remove + PT_fail
		WL = D1_listed + D2_listed + WL_stay + H_preemptive
		PT = WL_tx + PT_stay

		total_ly += (H + D1 + D2 + WL + PT) * CYCLE_YRS
		yr += 1
		if yr > 300:
			break

	return total_ly


# ── RUN BASE CASE ─────────────────────────────────────────────────────────────
def run_base_case(age_at_donation=40, n=N_SIM):
	print(f"\nRunning base case — age at donation {age_at_donation}, n={n:,}")
	p = BASE.copy()

	ly_donor    = simulate_cohort(p, n, age_at_donation, donor=True,
								  rng=np.random.default_rng(99))
	ly_nondonor = simulate_cohort(p, n, age_at_donation, donor=False,
								  rng=np.random.default_rng(99_000))

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


# ── PROBABILISTIC SENSITIVITY ANALYSIS ───────────────────────────────────────
def run_psa(age_at_donation=40, n_draws=N_DRAWS):
	print(f"\nRunning PSA (analytic) — {n_draws} parameter draws...")
	rng = np.random.default_rng(77)
	diffs = []

	for i in range(n_draws):
		if (i+1) % 50 == 0:
			print(f"  PSA draw {i+1}/{n_draws}")
		p = sample_params(rng)
		diffs.append(run_arm_analytic(p, age_at_donation, donor=True)
					 - run_arm_analytic(p, age_at_donation, donor=False))

	diffs = np.array(diffs)
	print(f"  PSA ΔLE: mean={diffs.mean():+.3f}, "
		  f"95% CrI [{np.percentile(diffs,2.5):+.3f}, {np.percentile(diffs,97.5):+.3f}]")
	print(f"  P(donation beneficial): {(diffs > 0).mean():.1%}")
	return diffs


# ── ONE-WAY SENSITIVITY ANALYSIS ─────────────────────────────────────────────
def run_owsa(age_at_donation=40):
	print("\nRunning one-way sensitivity analysis (analytic)...")

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
		# Post-Tx survival: LDKT quality vs base-case DDKT (SRTR 2023 ADR Fig KI 76)
		"Post-Tx: LDKT quality":         {
			"posttx_annual_mort_age1834": BASE.get("posttx_ld_annual_mort_age1834", 0.0042),
			"posttx_annual_mort_age3549": BASE.get("posttx_ld_annual_mort_age3549", 0.0079),
			"posttx_annual_mort_age5064": BASE.get("posttx_ld_annual_mort_age5064", 0.0172),
			"posttx_annual_mort_age65p":  BASE.get("posttx_ld_annual_mort_age65p",  0.0392),
			"posttx_annual_mort":         BASE.get("posttx_ld_annual_mort",          0.0175),
		},
	}

	results = {}
	for label, overrides in scenarios.items():
		p = BASE.copy()
		p.update(overrides)
		diff = run_arm_analytic(p, age_at_donation, donor=True) \
			 - run_arm_analytic(p, age_at_donation, donor=False)
		results[label] = diff
		print(f"  {label:<44} ΔLE = {diff:+.3f} yr  ({diff*365.25:+.1f} days)")

	return results


# ── AGE SUBGROUP ANALYSIS ─────────────────────────────────────────────────────
def run_age_subgroups():
	print("\nRunning age-at-donation subgroup analysis (analytic)...")
	ages = [25, 35, 40, 45, 55]
	results = {}
	for age in ages:
		p    = BASE.copy()
		p["esrd_15yr_donor_overall"] = _age_adjust_esrd_nonblack(
			BASE["esrd_15yr_donor_overall"], age)
		le_d  = run_arm_analytic(p, age, donor=True)
		le_nd = run_arm_analytic(p, age, donor=False)
		diff  = le_d - le_nd
		results[age] = {"le_donor": le_d, "le_nondonor": le_nd, "diff": diff}
		print(f"  Age {age}: ΔLE = {diff:+.3f} yr  "
			  f"(donor {le_d:.1f} yr, non-donor {le_nd:.1f} yr)")
	return results


# ── SHARED SUBGROUP HELPERS ───────────────────────────────────────────────────
def _scale_posttx_mort(base_params: dict, target_overall: float) -> dict:
	"""Scale age-specific post-tx mortality rates proportionally to target_overall."""
	base_overall = base_params.get("posttx_annual_mort", 0.036)
	if base_overall <= 0:
		return {"posttx_annual_mort": target_overall}
	scale = target_overall / base_overall
	return {
		"posttx_annual_mort":         target_overall,
		"posttx_annual_mort_age1834": base_params.get("posttx_annual_mort_age1834", 0.009) * scale,
		"posttx_annual_mort_age3549": base_params.get("posttx_annual_mort_age3549", 0.018) * scale,
		"posttx_annual_mort_age5064": base_params.get("posttx_annual_mort_age5064", 0.039) * scale,
		"posttx_annual_mort_age65p":  base_params.get("posttx_annual_mort_age65p",  0.068) * scale,
	}


def _age_adjust_esrd_nonblack(base_rate: float, age: int, reference_age: int = 40) -> float:
	"""Scale non-Black donor ESRD rate by Massie 2017 HR per decade (1.40/decade)."""
	hr = float(BASE.get("hr_age_per_decade_nonblack", 1.40))
	return base_rate * hr ** ((age - reference_age) / 10)


# Grams 2016 direct sex×race non-donor baselines (Table 2, NEJM 374:411).
# Option C sourcing rule: use this table wherever sex or race is stratified.
# For race-unspecified ("Overall") sex analyses we use the White sex-specific
# rates as a conservative approximation — White donors dominate the SRTR pool
# (~70%) and Grams White female (0.04%) ≈ Muzaale matched-control overall
# (0.039%), so the overall and sex-only panels remain comparable.
_GRAMS_NONDONOR = {
	("Black",   "Female"): 0.0015,   # 0.15%
	("Black",   "Male"):   0.0024,   # 0.24%
	("White",   "Female"): 0.0004,   # 0.04%
	("White",   "Male"):   0.0006,   # 0.06%
	# Sex-stratified, race-unspecified: White rates as donor-pool approximation
	("Overall", "Female"): 0.0004,
	("Overall", "Male"):   0.0006,
}


# ── RACE SUBGROUP ANALYSIS ────────────────────────────────────────────────────
def run_race_subgroups(age_at_donation=40, n=500_000):
	print("\nRunning race subgroup analysis...")
	scenarios = {
		"White donor": {
			"esrd_15yr_donor_overall": 0.00227,
			"esrd_15yr_nondonor":      0.00050,
			"wl_mort_per_100py":       BASE.get("wl_mort_white_per_100py", 5.71),
			**_scale_posttx_mort(BASE, BASE.get("posttx_annual_mort_white", 0.038)),
		},
		"Overall (base)": {},
		"Black donor": {
			"esrd_15yr_donor_overall": 0.00747,
			"esrd_15yr_nondonor":      0.00195,
			"wl_mort_per_100py":       BASE.get("wl_mort_black_per_100py", 4.62),
			**_scale_posttx_mort(BASE, BASE.get("posttx_annual_mort_black", 0.035)),
		},
	}
	results = {}
	for i, (label, overrides) in enumerate(scenarios.items()):
		p = BASE.copy()
		p.update(overrides)
		# Analytic ΔLE for panel E
		diff = run_arm_analytic(p, age_at_donation, donor=True) \
			 - run_arm_analytic(p, age_at_donation, donor=False)
		# Simulation draws for the separate race histogram figure
		ld  = simulate_cohort(p, n, age_at_donation, donor=True,
							  rng=np.random.default_rng(42  + i * 100))
		lnd = simulate_cohort(p, n, age_at_donation, donor=False,
							  rng=np.random.default_rng(4200 + i * 100))
		results[label] = {
			"le_donor":    ld.mean(),
			"le_nondonor": lnd.mean(),
			"diff":        diff,
			"ly_donor":    ld,
			"ly_nondonor": lnd,
		}
		print(f"  {label}: ΔLE = {diff:+.3f} yr  ({diff*365.25:+.1f} days)")
	return results


# ── SEX SUBGROUP ANALYSIS ─────────────────────────────────────────────────────
def run_sex_subgroups(age_at_donation=40, n=500_000):
	print("\nRunning sex subgroup analysis...")

	lt_m_path = DATA_PROC / "lifetable_male_2021.csv"
	lt_f_path = DATA_PROC / "lifetable_female_2021.csv"
	lt_male   = load_life_table(lt_m_path if lt_m_path.exists() else None)
	lt_female = load_life_table(lt_f_path if lt_f_path.exists() else None)

	# Donor ESRD rates: scaled via Massie 2017 within-donor sex HR (1.88).
	# Sex mix: 60% female — SRTR 2023 ADR Figure KI 7 (59–61%, 2019–2023).
	# Non-donor rates: Grams 2016 sex-specific direct values (Option C).
	# Using Massie HR to scale the non-donor arm is inappropriate because that
	# HR was estimated within the donor cohort, not the general population.
	f_female = 0.60
	f_male   = 1.0 - f_female
	hr_male  = float(BASE.get("hr_male_sex", 1.88))
	denom    = f_female + f_male * hr_male

	donor_f = BASE["esrd_15yr_donor_overall"] / denom
	donor_m = donor_f * hr_male

	scenarios = {
		"Female donor": {
			"esrd_15yr_donor_overall": donor_f,
			"esrd_15yr_nondonor":      _GRAMS_NONDONOR[("Overall", "Female")],
			"_lt": lt_female,
		},
		"Overall (base)": {"_lt": LIFE_TABLE_QX},
		"Male donor": {
			"esrd_15yr_donor_overall": donor_m,
			"esrd_15yr_nondonor":      _GRAMS_NONDONOR[("Overall", "Male")],
			"_lt": lt_male,
		},
	}

	results = {}
	for i, (label, spec) in enumerate(scenarios.items()):
		lt       = spec.get("_lt")
		overrides = {k: v for k, v in spec.items() if k != "_lt"}
		p = BASE.copy()
		p.update(overrides)
		diff = run_arm_analytic(p, age_at_donation, donor=True,  life_table=lt) \
			 - run_arm_analytic(p, age_at_donation, donor=False, life_table=lt)
		ld  = simulate_cohort(p, n, age_at_donation, donor=True,
							  rng=np.random.default_rng(52  + i * 100),
							  life_table=lt)
		lnd = simulate_cohort(p, n, age_at_donation, donor=False,
							  rng=np.random.default_rng(5200 + i * 100),
							  life_table=lt)
		results[label] = {
			"le_donor":    ld.mean(),
			"le_nondonor": lnd.mean(),
			"diff":        diff,
			"ly_donor":    ld,
			"ly_nondonor": lnd,
		}
		print(f"  {label}: ΔLE = {diff:+.3f} yr  ({diff*365.25:+.1f} days)")

	return results


# ── PLOTTING HELPERS ──────────────────────────────────────────────────────────
_TEAL   = "#1D9E75"
_CORAL  = "#D85A30"
_PURPLE = "#534AB7"
_AMBER  = "#BA7517"
_GRAY   = "#888780"
_LIGHT  = "#F1EFE8"


def _style_ax(ax, title, fontsize=10):
	ax.set_facecolor(_LIGHT)
	ax.spines[["top", "right"]].set_visible(False)
	ax.spines[["left", "bottom"]].set_color("#B4B2A9")
	ax.tick_params(colors="#5F5E5A", labelsize=9)
	ax.set_title(title, fontsize=fontsize, fontweight="bold", color="#2C2C2A", pad=7)


# ── (A) LE DISTRIBUTIONS ──────────────────────────────────────────────────────
def make_fig_distributions(base_res):
	age = base_res["age"]
	fig, ax = plt.subplots(figsize=(10, 4))
	fig.patch.set_facecolor("#FAFAF8")
	bins = np.arange(0, MAX_AGE - age + 2, 2)
	ax.hist(base_res["ly_nondonor"], bins=bins, density=True,
			color=_TEAL,  alpha=0.55, label=f"Non-donor  μ={base_res['le_nondonor']:.1f} yr")
	ax.hist(base_res["ly_donor"],    bins=bins, density=True,
			color=_CORAL, alpha=0.55, label=f"Donor      μ={base_res['le_donor']:.1f} yr")
	ax.axvline(base_res["le_nondonor"], color=_TEAL,  lw=1.8, ls="--")
	ax.axvline(base_res["le_donor"],    color=_CORAL, lw=1.8, ls="--")
	ax.legend(fontsize=9, frameon=False)
	ax.set_xlabel("Remaining life-years from donation age", fontsize=9)
	ax.set_ylabel("Density", fontsize=9)
	_style_ax(ax, f"Lifetime distributions — base case (age {age}, n={N_SIM//1_000_000}M/arm)",
			  fontsize=11)
	fig.tight_layout()
	return fig


# ── (B) PSA ───────────────────────────────────────────────────────────────────
def make_fig_psa(psa_diffs):
	fig, ax = plt.subplots(figsize=(6, 4))
	fig.patch.set_facecolor("#FAFAF8")
	ax.hist(psa_diffs * 365.25, bins=40, color=_PURPLE, alpha=0.75,
			edgecolor="white", linewidth=0.3)
	ax.axvline(0, color=_GRAY, lw=1.5, ls="--")
	lo, hi = np.percentile(psa_diffs * 365.25, [2.5, 97.5])
	ax.axvline(lo, color=_AMBER, lw=1, ls=":")
	ax.axvline(hi, color=_AMBER, lw=1, ls=":")
	p_ben = (psa_diffs > 0).mean()
	ax.set_xlabel("ΔLE (days): donor − non-donor", fontsize=9)
	ax.set_ylabel("Count", fontsize=9)
	_style_ax(ax, f"PSA ({N_DRAWS} draws)\nP(beneficial)={p_ben:.0%}, 95% CrI [{lo:.0f}, {hi:.0f}] d",
			  fontsize=10)
	fig.tight_layout()
	return fig


# ── (C) TORNADO / OWSA ────────────────────────────────────────────────────────
def make_fig_tornado(owsa_res):
	labels   = list(owsa_res.keys())
	vals     = [owsa_res[lbl] * 365.25 for lbl in labels]
	colors_c = [_TEAL if v >= 0 else _CORAL for v in vals]
	sorted_idx = np.argsort(np.abs(vals))
	labels_s   = [labels[i]   for i in sorted_idx]
	vals_s     = [vals[i]     for i in sorted_idx]
	colors_s   = [colors_c[i] for i in sorted_idx]
	y_pos      = np.arange(len(labels_s))

	abs_s         = np.abs(vals_s)
	dominant_val  = vals_s[-1]
	dominant_sign = int(np.sign(dominant_val))
	break_abs     = abs_s[-2] * 1.35
	break_val     = dominant_sign * break_abs

	fig = plt.figure(figsize=(14, 5))
	fig.patch.set_facecolor("#FAFAF8")
	gs = gridspec.GridSpec(1, 2, figure=fig,
						   width_ratios=[1, 3] if dominant_sign < 0 else [3, 1],
						   wspace=0.03)

	if dominant_sign < 0:
		ax_ext  = fig.add_subplot(gs[0])
		ax_main = fig.add_subplot(gs[1])
		ext_xlim  = (dominant_val * 1.08, break_val)
		main_xlim = (break_val, 0)
	else:
		ax_main = fig.add_subplot(gs[0])
		ax_ext  = fig.add_subplot(gs[1])
		neg_vals = [v for v in vals_s if v < 0]
		ext_xlim  = (break_val, dominant_val * 1.08)
		main_xlim = (min(neg_vals) * 1.25 if neg_vals else -break_abs * 0.3, break_val)

	for ax, xlim in ((ax_ext, ext_xlim), (ax_main, main_xlim)):
		ax.barh(y_pos, vals_s, color=colors_s, alpha=0.8, height=0.65)
		ax.set_xlim(*xlim)
		ax.set_facecolor(_LIGHT)
		ax.spines["top"].set_visible(False)
		ax.spines["bottom"].set_color("#B4B2A9")
		ax.tick_params(colors="#5F5E5A", labelsize=9)

	ax_main.set_yticks(y_pos)
	ax_main.set_yticklabels(labels_s, fontsize=9)
	ax_main.yaxis.tick_right()
	ax_ext.set_yticks([])

	if dominant_sign < 0:
		ax_ext.spines["right"].set_visible(False)
		ax_ext.spines["left"].set_color("#B4B2A9")
		ax_main.spines["left"].set_visible(False)
		ax_main.spines["right"].set_color(_GRAY)
	else:
		ax_main.spines["right"].set_visible(False)
		ax_main.spines["left"].set_color("#B4B2A9")
		ax_ext.spines["left"].set_visible(False)
		ax_ext.spines["right"].set_color("#B4B2A9")

	d = 0.03
	bkw = dict(color="#5F5E5A", clip_on=False, lw=1.5, zorder=5)
	if dominant_sign < 0:
		for y0 in (0, 1):
			ax_ext.plot((1 - d, 1 + d), (y0 - d, y0 + d),
						transform=ax_ext.transAxes, **bkw)
			ax_main.plot((-d, +d), (y0 - d, y0 + d),
						 transform=ax_main.transAxes, **bkw)
	else:
		for y0 in (0, 1):
			ax_main.plot((1 - d, 1 + d), (y0 - d, y0 + d),
						 transform=ax_main.transAxes, **bkw)
			ax_ext.plot((-d, +d), (y0 - d, y0 + d),
						transform=ax_ext.transAxes, **bkw)

	ax_ext.text(
		(break_val + dominant_val) / 2, len(labels_s) - 1,
		f"{dominant_val:+.0f} d",
		va="center", ha="center", fontsize=8.5, color="white", fontweight="bold"
	)

	ax_ext.set_xlabel("ΔLE (days)", fontsize=9)
	ax_ext.set_facecolor("#EDE9E0")
	ax_main.set_xlabel("ΔLE (days): donor − non-donor", fontsize=9)
	ax_main.set_title("One-way sensitivity analysis — tornado chart",
					  fontsize=11, fontweight="bold", color="#2C2C2A", pad=7)
	return fig


# ── (D) AGE SUBGROUP ─────────────────────────────────────────────────────────
def make_fig_age_subgroup(age_res):
	fig, ax = plt.subplots(figsize=(8, 4))
	fig.patch.set_facecolor("#FAFAF8")
	ages_s  = sorted(age_res.keys())
	diffs_d = [age_res[a]["diff"] * 365.25 for a in ages_s]
	bar_cols = [_TEAL if d >= 0 else _CORAL for d in diffs_d]
	ax.bar(range(len(ages_s)), diffs_d, color=bar_cols, alpha=0.8, width=0.6)
	ax.axhline(0, color=_GRAY, lw=1.2)
	ax.set_xticks(range(len(ages_s)))
	ax.set_xticklabels([f"Age {a}" for a in ages_s], fontsize=9)
	ax.set_ylabel("ΔLE (days)", fontsize=9)
	_style_ax(ax, "ΔLE by age at donation", fontsize=11)
	fig.tight_layout()
	return fig


# ── (E) RACE SUBGROUP ────────────────────────────────────────────────────────
def make_fig_race_subgroup(race_res):
	fig, ax = plt.subplots(figsize=(6, 4))
	fig.patch.set_facecolor("#FAFAF8")
	race_labels = list(race_res.keys())
	race_diffs  = [race_res[r]["diff"] * 365.25 for r in race_labels]
	rc = [_TEAL if d >= 0 else _CORAL for d in race_diffs]
	ax.bar(range(len(race_labels)), race_diffs, color=rc, alpha=0.8, width=0.5)
	ax.axhline(0, color=_GRAY, lw=1.2)
	ax.set_xticks(range(len(race_labels)))
	ax.set_xticklabels(race_labels, fontsize=8.5)
	ax.set_ylabel("ΔLE (days)", fontsize=9)
	_style_ax(ax, "ΔLE by race", fontsize=11)
	fig.tight_layout()
	return fig


def make_race_figure(race_res):
	"""LE distributions by race group — one row per race, donor vs non-donor."""
	TEAL  = "#1D9E75"
	CORAL = "#D85A30"
	LIGHT = "#F1EFE8"
	GRAY  = "#888780"

	labels = list(race_res.keys())
	fig, axes = plt.subplots(1, len(labels), figsize=(5 * len(labels), 5),
							 sharey=False)
	fig.patch.set_facecolor("#FAFAF8")

	for ax, label in zip(axes, labels):
		res  = race_res[label]
		bins = np.arange(0, MAX_AGE - 40 + 2, 2)
		ax.hist(res["ly_nondonor"], bins=bins, density=True,
				color=TEAL,  alpha=0.55,
				label=f"Non-donor  μ={res['le_nondonor']:.1f} yr")
		ax.hist(res["ly_donor"], bins=bins, density=True,
				color=CORAL, alpha=0.55,
				label=f"Donor      μ={res['le_donor']:.1f} yr")
		ax.axvline(res["le_nondonor"], color=TEAL,  lw=1.8, ls="--")
		ax.axvline(res["le_donor"],    color=CORAL, lw=1.8, ls="--")
		diff_d = res["diff"] * 365.25
		ax.set_title(f"{label}\nΔLE = {diff_d:+.1f} days",
					 fontsize=10, fontweight="bold", color="#2C2C2A")
		ax.set_xlabel("Remaining life-years from age 40", fontsize=9)
		ax.set_ylabel("Density", fontsize=9)
		ax.legend(fontsize=8, frameon=False)
		ax.set_facecolor(LIGHT)
		ax.spines[["top","right"]].set_visible(False)
		ax.spines[["left","bottom"]].set_color("#B4B2A9")
		ax.tick_params(colors="#5F5E5A", labelsize=9)

	plt.suptitle(
		"Living Kidney Donation — ΔLE by Race (age 40)\n"
		"Donor vs Matched Non-Donor Life Expectancy",
		fontsize=12, fontweight="bold", color="#2C2C2A", y=1.02)
	plt.tight_layout()
	return fig


def make_sex_figure(sex_res):
	"""Bar chart of ΔLE by sex — cohort simulation means, styled like panel D."""
	TEAL  = "#1D9E75"
	CORAL = "#D85A30"
	LIGHT = "#F1EFE8"
	GRAY  = "#888780"

	labels = list(sex_res.keys())
	# Use cohort simulation means, not the analytic diff
	diffs  = [(sex_res[s]["le_donor"] - sex_res[s]["le_nondonor"]) * 365.25
			  for s in labels]
	colors = [TEAL if d >= 0 else CORAL for d in diffs]

	fig, ax = plt.subplots(figsize=(6, 4))
	fig.patch.set_facecolor("#FAFAF8")
	ax.set_facecolor(LIGHT)

	bars = ax.bar(range(len(labels)), diffs, color=colors, alpha=0.8, width=0.5)
	ax.axhline(0, color=GRAY, lw=1.2)

	ax.set_xticks(range(len(labels)))
	ax.set_xticklabels(labels, fontsize=10)
	ax.set_ylabel("ΔLE (days): donor − non-donor", fontsize=9)
	ax.spines[["top", "right"]].set_visible(False)
	ax.spines[["left", "bottom"]].set_color("#B4B2A9")
	ax.tick_params(colors="#5F5E5A", labelsize=9)
	ax.set_title(
		"ΔLE by sex at donation (age 40)\n"
		"Cohort Markov simulation  |  donor vs matched non-donor",
		fontsize=10, fontweight="bold", color="#2C2C2A", pad=8)

	fig.tight_layout()
	return fig


# ── AGE × RACE MATRIX FIGURE ─────────────────────────────────────────────────
def make_age_race_matrix():
	"""
	Heatmap of ΔLE (days) across a grid of donation ages × race groups,
	computed analytically. Each cell is exact — no Monte Carlo noise.
	"""
	ages = [25, 30, 35, 40, 45, 50, 55]

	race_params = {
		"White": {
			"esrd_15yr_donor_overall": 0.00227,
			"esrd_15yr_nondonor":      0.00050,
			"wl_mort_per_100py":       BASE.get("wl_mort_white_per_100py", 5.71),
			**_scale_posttx_mort(BASE, BASE.get("posttx_annual_mort_white", 0.038)),
		},
		"Overall": {},
		"Black": {
			"esrd_15yr_donor_overall": 0.00747,
			"esrd_15yr_nondonor":      0.00195,
			"wl_mort_per_100py":       BASE.get("wl_mort_black_per_100py", 4.62),
			**_scale_posttx_mort(BASE, BASE.get("posttx_annual_mort_black", 0.035)),
		},
	}

	races = list(race_params.keys())
	matrix = np.zeros((len(races), len(ages)))   # rows = race, cols = age

	for r, race in enumerate(races):
		p = BASE.copy()
		p.update(race_params[race])
		base_donor = p["esrd_15yr_donor_overall"]
		for a, age in enumerate(ages):
			if race != "Black":
				p["esrd_15yr_donor_overall"] = _age_adjust_esrd_nonblack(base_donor, age)
			else:
				p["esrd_15yr_donor_overall"] = base_donor
			le_d  = run_arm_analytic(p, age, donor=True)
			le_nd = run_arm_analytic(p, age, donor=False)
			matrix[r, a] = (le_d - le_nd) * 365.25

	# ── Plot ──────────────────────────────────────────────────────────────────
	fig, ax = plt.subplots(figsize=(10, 4))
	fig.patch.set_facecolor("#FAFAF8")
	ax.set_facecolor("#FAFAF8")

	abs_max = np.abs(matrix).max()
	im = ax.imshow(matrix, cmap="RdYlGn", vmin=-abs_max, vmax=abs_max,
				   aspect="auto")

	# Cell labels
	for r in range(len(races)):
		for a in range(len(ages)):
			val = matrix[r, a]
			color = "black" if abs(val) < 0.6 * abs_max else "white"
			ax.text(a, r, f"{val:+.0f} d", ha="center", va="center",
					fontsize=11, fontweight="bold", color=color)

	ax.set_xticks(range(len(ages)))
	ax.set_xticklabels([f"Age {a}" for a in ages], fontsize=10)
	ax.set_yticks(range(len(races)))
	ax.set_yticklabels(races, fontsize=11, fontweight="bold")
	ax.set_xlabel("Age at donation", fontsize=11)

	cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
	cbar.set_label("ΔLE (days): donor − non-donor", fontsize=9)

	ax.set_title(
		"Life-expectancy impact of kidney donation by age and race\n"
		"Donor vs matched non-donor  |  analytic cohort Markov  |  green = donor benefit, red = donor harm",
		fontsize=11, fontweight="bold", color="#2C2C2A", pad=10)

	ax.spines[:].set_visible(False)
	ax.tick_params(length=0)
	plt.tight_layout()
	return fig


# ── AGE × SEX MATRIX FIGURE ──────────────────────────────────────────────────
def make_age_sex_matrix():
	"""Heatmap of ΔLE (days) across donation ages × sex groups."""
	ages = [25, 30, 35, 40, 45, 50, 55]

	lt_m_path = DATA_PROC / "lifetable_male_2021.csv"
	lt_f_path = DATA_PROC / "lifetable_female_2021.csv"
	lt_male   = load_life_table(lt_m_path if lt_m_path.exists() else None)
	lt_female = load_life_table(lt_f_path if lt_f_path.exists() else None)

	# Donor rates: Massie 2017 within-donor sex HR. Non-donor rates: Grams 2016
	# sex-specific direct values (Option C — same sourcing as sex×race matrix).
	f_female = 0.60  # SRTR 2023 ADR Figure KI 7: 59-61% of LD donors female
	hr_male  = float(BASE.get("hr_male_sex", 1.88))
	denom    = f_female + (1.0 - f_female) * hr_male

	sex_params = {
		"Female": {
			"esrd_15yr_donor_overall": BASE["esrd_15yr_donor_overall"] / denom,
			"esrd_15yr_nondonor":      _GRAMS_NONDONOR[("Overall", "Female")],
			"_lt": lt_female,
		},
		"Overall": {"_lt": LIFE_TABLE_QX},
		"Male": {
			"esrd_15yr_donor_overall": BASE["esrd_15yr_donor_overall"] * hr_male / denom,
			"esrd_15yr_nondonor":      _GRAMS_NONDONOR[("Overall", "Male")],
			"_lt": lt_male,
		},
	}

	sexes  = list(sex_params.keys())
	matrix = np.zeros((len(sexes), len(ages)))

	for s, sex in enumerate(sexes):
		spec = sex_params[sex]
		lt   = spec.get("_lt")
		overrides = {k: v for k, v in spec.items() if k != "_lt"}
		p = BASE.copy()
		p.update(overrides)
		base_donor = p["esrd_15yr_donor_overall"]
		for a, age in enumerate(ages):
			p["esrd_15yr_donor_overall"] = _age_adjust_esrd_nonblack(base_donor, age)
			le_d  = run_arm_analytic(p, age, donor=True,  life_table=lt)
			le_nd = run_arm_analytic(p, age, donor=False, life_table=lt)
			matrix[s, a] = (le_d - le_nd) * 365.25

	fig, ax = plt.subplots(figsize=(10, 3))
	fig.patch.set_facecolor("#FAFAF8")
	ax.set_facecolor("#FAFAF8")

	abs_max = np.abs(matrix).max()
	im = ax.imshow(matrix, cmap="RdYlGn", vmin=-abs_max, vmax=abs_max, aspect="auto")

	for s in range(len(sexes)):
		for a in range(len(ages)):
			val = matrix[s, a]
			color = "black" if abs(val) < 0.6 * abs_max else "white"
			ax.text(a, s, f"{val:+.0f} d", ha="center", va="center",
					fontsize=11, fontweight="bold", color=color)

	ax.set_xticks(range(len(ages)))
	ax.set_xticklabels([f"Age {a}" for a in ages], fontsize=10)
	ax.set_yticks(range(len(sexes)))
	ax.set_yticklabels(sexes, fontsize=11, fontweight="bold")
	ax.set_xlabel("Age at donation", fontsize=11)

	cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
	cbar.set_label("ΔLE (days): donor − non-donor", fontsize=9)

	ax.set_title(
		"Life-expectancy impact of kidney donation by age and sex\n"
		"Donor vs matched non-donor  |  analytic cohort Markov  |  green = donor benefit, red = donor harm",
		fontsize=11, fontweight="bold", color="#2C2C2A", pad=10)

	ax.spines[:].set_visible(False)
	ax.tick_params(length=0)
	plt.tight_layout()
	return fig


# ── SEX × RACE MATRIX FIGURE ─────────────────────────────────────────────────
def make_sex_race_matrix():
	"""Heatmap of ΔLE (days) across sex groups × race groups at age 40."""
	lt_m_path = DATA_PROC / "lifetable_male_2021.csv"
	lt_f_path = DATA_PROC / "lifetable_female_2021.csv"
	lt_male   = load_life_table(lt_m_path if lt_m_path.exists() else None)
	lt_female = load_life_table(lt_f_path if lt_f_path.exists() else None)

	f_female = 0.60  # SRTR 2023 ADR Figure KI 7: 59-61% of LD donors female
	hr_male  = float(BASE.get("hr_male_sex", 1.88))
	denom    = f_female + (1.0 - f_female) * hr_male

	# Race-specific sex-averaged ESRD rates + other race overrides
	race_base = {
		"White": {
			"esrd_donor":    0.00227,
			"esrd_nondonor": 0.00050,
			"overrides": {
				"wl_mort_per_100py": BASE.get("wl_mort_white_per_100py", 5.71),
				**_scale_posttx_mort(BASE, BASE.get("posttx_annual_mort_white", 0.038)),
			},
		},
		"Overall": {
			"esrd_donor":    BASE["esrd_15yr_donor_overall"],
			"esrd_nondonor": BASE["esrd_15yr_nondonor"],
			"overrides": {},
		},
		"Black": {
			"esrd_donor":    0.00747,
			"esrd_nondonor": 0.00195,
			"overrides": {
				"wl_mort_per_100py": BASE.get("wl_mort_black_per_100py", 4.62),
				**_scale_posttx_mort(BASE, BASE.get("posttx_annual_mort_black", 0.035)),
			},
		},
	}

	# Sex scaling factors applied to each race's ESRD rates
	sex_specs = {
		"Female":  {"scale": 1.0 / denom,          "_lt": lt_female},
		"Overall": {"scale": 1.0,                   "_lt": LIFE_TABLE_QX},
		"Male":    {"scale": hr_male / denom,        "_lt": lt_male},
	}

	races  = list(race_base.keys())
	sexes  = list(sex_specs.keys())
	age    = 40
	matrix = np.zeros((len(races), len(sexes)))  # rows = race, cols = sex

	for r, race in enumerate(races):
		rb = race_base[race]
		for s, sex in enumerate(sexes):
			ss = sex_specs[sex]
			p  = BASE.copy()
			p.update(rb["overrides"])
			p["esrd_15yr_donor_overall"] = rb["esrd_donor"] * ss["scale"]
			# Use Grams 2016 direct value where published; otherwise derive from sex HR
			p["esrd_15yr_nondonor"] = _GRAMS_NONDONOR.get(
				(race, sex), rb["esrd_nondonor"] * ss["scale"])
			lt = ss["_lt"]
			le_d  = run_arm_analytic(p, age, donor=True,  life_table=lt)
			le_nd = run_arm_analytic(p, age, donor=False, life_table=lt)
			matrix[r, s] = (le_d - le_nd) * 365.25

	fig, ax = plt.subplots(figsize=(6, 4))
	fig.patch.set_facecolor("#FAFAF8")
	ax.set_facecolor("#FAFAF8")

	abs_max = np.abs(matrix).max()
	im = ax.imshow(matrix, cmap="RdYlGn", vmin=-abs_max, vmax=abs_max, aspect="auto")

	for r in range(len(races)):
		for s in range(len(sexes)):
			val = matrix[r, s]
			color = "black" if abs(val) < 0.6 * abs_max else "white"
			ax.text(s, r, f"{val:+.0f} d", ha="center", va="center",
					fontsize=11, fontweight="bold", color=color)

	ax.set_xticks(range(len(sexes)))
	ax.set_xticklabels(sexes, fontsize=11, fontweight="bold")
	ax.set_yticks(range(len(races)))
	ax.set_yticklabels(races, fontsize=11, fontweight="bold")
	ax.set_xlabel("Sex", fontsize=11)

	cbar = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.02)
	cbar.set_label("ΔLE (days): donor − non-donor", fontsize=9)

	ax.set_title(
		"Life-expectancy impact of kidney donation by race and sex (age 40)\n"
		"Donor vs matched non-donor  |  analytic cohort Markov  |  green = donor benefit, red = donor harm",
		fontsize=11, fontweight="bold", color="#2C2C2A", pad=10)

	ax.spines[:].set_visible(False)
	ax.tick_params(length=0)
	plt.tight_layout()
	return fig


# ── AGE × RACE × SEX FACETED MATRIX ──────────────────────────────────────────
def make_age_race_sex_matrix():
	"""
	Three age×race heatmaps side by side — one panel per sex group.
	All panels share the same colorscale for direct comparison.
	"""
	ages = [25, 30, 35, 40, 45, 50, 55]

	lt_m_path = DATA_PROC / "lifetable_male_2021.csv"
	lt_f_path = DATA_PROC / "lifetable_female_2021.csv"
	lt_male   = load_life_table(lt_m_path if lt_m_path.exists() else None)
	lt_female = load_life_table(lt_f_path if lt_f_path.exists() else None)

	f_female = 0.60  # SRTR 2023 ADR Figure KI 7: 59-61% of LD donors female
	hr_male  = float(BASE.get("hr_male_sex", 1.88))
	denom    = f_female + (1.0 - f_female) * hr_male

	race_params = {
		"White": {
			"esrd_15yr_donor_overall": 0.00227,
			"esrd_15yr_nondonor":      0.00050,
			"wl_mort_per_100py":       BASE.get("wl_mort_white_per_100py", 5.71),
			**_scale_posttx_mort(BASE, BASE.get("posttx_annual_mort_white", 0.038)),
		},
		"Overall": {},
		"Black": {
			"esrd_15yr_donor_overall": 0.00747,
			"esrd_15yr_nondonor":      0.00195,
			"wl_mort_per_100py":       BASE.get("wl_mort_black_per_100py", 4.62),
			**_scale_posttx_mort(BASE, BASE.get("posttx_annual_mort_black", 0.035)),
		},
	}

	sex_specs = {
		"Female":  {"scale": 1.0 / denom,        "_lt": lt_female},
		"Overall": {"scale": 1.0,                 "_lt": LIFE_TABLE_QX},
		"Male":    {"scale": hr_male / denom,      "_lt": lt_male},
	}

	races = list(race_params.keys())
	sexes = list(sex_specs.keys())

	# Compute all matrices first so colorscale spans all three panels
	matrices = {}
	for sex in sexes:
		ss  = sex_specs[sex]
		lt  = ss["_lt"]
		mat = np.zeros((len(races), len(ages)))
		for r, race in enumerate(races):
			p = BASE.copy()
			p.update(race_params[race])
			base_donor    = p["esrd_15yr_donor_overall"] * ss["scale"]
			base_nondonor = p["esrd_15yr_nondonor"]      * ss["scale"]
			for a, age in enumerate(ages):
				# Age HR for non-Black donors (Massie 2017)
				if race != "Black":
					p["esrd_15yr_donor_overall"] = _age_adjust_esrd_nonblack(base_donor, age)
				else:
					p["esrd_15yr_donor_overall"] = base_donor
				# Grams 2016 direct non-donor baseline where published
				p["esrd_15yr_nondonor"] = _GRAMS_NONDONOR.get(
					(race, sex), base_nondonor)
				le_d  = run_arm_analytic(p, age, donor=True,  life_table=lt)
				le_nd = run_arm_analytic(p, age, donor=False, life_table=lt)
				mat[r, a] = (le_d - le_nd) * 365.25
		matrices[sex] = mat

	abs_max = max(np.abs(m).max() for m in matrices.values())

	fig, axes = plt.subplots(1, 3, figsize=(20, 4))
	fig.patch.set_facecolor("#FAFAF8")

	im = None
	for i, (ax, sex) in enumerate(zip(axes, sexes)):
		mat = matrices[sex]
		ax.set_facecolor("#FAFAF8")
		im = ax.imshow(mat, cmap="RdYlGn", vmin=-abs_max, vmax=abs_max, aspect="auto")

		for r in range(len(races)):
			for a in range(len(ages)):
				val = mat[r, a]
				color = "black" if abs(val) < 0.6 * abs_max else "white"
				ax.text(a, r, f"{val:+.0f} d", ha="center", va="center",
						fontsize=9, fontweight="bold", color=color)

		ax.set_xticks(range(len(ages)))
		ax.set_xticklabels([f"Age {a}" for a in ages], fontsize=9)
		ax.set_yticks(range(len(races)))
		ax.set_yticklabels(races if i == 0 else [], fontsize=10, fontweight="bold")
		ax.set_xlabel("Age at donation", fontsize=10)
		ax.set_title(sex, fontsize=13, fontweight="bold", color="#2C2C2A", pad=6)
		ax.spines[:].set_visible(False)
		ax.tick_params(length=0)

	cbar = fig.colorbar(im, ax=axes[2], fraction=0.05, pad=0.02)
	cbar.set_label("ΔLE (days): donor − non-donor", fontsize=9)

	fig.suptitle(
		"Life-expectancy impact of kidney donation by age, race, and sex\n"
		"Donor vs matched non-donor  |  analytic cohort Markov  |  green = donor benefit, red = donor harm",
		fontsize=11, fontweight="bold", color="#2C2C2A", y=1.05)

	plt.tight_layout()
	return fig


# ── RESULTS TABLE ─────────────────────────────────────────────────────────────
def make_results_table(base_res, psa_diffs, owsa_res, age_res, race_res, sex_res=None):
	rows = []

	rows.append({
		"Analysis": "Base case (age 40, overall)",
		"LE Donor (yr)": f"{base_res['le_donor']:.2f}",
		"LE Non-donor (yr)": f"{base_res['le_nondonor']:.2f}",
		"ΔLE (days)": f"{base_res['diff']*365.25:+.1f}",
		"Note": f"n={N_SIM//1_000_000}M/arm"
	})

	rows.append({
		"Analysis": "PSA (probabilistic)",
		"LE Donor (yr)": "—",
		"LE Non-donor (yr)": "—",
		"ΔLE (days)": (f"{psa_diffs.mean()*365.25:+.1f} "
					   f"[{np.percentile(psa_diffs*365.25,2.5):+.1f}, "
					   f"{np.percentile(psa_diffs*365.25,97.5):+.1f}]"),
		"Note": f"P(beneficial)={(psa_diffs>0).mean():.0%}"
	})

	for a, v in sorted(age_res.items()):
		rows.append({
			"Analysis": f"Age at donation: {a}",
			"LE Donor (yr)": f"{v['le_donor']:.2f}",
			"LE Non-donor (yr)": f"{v['le_nondonor']:.2f}",
			"ΔLE (days)": f"{v['diff']*365.25:+.1f}",
			"Note": ""
		})

	for r, v in race_res.items():
		rows.append({
			"Analysis": f"Race: {r}",
			"LE Donor (yr)": f"{v['le_donor']:.2f}",
			"LE Non-donor (yr)": f"{v['le_nondonor']:.2f}",
			"ΔLE (days)": f"{v['diff']*365.25:+.1f}",
			"Note": "age 40"
		})

	if sex_res:
		for s, v in sex_res.items():
			rows.append({
				"Analysis": f"Sex: {s}",
				"LE Donor (yr)": f"{v['le_donor']:.2f}",
				"LE Non-donor (yr)": f"{v['le_nondonor']:.2f}",
				"ΔLE (days)": f"{v['diff']*365.25:+.1f}",
				"Note": "age 40"
			})

	for lbl in ["No priority (PLD=standard)", "Donor mort HR=1.30 (Mjøen)",
				"ESRD RR ×11 (Mjøen upper)", "PLD wait 50 days (optimistic)",
				"Post-Tx: LDKT quality"]:
		rows.append({
			"Analysis": f"Sensitivity: {lbl}",
			"LE Donor (yr)": "—",
			"LE Non-donor (yr)": "—",
			"ΔLE (days)": f"{owsa_res[lbl]*365.25:+.1f}",
			"Note": "OWSA"
		})

	return pd.DataFrame(rows)


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main():
	print("="*60)
	print("LIVING KIDNEY DONATION MARKOV MODEL")
	print("="*60)

	base_res  = run_base_case(age_at_donation=40, n=N_SIM)
	psa_diffs = run_psa(age_at_donation=40, n_draws=N_DRAWS)
	owsa_res  = run_owsa(age_at_donation=40)
	age_res   = run_age_subgroups()
	race_res  = run_race_subgroups(age_at_donation=40)
	sex_res   = run_sex_subgroups(age_at_donation=40)

	# Save individual panel figures
	panels = [
		(make_fig_distributions, "kidney_model_A_distributions.png", (base_res,)),
		(make_fig_psa,           "kidney_model_B_psa.png",           (psa_diffs,)),
		(make_fig_tornado,       "kidney_model_C_tornado.png",       (owsa_res,)),
		(make_fig_age_subgroup,  "kidney_model_D_age_subgroup.png",  (age_res,)),
		(make_fig_race_subgroup, "kidney_model_E_race_subgroup.png", (race_res,)),
	]
	for make_fn, fname, args in panels:
		fig = make_fn(*args)
		path = RESULTS / fname
		fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
		plt.close()
		print(f"Figure saved: {path}")

	# Save age × race matrix figure
	matrix_fig = make_age_race_matrix()
	matrix_path = RESULTS / "kidney_model_age_race_matrix.png"
	matrix_fig.savefig(matrix_path, dpi=150, bbox_inches="tight",
					   facecolor=matrix_fig.get_facecolor())
	plt.close()
	print(f"Age×race matrix saved: {matrix_path}")

	# Save age × sex matrix figure
	age_sex_fig  = make_age_sex_matrix()
	age_sex_path = RESULTS / "kidney_model_age_sex_matrix.png"
	age_sex_fig.savefig(age_sex_path, dpi=150, bbox_inches="tight",
						facecolor=age_sex_fig.get_facecolor())
	plt.close()
	print(f"Age×sex matrix saved: {age_sex_path}")

	# Save sex × race matrix figure
	sex_race_fig  = make_sex_race_matrix()
	sex_race_path = RESULTS / "kidney_model_sex_race_matrix.png"
	sex_race_fig.savefig(sex_race_path, dpi=150, bbox_inches="tight",
						 facecolor=sex_race_fig.get_facecolor())
	plt.close()
	print(f"Sex×race matrix saved: {sex_race_path}")

	# Save age × race × sex faceted matrix figure
	ars_fig  = make_age_race_sex_matrix()
	ars_path = RESULTS / "kidney_model_age_race_sex_matrix.png"
	ars_fig.savefig(ars_path, dpi=150, bbox_inches="tight",
					facecolor=ars_fig.get_facecolor())
	plt.close()
	print(f"Age×race×sex matrix saved: {ars_path}")

	# Save race figure
	race_fig = make_race_figure(race_res)
	race_fig_path = RESULTS / "kidney_model_race.png"
	race_fig.savefig(race_fig_path, dpi=150, bbox_inches="tight",
					 facecolor=race_fig.get_facecolor())
	plt.close()
	print(f"Race figure saved: {race_fig_path}")

	# Save sex figure
	sex_fig = make_sex_figure(sex_res)
	sex_fig_path = RESULTS / "kidney_model_sex.png"
	sex_fig.savefig(sex_fig_path, dpi=150, bbox_inches="tight",
					facecolor=sex_fig.get_facecolor())
	plt.close()
	print(f"Sex figure saved: {sex_fig_path}")

	# Save results table
	df = make_results_table(base_res, psa_diffs, owsa_res, age_res, race_res, sex_res)
	csv_path = RESULTS / "kidney_model_results.csv"
	df.to_csv(csv_path, index=False)
	print(f"Results table saved: {csv_path}")
	print("\n" + df.to_string(index=False))
	print("\nDone.")


if __name__ == "__main__":
	main()
