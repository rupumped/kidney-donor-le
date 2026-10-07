"""
tests/test_model.py
────────────────────
Validation checks for the Markov simulation.

Run with: python -m pytest tests/ -v
Or:        python tests/test_model.py
"""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from utils import (load_life_table, _hardcoded_base_params,
                   weibull_scale_from_cumrisk, weibull_scale_from_cumrisk_competing,
                   weibull_annual_prob, median_to_annual_tx_prob,
                   dialysis_annual_mort, posttx_annual_mort, waitlist_annual_mort_at,
                   esrd_onset_mortality_curve, grams_nondonor_age_scale)

# Import model functions without running main
import importlib.util, types

def load_script(filename):
    spec = importlib.util.spec_from_file_location(
        filename.split("_")[0], Path(__file__).resolve().parent.parent / "src" / filename
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

SIM = load_script("06_markov_simulation.py")
BASE = _hardcoded_base_params()


# ────────────────────────────────────────────────────────────────────────────
class TestLifeTable:
    def test_length(self):
        qx = load_life_table()
        assert len(qx) == 101, "Life table should cover ages 0–100"

    def test_range(self):
        qx = load_life_table()
        assert (qx >= 0).all() and (qx <= 1).all(), "qx must be in [0,1]"

    def test_monotone_at_old_ages(self):
        qx = load_life_table()
        # qx should be increasing for ages 40–99
        assert all(qx[i+1] >= qx[i] for i in range(40, 98)), \
            "qx should be non-decreasing at old ages"


# ────────────────────────────────────────────────────────────────────────────
class TestWeibullCalibration:
    """
    Annual transition probabilities must sum over 15 yrs to exactly reproduce
    Muzaale 2014 CIF (the sum telescopes to 1 - S(15) by construction).
    """

    def _integrate_cif(self, cum_risk_15, k, years=15):
        lam = weibull_scale_from_cumrisk(cum_risk_15, k)
        # Sum of annual probabilities telescopes exactly to 1 - S(years)
        return sum(weibull_annual_prob(t, lam, k) for t in range(years))

    def test_donor_overall_calibration(self):
        target = BASE["esrd_15yr_donor_overall"]
        achieved = self._integrate_cif(target, BASE["weibull_shape"])
        assert abs(achieved - target) / target < 0.001, \
            f"Donor calibration off: target={target:.4%}, achieved={achieved:.4%}"

    def test_nondonor_calibration(self):
        target = BASE["esrd_15yr_nondonor"]
        achieved = self._integrate_cif(target, BASE["weibull_shape"])
        assert abs(achieved - target) / target < 0.001, \
            f"Non-donor calibration off: target={target:.4%}, achieved={achieved:.4%}"

    def test_accelerating_probability(self):
        """Weibull with k>1 should give higher annual probability at t=10 than t=1."""
        k = BASE["weibull_shape"]
        lam = weibull_scale_from_cumrisk(BASE["esrd_15yr_donor_overall"], k)
        p1  = weibull_annual_prob(1, lam, k)
        p10 = weibull_annual_prob(10, lam, k)
        assert p10 > p1, f"Annual probability should increase over time (k={k}>1)"

    def test_zero_before_donation(self):
        k = BASE["weibull_shape"]
        lam = weibull_scale_from_cumrisk(BASE["esrd_15yr_donor_overall"], k)
        assert weibull_annual_prob(-1, lam, k) == 0.0, "Probability must be 0 before donation"

    def test_donor_higher_than_nondonor(self):
        """Donor annual probability should be ~8× non-donor at any given time point."""
        k = BASE["weibull_shape"]
        lam_d  = weibull_scale_from_cumrisk(BASE["esrd_15yr_donor_overall"], k)
        lam_nd = weibull_scale_from_cumrisk(BASE["esrd_15yr_nondonor"], k)
        p_d  = weibull_annual_prob(5, lam_d, k)
        p_nd = weibull_annual_prob(5, lam_nd, k)
        rr = p_d / p_nd
        assert 6 < rr < 11, f"Donor:nondonor probability ratio at t=5 should be ~8×, got {rr:.1f}×"

    def test_competing_risk_calibration(self):
        """weibull_scale_from_cumrisk_competing must reproduce the target competing-risk CIF.

        This is the function actually used in the simulation; the simpler
        weibull_scale_from_cumrisk (tested above) is not used at runtime.
        """
        qx     = load_life_table()
        target = BASE["esrd_15yr_donor_overall"]
        k      = BASE["weibull_shape"]
        age0   = 40
        lam    = weibull_scale_from_cumrisk_competing(target, k, qx, age0, bg_hr=1.0)
        # Replicate the competing-risk integration from weibull_scale_from_cumrisk_competing
        cif, s_bg = 0.0, 1.0
        for t in range(15):
            cif  += s_bg * weibull_annual_prob(t, lam, k)
            s_bg *= 1.0 - qx[min(age0 + t, 100)]
        assert abs(cif - target) / target < 0.001, (
            f"Competing-risk calibration off: target={target:.4%}, achieved={cif:.4%}"
        )

    def test_competing_risk_bg_hr_effect(self):
        """Higher bg_hr requires a smaller Weibull lambda to hit the same CIF.

        Higher background mortality depletes the at-risk pool faster (stronger
        competing risk), which lowers the observed ESRD CIF.  To compensate and
        still match the target CIF, the intrinsic ESRD hazard must be higher —
        and in the Weibull parameterisation S(t)=exp(-(t/λ)^k), smaller λ means
        higher hazard, so λ must decrease.
        """
        qx     = load_life_table()
        target = BASE["esrd_15yr_donor_overall"]
        k      = BASE["weibull_shape"]
        lam_hr1  = weibull_scale_from_cumrisk_competing(target, k, qx, 40, bg_hr=1.0)
        lam_hr13 = weibull_scale_from_cumrisk_competing(target, k, qx, 40, bg_hr=1.30)
        assert lam_hr13 < lam_hr1, (
            f"Higher bg_hr needs smaller lambda (stronger hazard); "
            f"got lam_hr1={lam_hr1:.2f}, lam_hr13={lam_hr13:.2f}"
        )


# ────────────────────────────────────────────────────────────────────────────
class TestWaitlistTransitions:
    def test_priority_faster_than_standard(self):
        p_pri = median_to_annual_tx_prob(100)
        p_std = median_to_annual_tx_prob(985)
        assert p_pri > p_std, "Priority annual Tx prob should exceed standard"

    def test_priority_near_certain_in_year1(self):
        """100-day median → ~92%/yr annual Tx probability."""
        p_pri = median_to_annual_tx_prob(100)
        assert 0.85 < p_pri < 0.99, f"Priority Tx prob should be ~0.92, got {p_pri:.3f}"

    def test_standard_about_22pct_per_year(self):
        """985-day (~32.8 month) median → ~22%/yr annual probability."""
        p_std = median_to_annual_tx_prob(985)
        assert 0.18 < p_std < 0.28, f"Standard Tx prob should be ~0.22, got {p_std:.3f}"

    def test_median_consistency(self):
        """Verify: simulating the exponential model should give ~median at 50th pct."""
        rng = np.random.default_rng(0)
        n = 100_000
        annual_p = median_to_annual_tx_prob(985)
        wait_years = []
        for _ in range(n):
            for yr in range(1, 50):
                if rng.random() < annual_p:
                    wait_years.append(yr)
                    break
        median_sim = np.median(wait_years) * 365.25
        # Should be within 15% of 985 days
        assert abs(median_sim - 985) / 985 < 0.15, \
            f"Simulated median {median_sim:.0f}d vs target 985d"


# ────────────────────────────────────────────────────────────────────────────
class TestSimulation:
    """Smoke tests on the full simulation (small n for speed)."""

    def test_healthy_arm_le_reasonable(self):
        """LE from age 40 should be 30–40 years."""
        p = BASE.copy()
        rng = np.random.default_rng(0)
        ly = SIM.simulate_cohort(p, 10_000, 40, donor=False, rng=rng)
        le = ly.mean()
        assert 28 < le < 42, f"LE from age 40 should be ~33 yr, got {le:.1f}"

    def test_donor_le_lower_than_nondonor(self):
        """Under the base-case time-varying donor mortality HR (flat at 1.0
        through year 10, ramping to 1.30 by year 15, held thereafter), donor
        LE should be lower than non-donor LE by roughly 2-3 years -- this is
        now the dominant, intentional effect in the model (see design.tex),
        not a small ESRD-driven perturbation, so a tight "should be close"
        tolerance no longer applies."""
        p = BASE.copy()
        le_nd = SIM.simulate_cohort(p, 20_000, 40, False, rng=np.random.default_rng(1)).mean()
        le_d  = SIM.simulate_cohort(p, 20_000, 40, True,  rng=np.random.default_rng(1)).mean()
        diff = le_nd - le_d
        assert 1.5 < diff < 4.0, \
            f"Donor LE should be ~2-3 yr lower than non-donor under base-case HR, diff={diff:.2f}"

    def test_esrd_cost_negative(self):
        """Donation without priority should reduce LE vs non-donation."""
        p = BASE.copy()
        p["wl_pld_mean_days"] = p["wl_std_mean_days"]  # strip priority
        le_nd = SIM.simulate_cohort(p, 50_000, 40, False, rng=np.random.default_rng(2)).mean()
        le_d  = SIM.simulate_cohort(p, 50_000, 40, True,  rng=np.random.default_rng(2)).mean()
        assert le_d < le_nd, "Donation without priority should reduce LE"

    def test_priority_helps_esrd_patients(self):
        """Priority wait should give higher LE than standard wait for ESRD patients."""
        p_pri = BASE.copy()
        p_std = BASE.copy()
        p_std["wl_pld_mean_days"] = BASE["wl_std_mean_days"]
        # Start cohort in ESRD state directly
        def sim_from_esrd(params, priority, seed):
            state = np.ones(20_000, dtype=np.int8)
            ly    = np.zeros(20_000)
            alive = np.ones(20_000, dtype=bool)
            rng   = np.random.default_rng(seed)
            wl_tx = SIM.waitlist_annual_tx_prob(params, priority)
            wl_m  = SIM.waitlist_annual_mort(params)
            d_m   = dialysis_annual_mort(params, 50)
            pt_m  = posttx_annual_mort(params, 50)
            for yr in range(50):
                a = 50 + yr
                if a >= 100: break
                u = rng.random((20_000, 3))
                ns = state.copy()
                m1 = (state==1)&alive
                if m1.any():
                    idx=np.where(m1)[0]; die=u[m1,0]<d_m
                    ns[idx[die]]=5; ns[idx[~die]] = 3 if priority else 2
                wl_s = 3 if priority else 2
                mwl = (state==wl_s)&alive
                if mwl.any():
                    idx=np.where(mwl)[0]; dw=u[mwl,0]<wl_m; gt=~dw&(u[mwl,1]<wl_tx)
                    ns[idx[dw]]=5; ns[idx[gt]]=4
                m4=(state==4)&alive
                if m4.any():
                    idx=np.where(m4)[0]; dp=u[m4,0]<pt_m; ns[idx[dp]]=5
                surv=ns!=5; ly[surv&alive]+=1; alive=surv; state=ns
            return ly.mean()

        le_pri = sim_from_esrd(BASE, True,  10)
        le_std = sim_from_esrd(BASE, False, 10)
        assert le_pri > le_std, \
            f"Priority should give higher LE for ESRD patients: {le_pri:.2f} vs {le_std:.2f}"

    def test_older_donation_less_harmful(self):
        """Net harm from donation should be smaller at older donation ages."""
        p = BASE.copy()
        diffs = {}
        for age in [25, 40, 55]:
            le_nd = SIM.simulate_cohort(p, 30_000, age, False, rng=np.random.default_rng(age)).mean()
            le_d  = SIM.simulate_cohort(p, 30_000, age, True,  rng=np.random.default_rng(age)).mean()
            diffs[age] = le_d - le_nd
        assert diffs[25] < diffs[55], \
            f"Harm should be smaller at age 55 than 25: {diffs[25]:.3f} vs {diffs[55]:.3f}"


# ────────────────────────────────────────────────────────────────────────────
class TestNonDonorAgeGradient:
    """Non-donor ESRD baseline by age follows Grams 2016 projections."""

    def test_reproduces_grams_cells(self):
        for race, sex, a, expected in [("Black", "Male", 20, 0.0008), ("Black", "Male", 60, 0.0032),
                                       ("White", "Female", 20, 0.0001), ("White", "Female", 60, 0.0008)]:
            at40 = BASE[f"grams_nondonor_15yr_{race.lower()}_{sex.lower()}_age40"]
            got = at40 * grams_nondonor_age_scale(BASE, a, race, sex)
            assert abs(got - expected) < 1e-12, f"{race} {sex} age {a}: {got} vs {expected}"

    def test_black_gradient_flatter_after_40(self):
        black = grams_nondonor_age_scale(BASE, 60, "Black")
        white = grams_nondonor_age_scale(BASE, 60, "White")
        assert black < white, f"Grams: Black 40→60 rise ({black:.2f}×) < White ({white:.2f}×)"


# ────────────────────────────────────────────────────────────────────────────
class TestPostESRDMortality:
    """Age-dependent post-ESRD mortality and the life-table floor."""

    def test_never_below_background(self):
        """Post-ESRD state mortality must never fall below general-population
        qx at the same age, or ESRD becomes protective at old ages."""
        qx = load_life_table()
        for age in range(18, 101):
            for donor in (False, True):
                assert dialysis_annual_mort(BASE, age, donor=donor, q_floor=qx[age]) >= qx[age]
                assert posttx_annual_mort(BASE, age, donor=donor, q_floor=qx[age]) >= qx[age]
            assert waitlist_annual_mort_at(BASE, age, q_floor=qx[age]) >= qx[age]

    def test_dialysis_mortality_rises_with_age(self):
        ages = [30, 50, 70, 80]
        q = [dialysis_annual_mort(BASE, a) for a in ages]
        assert q == sorted(q) and q[0] < q[-1], f"Dialysis mortality should rise with age: {q}"

    def test_no_jumps_at_band_edges(self):
        """Interpolated age curves must not step at registry band edges
        (steps produced a spurious bump in the ESRD-conditional age sweep)."""
        for f in (dialysis_annual_mort, posttx_annual_mort, waitlist_annual_mort_at):
            for age in range(20, 95):
                a, b = f(BASE, age), f(BASE, age + 1)
                assert abs(b / a - 1) < 0.10, f"{f.__name__} jumps {a:.4f}→{b:.4f} at {age}→{age+1}"

    def test_donor_dialysis_mortality_lower(self):
        """Muzaale 2016: prior donors with ESRD have lower mortality than
        general ESRD patients of the same age."""
        assert dialysis_annual_mort(BASE, 50, donor=True) < dialysis_annual_mort(BASE, 50)

    def test_more_esrd_never_helps(self):
        """With no excess all-cause mortality, raising donor ESRD risk must
        lower donor LE. Regression test: without the background floor, even
        donor-specific (lower) ESRD mortality made ESRD net-protective."""
        # Worst case for the paradox: older entry, every ESRD patient listed
        # preemptively and transplanted fast, so most post-ESRD time is spent
        # post-Tx at old ages where unfloored post-Tx mortality < qx.
        p = dict(BASE, donor_mort_hr_late=1.0, donor_esrd_preemptive_prob=1.0)
        le = [SIM.run_arm_analytic(dict(p, esrd_15yr_donor_overall=r), 70, donor=True)
              for r in (0.001, 0.01, 0.05)]
        assert le[0] > le[1] > le[2], f"LE should fall as ESRD risk rises: {le}"

    def test_curve_matches_esrd_conditional_model(self):
        """utils.esrd_onset_mortality_curve (used for calibration) must match
        script 10's ESRD-conditional cohort for the general pathway."""
        m10 = load_script("10_esrd_conditional_cohort_markov.py")
        p = dict(m10.BASE)
        n = 1e6   # run_arm stops below 0.5 survivors, so n must be large
        _, trace, _ = m10.run_arm(p, n, 50, priority=False)
        curve = esrd_onset_mortality_curve(p, 50, False, [1, 3, 5, 10], m10.LIFE_TABLE_QX)
        for t in (1, 3, 5, 10):
            mort10 = 1 - trace[t - 1]["alive"] / n
            assert abs(mort10 - curve[t]) < 1e-9, \
                f"Mortality at {t} yr differs: script 10 {mort10:.6f} vs utils {curve[t]:.6f}"


# ────────────────────────────────────────────────────────────────────────────
def run_all_tests():
    """Run all tests without pytest."""
    import traceback
    classes = [TestLifeTable, TestWeibullCalibration,
               TestWaitlistTransitions, TestSimulation, TestPostESRDMortality,
               TestNonDonorAgeGradient]
    passed = failed = 0
    for cls in classes:
        obj = cls()
        methods = [m for m in dir(obj) if m.startswith("test_")]
        for meth in methods:
            try:
                getattr(obj, meth)()
                print(f"  PASS  {cls.__name__}.{meth}")
                passed += 1
            except AssertionError as e:
                print(f"  FAIL  {cls.__name__}.{meth}: {e}")
                failed += 1
            except Exception as e:
                print(f"  ERROR {cls.__name__}.{meth}: {e}")
                traceback.print_exc()
                failed += 1
    print(f"\n{passed} passed, {failed} failed")
    return failed == 0


if __name__ == "__main__":
    print("Running model validation tests...\n")
    ok = run_all_tests()
    sys.exit(0 if ok else 1)
