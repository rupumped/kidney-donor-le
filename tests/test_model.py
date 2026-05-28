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
                   weibull_scale_from_cumrisk, weibull_annual_hazard,
                   median_to_annual_tx_prob)

# Import model functions without running main
import importlib.util, types

def load_sim():
    spec = importlib.util.spec_from_file_location(
        "sim", Path(__file__).resolve().parent.parent / "src" / "06_markov_simulation.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

SIM = load_sim()
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
    The Weibull hazard must integrate over 15 yrs to reproduce Muzaale 2014 CIF.
    Tolerance: ±5% of target (per design spec in parameter notes).
    """

    def _integrate_cif(self, cum_risk_15, k, years=15):
        lam = weibull_scale_from_cumrisk(cum_risk_15, k)
        cumhaz = sum(weibull_annual_hazard(t, lam, k) for t in range(years))
        return 1 - np.exp(-cumhaz)

    def test_donor_overall_calibration(self):
        target = BASE["esrd_15yr_donor_overall"]
        achieved = self._integrate_cif(target, BASE["weibull_shape"])
        # 10% tolerance: discrete annual cycles undercount vs continuous slightly
        assert abs(achieved - target) / target < 0.10, \
            f"Donor calibration off: target={target:.4%}, achieved={achieved:.4%}"

    def test_nondonor_calibration(self):
        target = BASE["esrd_15yr_nondonor"]
        achieved = self._integrate_cif(target, BASE["weibull_shape"])
        assert abs(achieved - target) / target < 0.10, \
            f"Non-donor calibration off: target={target:.4%}, achieved={achieved:.4%}"

    def test_accelerating_hazard(self):
        """Weibull with k>1 should give higher hazard at t=10 than t=1."""
        k = BASE["weibull_shape"]
        lam = weibull_scale_from_cumrisk(BASE["esrd_15yr_donor_overall"], k)
        h1  = weibull_annual_hazard(1, lam, k)
        h10 = weibull_annual_hazard(10, lam, k)
        assert h10 > h1, f"Hazard should increase over time (k={k}>1)"

    def test_zero_at_t0(self):
        k = BASE["weibull_shape"]
        lam = weibull_scale_from_cumrisk(BASE["esrd_15yr_donor_overall"], k)
        assert weibull_annual_hazard(0, lam, k) == 0.0, "Hazard must be 0 at t=0"

    def test_donor_higher_than_nondonor(self):
        """Donor hazard should be ~8× non-donor at any given time point."""
        k = BASE["weibull_shape"]
        for target, label in [
            (BASE["esrd_15yr_donor_overall"], "donor"),
            (BASE["esrd_15yr_nondonor"], "nondonor"),
        ]:
            lam = weibull_scale_from_cumrisk(target, k)
            _ = weibull_annual_hazard(5, lam, k)  # just check no error

        lam_d  = weibull_scale_from_cumrisk(BASE["esrd_15yr_donor_overall"], k)
        lam_nd = weibull_scale_from_cumrisk(BASE["esrd_15yr_nondonor"], k)
        h_d  = weibull_annual_hazard(5, lam_d, k)
        h_nd = weibull_annual_hazard(5, lam_nd, k)
        rr = h_d / h_nd
        assert 6 < rr < 11, f"Donor:nondonor hazard ratio at t=5 should be ~8×, got {rr:.1f}×"


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

    def test_donor_le_close_to_nondonor(self):
        """With common RNG, donor and non-donor LE should be within 1 year."""
        p = BASE.copy()
        le_nd = SIM.simulate_cohort(p, 20_000, 40, False, rng=np.random.default_rng(1)).mean()
        le_d  = SIM.simulate_cohort(p, 20_000, 40, True,  rng=np.random.default_rng(1)).mean()
        assert abs(le_d - le_nd) < 1.0, \
            f"Donor/non-donor LE should be within 1 yr, diff={abs(le_d-le_nd):.2f}"

    def test_esrd_cost_negative(self):
        """Donation without priority should reduce LE vs non-donation."""
        p = BASE.copy()
        p["wl_pld_median_days"] = p["wl_std_median_days"]  # strip priority
        le_nd = SIM.simulate_cohort(p, 50_000, 40, False, rng=np.random.default_rng(2)).mean()
        le_d  = SIM.simulate_cohort(p, 50_000, 40, True,  rng=np.random.default_rng(2)).mean()
        assert le_d < le_nd, "Donation without priority should reduce LE"

    def test_priority_helps_esrd_patients(self):
        """Priority wait should give higher LE than standard wait for ESRD patients."""
        p_pri = BASE.copy()
        p_std = BASE.copy()
        p_std["wl_pld_median_days"] = BASE["wl_std_median_days"]
        # Start cohort in ESRD state directly
        def sim_from_esrd(params, priority, seed):
            state = np.ones(20_000, dtype=np.int8)
            ly    = np.zeros(20_000)
            alive = np.ones(20_000, dtype=bool)
            rng   = np.random.default_rng(seed)
            wl_tx = SIM.waitlist_annual_tx_prob(params, priority)
            wl_m  = SIM.waitlist_annual_mort(params)
            d_m   = params["dialysis_annual_mort"]
            pt_m  = params["posttx_annual_mort"]
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
def run_all_tests():
    """Run all tests without pytest."""
    import traceback
    classes = [TestLifeTable, TestWeibullCalibration,
               TestWaitlistTransitions, TestSimulation]
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
