"""
Test whether a weighted-inputs "Overall" row sits between the White and Black
rows in the age × race × sex matrix.

Weights: SRTR 2023 ADR Fig KI 105 (US living kidney donors, 2023)
  f_black = 0.072; all non-Black treated as White-equivalent for ESRD.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))
spec = importlib.util.spec_from_file_location("mks", SRC / "06_markov_simulation.py")
mks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mks)

BASE = mks.BASE
_scale_posttx_mort     = mks._scale_posttx_mort
_age_adjust_esrd_nonblack = mks._age_adjust_esrd_nonblack
_age_adjust_esrd_nondonor = mks._age_adjust_esrd_nondonor
_GRAMS_NONDONOR        = mks._GRAMS_NONDONOR
run_arm_analytic       = mks.run_arm_analytic
LIFE_TABLE_QX          = mks.LIFE_TABLE_QX
load_life_table        = mks.load_life_table
DATA_PROC              = mks.DATA_PROC

F_BLACK = 0.072   # SRTR 2023 ADR KI 105

ages = [25, 30, 35, 40, 45, 50, 55]

lt_m = load_life_table(DATA_PROC / "lifetable_male_2021.csv")
lt_f = load_life_table(DATA_PROC / "lifetable_female_2021.csv")

f_female = 0.635
hr_male  = float(BASE.get("hr_male_sex", 1.88))
denom    = f_female + (1.0 - f_female) * hr_male

white = {
    "esrd_15yr_donor_overall": 0.00227,
    "esrd_15yr_nondonor":      0.00050,
    "wl_mort_race_mult": BASE["wl_mort_race_mult_white"],
    "posttx_mort_target": BASE.get("posttx_annual_mort_white", 0.038),
}
black = {
    "esrd_15yr_donor_overall": 0.00747,
    "esrd_15yr_nondonor":      0.00195,
    "wl_mort_race_mult": BASE["wl_mort_race_mult_black"],
    "posttx_mort_target": BASE.get("posttx_annual_mort_black", 0.035),
}

def weighted_overall(w_black: float) -> dict:
    wb, ww = w_black, 1.0 - w_black
    return {
        "esrd_15yr_donor_overall": ww * white["esrd_15yr_donor_overall"] + wb * black["esrd_15yr_donor_overall"],
        "esrd_15yr_nondonor":      ww * white["esrd_15yr_nondonor"]      + wb * black["esrd_15yr_nondonor"],
        "wl_mort_race_mult":       ww * white["wl_mort_race_mult"]       + wb * black["wl_mort_race_mult"],
        "posttx_mort_target":      ww * white["posttx_mort_target"]      + wb * black["posttx_mort_target"],
    }

race_defs = {
    "White":   {**white,                      "race_label": "White"},
    "Overall": {**weighted_overall(F_BLACK), "race_label": "Overall"},
    "Black":   {**black,                      "race_label": "Black"},
}

sex_specs = {
    "Female":  {"scale": 1.0 / denom,       "_lt": lt_f},
    "Overall": {"scale": 1.0,                "_lt": LIFE_TABLE_QX},
    "Male":    {"scale": hr_male / denom,    "_lt": lt_m},
}

def grams_for(race_label: str, sex: str, w_black: float, fallback: float) -> float:
    """For the weighted Overall race, interpolate the Grams non-donor baseline
    the same way as the inputs. Falls back to `fallback` (the race's own
    per-sex-scaled non-donor baseline) when no (race, sex) entry exists —
    matching the original make_age_race_sex_matrix behavior."""
    if race_label == "Overall":
        ww, wb = 1.0 - w_black, w_black
        w_val = _GRAMS_NONDONOR.get(("White", sex), white["esrd_15yr_nondonor"])
        b_val = _GRAMS_NONDONOR.get(("Black", sex), black["esrd_15yr_nondonor"])
        return ww * w_val + wb * b_val
    return _GRAMS_NONDONOR.get((race_label, sex), fallback)

def simulate(races, sex_specs, w_black):
    mats = {}
    for sex, ss in sex_specs.items():
        lt = ss["_lt"]
        mat = np.zeros((len(races), len(ages)))
        for r, (race_name, rdef) in enumerate(races.items()):
            p = BASE.copy()
            p["esrd_15yr_donor_overall"] = rdef["esrd_15yr_donor_overall"]
            p["esrd_15yr_nondonor"]      = rdef["esrd_15yr_nondonor"]
            p["wl_mort_race_mult"]       = rdef["wl_mort_race_mult"]
            p.update(_scale_posttx_mort(BASE, rdef["posttx_mort_target"]))
            base_donor    = p["esrd_15yr_donor_overall"] * ss["scale"]
            base_nondonor = p["esrd_15yr_nondonor"]      * ss["scale"]
            for a, age in enumerate(ages):
                # Use the non-Black age HR path for White/Overall; direct for Black
                if race_name != "Black":
                    p["esrd_15yr_donor_overall"] = _age_adjust_esrd_nonblack(base_donor, age)
                else:
                    p["esrd_15yr_donor_overall"] = base_donor
                nondonor_ref = grams_for(rdef["race_label"], sex, w_black, base_nondonor) \
                               if race_name == "Overall" \
                               else _GRAMS_NONDONOR.get((race_name, sex), base_nondonor)
                p["esrd_15yr_nondonor"] = _age_adjust_esrd_nondonor(
                    nondonor_ref, age, race=race_name, sex=sex
                )
                le_d  = run_arm_analytic(p, age, donor=True,  life_table=lt)
                le_nd = run_arm_analytic(p, age, donor=False, life_table=lt)
                mat[r, a] = (le_d - le_nd) * 365.25
        mats[sex] = mat
    return mats

def report(mats, label):
    print(f"\n── {label} ─────────────────────────────────────────────────────")
    header = "race   " + "  ".join(f"A{a:>2}" for a in ages)
    out_of_range = 0
    for sex, mat in mats.items():
        print(f"\n{sex}")
        print(header)
        for r, race in enumerate(race_defs.keys()):
            cells = "  ".join(f"{mat[r,a]:+5.0f}" for a in range(len(ages)))
            print(f"{race:<6} {cells}")
        # Check Overall between White and Black per cell
        w, o, b = mat[0], mat[1], mat[2]
        for a in range(len(ages)):
            lo, hi = min(w[a], b[a]), max(w[a], b[a])
            if not (lo - 1e-9 <= o[a] <= hi + 1e-9):
                out_of_range += 1
                print(f"  OUT-OF-RANGE at age {ages[a]}: W={w[a]:+.1f} O={o[a]:+.1f} B={b[a]:+.1f}")
    print(f"\nTotal cells where Overall is NOT between White and Black: {out_of_range}")
    return out_of_range

# Main test at f_black = 0.072
mats = simulate(race_defs, sex_specs, F_BLACK)
total_out = report(mats, f"Weighted inputs, f_black = {F_BLACK:.3f}")

# Sweep
print("\n\n════ Sweep over f_black ════")
for fb in [0.00, 0.05, 0.072, 0.10, 0.15, 0.20, 0.50, 1.00]:
    rd = {
        "White":   {**white,                 "race_label": "White"},
        "Overall": {**weighted_overall(fb), "race_label": "Overall"},
        "Black":   {**black,                 "race_label": "Black"},
    }
    ms = simulate(rd, sex_specs, fb)
    oor = 0
    for sex, mat in ms.items():
        w, o, b = mat[0], mat[1], mat[2]
        for a in range(len(ages)):
            lo, hi = min(w[a], b[a]), max(w[a], b[a])
            if not (lo - 1e-9 <= o[a] <= hi + 1e-9):
                oor += 1
    print(f"  f_black={fb:.3f}: cells out-of-range = {oor}/{3*len(ages)}")
