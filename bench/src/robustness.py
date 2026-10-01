"""
Robustness harness.

Sweeps every registered free parameter across its range and reports which
conclusions survive the whole space and which are artefacts of one setting.

A conclusion that fails anywhere in the space may not be reported as a number.
This exists because two conclusions flipped in a single day when one unsourced
constant was corrected.

Only the latent physiology and the reference channel are needed for the oracle
floors, so the sensor model and the candidate estimators are skipped: the floors
are what the BEST estimator could reach from that information, not what mine does.
"""
import os, sys, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import tbgen as tb, params as pr

N_SUB = 12
SEED = 20260918
REPORTED = {"tp_cycle": 0.795, "tp_posture": 0.799}

# parameters that reach the latent drivers or the reference channel
PHYS_KEYS = ["co2_react_hypo_ratio", "autoreg_gain", "mcav_supine_extra",
             "mcav_co_coupling", "hydrostatic", "supine_sv_rise", "orthostatic_hr",
             "map_exercise_gain", "co_exercise_gain", "co2_rise_slope",
             "hypervent_threshold", "hypervent_slope", "skin_supine_rise",
             "skin_co_coupling", "skin_thermal_gain", "tcd_noise", "tcd_drift_sd",
             "tcd_dropout", "w_max_mean"]


def _z(x):
    x = np.asarray(x, float); sd = np.nanstd(x)
    return (x - np.nanmean(x)) / (sd if sd > 0 else 1.0)


def _r(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    k = np.isfinite(a) & np.isfinite(b)
    if k.sum() < 30 or np.std(a[k]) == 0 or np.std(b[k]) == 0:
        return np.nan
    return float(np.corrcoef(a[k], b[k])[0, 1])


def measure(w_max_mean):
    """Oracle floors for both protocols under the current OVERRIDE."""
    master = np.random.default_rng(SEED)
    subs = [tb.draw_subject(i, master) for i in range(N_SUB)]
    wm = {s["sid"]: float(np.clip(master.normal(w_max_mean, 0.08), 0.30, 0.80)) for s in subs}
    out = {}
    for pname, cfg in tb.TEMPLE_PROTOCOLS.items():
        acc = {k: [] for k in ("ceiling", "scalp", "hr", "co")}
        for sub in subs:
            ph = tb.physiology(sub, cfg["kind"], cfg["dur"],
                               np.random.default_rng(sub["seed"] + 101),
                               w_max=wm[sub["sid"]])
            ref = tb.references(sub, ph, np.random.default_rng(sub["seed"] + 7717))
            n = len(ref["tcd"]); sl = (np.arange(n) * int(tb.FS)).astype(int)
            tz = _z(ref["tcd"])
            a_ec = ph["skin"] * (1 + 0.18 * (ph["hr_trend"] / sub["hr_rest"] - 1))
            acc["ceiling"].append(_r(_z(ph["mcav"][sl]), tz))
            acc["scalp"].append(_r(_z(a_ec[sl]), tz))
            acc["hr"].append(_r(_z(ref["ecg_hr"]), tz))
            acc["co"].append(_r(_z(ph["co"][sl]), tz))
        out[pname] = {k: float(np.nanmedian(v)) for k, v in acc.items()}
    return out


CLAIMS = {
  "E1 exercise: scalp alone reaches the reported 0.795":
      lambda m: abs(m["tp_cycle"]["scalp"]) >= REPORTED["tp_cycle"],
  "E2 exercise: heart rate alone reaches the reported 0.795":
      lambda m: abs(m["tp_cycle"]["hr"]) >= REPORTED["tp_cycle"],
  "P1 posture: scalp alone reaches the reported 0.799":
      lambda m: abs(m["tp_posture"]["scalp"]) >= REPORTED["tp_posture"],
  "P2 posture: cardiac output alone reaches the reported 0.799":
      lambda m: abs(m["tp_posture"]["co"]) >= REPORTED["tp_posture"],
  "P3 posture: the cardiac-output confound is POSITIVE (tracks MCAv)":
      lambda m: m["tp_posture"]["co"] > 0,
  "P4 posture: heart rate alone is NEGATIVE (the paper's own argument)":
      lambda m: m["tp_posture"]["hr"] < 0,
}


def main():
    settings = [("baseline", {})]
    for k in PHYS_KEYS:
        spec = pr.P[k]
        for tag, val in (("lo", spec["lo"]), ("hi", spec["hi"])):
            settings.append((f"{k}={val:g} ({tag})", {k: val}))

    results, held = [], {c: [] for c in CLAIMS}
    for label, ov in settings:
        tb.OVERRIDE = {k: v for k, v in ov.items() if k != "w_max_mean"}
        wmm = ov.get("w_max_mean", pr.P["w_max_mean"]["v"])
        m = measure(wmm)
        verdicts = {c: bool(f(m)) for c, f in CLAIMS.items()}
        for c, v in verdicts.items():
            held[c].append((label, v))
        results.append((label, m, verdicts))
    tb.OVERRIDE = {}

    print(f"{len(settings)} settings x {N_SUB} subjects x 2 protocols\n")
    print("CLAIM ROBUSTNESS")
    for c, rows in held.items():
        fails = [lab for lab, v in rows if not v]
        n_ok = len(rows) - len(fails)
        status = "SURVIVES ALL" if not fails else f"FAILS in {len(fails)}"
        print(f"\n  [{n_ok}/{len(rows)}] {c}\n      {status}")
        for lab in fails[:6]:
            print(f"        broken by: {lab}")
        if len(fails) > 6:
            print(f"        ... and {len(fails)-6} more")

    base = dict(results)[("baseline")] if False else results[0][1]
    print("\n\nBASELINE VALUES (all parameters at registered defaults)")
    for pname in ("tp_cycle", "tp_posture"):
        m = base[pname]
        lab = "exercise" if pname == "tp_cycle" else "posture "
        print(f"  {lab}  ceiling {m['ceiling']:+.3f}   scalp {m['scalp']:+.3f}   "
              f"HR {m['hr']:+.3f}   CO {m['co']:+.3f}   reported {REPORTED[pname]:+.3f}")

    # spread of each floor across the whole space
    print("\nRANGE OF EACH FLOOR ACROSS THE SWEPT SPACE")
    for pname in ("tp_cycle", "tp_posture"):
        for k in ("scalp", "co"):
            vals = [abs(m[pname][k]) for _, m, _ in results]
            lab = "exercise" if pname == "tp_cycle" else "posture "
            print(f"  {lab} {k:6s} |r| from {min(vals):.3f} to {max(vals):.3f}")
    return results


if __name__ == "__main__":
    main()
