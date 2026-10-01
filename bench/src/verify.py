"""
Verification: does the cohort reproduce physiology that was NOT set as a target?

Each check is against a published range. A generator that only reproduces what it
was tuned to reproduce proves nothing, so these are effects that emerge from the
coupling constants rather than being written in directly. Any failure is loud.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
CHECKS, FAILS = [], 0


def check(name, value, lo, hi, source):
    global FAILS
    ok = np.isfinite(value) and lo <= value <= hi
    FAILS += (not ok)
    CHECKS.append(dict(check=name, value=round(float(value), 3),
                       expected=f"{lo} to {hi}", source=source,
                       verdict="PASS" if ok else "FAIL"))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:52s} {value:8.3f}   expected {lo} to {hi}")


def main():
    df = pd.read_csv(os.path.join(DATA, "streams_1hz.csv"))
    d = df[(df.cond == "real")]
    cyc = d[d.protocol == "cycle"]
    pos = d[d.protocol == "posture"]

    def stage(g, a, b, col="mcav_true"):
        return g[(g.t >= a) & (g.t < b)].groupby("sid")[col].mean()

    base = stage(cyc, 60, 170)
    print("\nCerebral haemodynamics during incremental exercise")
    for lab, (a, b, lo, hi) in {
        "MCAv rise at 30% VO2max (%)": (300, 355, 4, 20),
        "MCAv rise at 50% VO2max (%)": (520, 595, 10, 28),
        "MCAv at 70% VO2max, back toward baseline (%)": (760, 835, -12, 10),
    }.items():
        v = float(np.median(100 * (stage(cyc, a, b) / base - 1)))
        check(lab, v, lo, hi, "Ogoh & Ainslie 2009; Smith & Ainslie 2017")

    print("\nPostural dissociation: the condition that separates flow from heart rate")
    st_m, su_m = stage(pos, 150, 290), stage(pos, 400, 590)
    check("supine MCAv above standing (%)",
          float(np.median(100 * (su_m / st_m - 1))), 6, 20,
          "supine vs upright MCAv, hydrostatic + venous return")
    st_h = stage(pos, 150, 290, "hr_true"); su_h = stage(pos, 400, 590, "hr_true")
    check("supine heart rate below standing (bpm)",
          float(np.median(su_h - st_h)), -22, -8, "orthostatic HR response")
    check("sign dissociation: subjects with MCAv up AND HR down (frac)",
          float(np.mean((su_m.values > st_m.values) & (su_h.values < st_h.values))),
          0.95, 1.0, "the dissociation must hold in essentially every subject")

    print("\nMetabolic and cardiac")
    v_end = stage(cyc, 780, 835, "vo2_true"); v_rec = stage(cyc, 900, 955, "vo2_true")
    check("VO2 fall from end-exercise to late recovery (%)",
          float(np.median(100 * (v_rec / v_end - 1))), -80, -50, "VO2 off-kinetics")
    h_end = stage(cyc, 780, 835, "hr_true"); h_rec = stage(cyc, 900, 955, "hr_true")
    check("HR fall from end-exercise to late recovery (bpm)",
          float(np.median(h_rec - h_end)), -100, -45, "HR recovery kinetics")
    rs = []
    for sid, g in cyc.groupby("sid"):
        g = g[g.t >= 130]
        rs.append(np.corrcoef(g.hr_true, g.vo2_true)[0, 1])
    check("within-session HR vs VO2 correlation (latent)",
          float(np.median(rs)), 0.85, 0.99,
          "HR-VO2 within-session r in the literature is high, not 0.55")

    print("\nSensor model")
    sub = pd.read_csv(os.path.join(DATA, "subjects.csv"))
    check("extracranial fraction beta, cohort median",
          float(sub.beta.median()), 0.35, 0.60, "scalp contamination of temple optics")
    check("beta spread across cohort (max - min)",
          float(sub.beta.max() - sub.beta.min()), 0.25, 0.60,
          "the cohort must span the contamination range, not sit at one value")

    out = pd.DataFrame(CHECKS)
    os.makedirs(os.path.join(DATA, "..", "results"), exist_ok=True)
    out.to_csv(os.path.join(DATA, "..", "results", "verification.csv"), index=False)
    print(f"\n{len(CHECKS) - FAILS}/{len(CHECKS)} checks passed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
