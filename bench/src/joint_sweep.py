"""
Correlated joint parameter sweep.

Replaces the one-at-a-time sweep, which probed 39 points on the axes through one
corner of a 19-dimensional space and, worse, held coupled parameters at baseline
while moving their partners. That constructs people who cannot exist: a subject
with a large orthostatic stroke-volume drop necessarily also has a large
compensatory tachycardia and strong cutaneous vasoconstriction, because all three
are the same reflex.

Here the parameters are grouped into physiological clusters. Each cluster draws
one latent, and every member maps monotonically onto its registered range from
that latent, so a draw is a coherent phenotype rather than an arbitrary corner.
Parameters with no physiological partner are drawn independently across their range.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import tbgen as tb, params as pr, robustness as rb

N_DRAWS = 400
N_SUB = 10
SEED = 4477

# cluster -> {parameter: direction}, +1 means the parameter rises with the latent
CLUSTERS = {
    # magnitude of the orthostatic reflex: bigger venous pooling on standing means
    # a bigger supine stroke-volume recovery, more compensatory tachycardia, and
    # stronger standing cutaneous vasoconstriction
    "orthostatic": {"supine_sv_rise": +1, "orthostatic_hr": +1, "skin_supine_rise": +1},
    # ventilatory responsiveness: an early hyperventilator has a LOW threshold and
    # a STEEP washout, so the two move opposite ways on one latent
    "ventilatory": {"hypervent_threshold": -1, "hypervent_slope": +1},
    # how hard the cycling block actually was, and the cardiovascular gain with it
    "exercise": {"w_max_mean": +1, "map_exercise_gain": +1, "co_exercise_gain": +1},
    # quality of one recording session
    "recording": {"tcd_noise": +1, "tcd_drift_sd": +1, "tcd_dropout": +1},
}
INDEPENDENT = ["co2_react_hypo_ratio", "autoreg_gain", "mcav_supine_extra",
               "mcav_co_coupling", "hydrostatic", "skin_co_coupling",
               "skin_thermal_gain", "co2_rise_slope"]


def draw(rng):
    ov, latents = {}, {}
    for cname, members in CLUSTERS.items():
        u = rng.random()
        latents[cname] = u
        for k, direction in members.items():
            spec = pr.P[k]
            f = u if direction > 0 else (1.0 - u)
            ov[k] = spec["lo"] + f * (spec["hi"] - spec["lo"])
    for k in INDEPENDENT:
        spec = pr.P[k]
        ov[k] = rng.uniform(spec["lo"], spec["hi"])
    return ov, latents


def main():
    rng = np.random.default_rng(SEED)
    rb.N_SUB = N_SUB
    rows = []
    for i in range(N_DRAWS):
        ov, latents = draw(rng)
        wmm = ov.pop("w_max_mean")
        tb.OVERRIDE = ov
        m = rb.measure(wmm)
        tb.OVERRIDE = {}
        rec = dict(latents)
        rec["w_max_mean"] = wmm
        for pn in ("tp_cycle", "tp_posture"):
            for k, v in m[pn].items():
                rec[f"{pn}.{k}"] = v
        rec["verdict"] = {c: bool(f(m)) for c, f in rb.CLAIMS.items()}
        rows.append(rec)
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{N_DRAWS}", flush=True)

    print(f"\n{N_DRAWS} coherent joint draws x {N_SUB} subjects x 2 protocols\n")
    print("CLAIM SURVIVAL ACROSS THE JOINT SPACE")
    for c in rb.CLAIMS:
        held = [r["verdict"][c] for r in rows]
        frac = float(np.mean(held))
        bar = "#" * int(round(frac * 40))
        print(f"  {frac*100:5.1f}%  {bar:<40}  {c}")

    print("\nDISTRIBUTION OF EACH FLOOR (|r| vs the simulated Doppler)")
    for pn, lab in (("tp_cycle", "exercise"), ("tp_posture", "posture ")):
        for k in ("scalp", "co", "hr", "ceiling"):
            v = np.abs([r[f"{pn}.{k}"] for r in rows])
            print(f"  {lab} {k:8s}  p5 {np.percentile(v,5):.3f}   median {np.median(v):.3f}"
                  f"   p95 {np.percentile(v,95):.3f}   min {v.min():.3f}")

    print("\nSIGN STABILITY (posture)")
    co = np.array([r["tp_posture.co"] for r in rows])
    hr = np.array([r["tp_posture.hr"] for r in rows])
    print(f"  cardiac output positive in {100*np.mean(co>0):.1f}% of draws"
          f"  (min {co.min():+.3f})")
    print(f"  heart rate negative in    {100*np.mean(hr<0):.1f}% of draws"
          f"  (max {hr.max():+.3f})")

    print("\nWHAT DRIVES FAILURE (latent correlation with claim holding)")
    for c in rb.CLAIMS:
        held = np.array([r["verdict"][c] for r in rows], float)
        if held.min() == held.max():
            continue
        drivers = []
        for lat in list(CLUSTERS) + ["w_max_mean"]:
            x = np.array([r[lat] for r in rows], float)
            drivers.append((abs(np.corrcoef(x, held)[0, 1]), lat,
                            np.corrcoef(x, held)[0, 1]))
        drivers.sort(reverse=True)
        top = ", ".join(f"{l} r={v:+.2f}" for _, l, v in drivers[:3])
        print(f"  {c[:52]:54s} {top}")

    json.dump([{k: (v if not isinstance(v, dict) else v) for k, v in r.items()}
               for r in rows], open('results/joint_sweep.json', 'w'), default=float)
    return rows


if __name__ == "__main__":
    main()
