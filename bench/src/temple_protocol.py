"""
The bench re-run against the published protocol.

Kumar, van Rosmalen, Gupta, Sharma, Gupta, Panda and Gupta,
"Concordance Between a Temple-Worn Optical Wearable and Transcranial Doppler
During Exercise and Postural Transitions in Healthy Adults",
medRxiv 10.64898/2026.09.02.26362022, posted 4 September 2026.

Their block structure, their sample size, their analysis rules. The floor is the
same signal with the intracranial term ablated and everything else held
bit-identical, run through the identical pipeline.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import tbgen as tb

N = 20                  # analysed per protocol in the paper
COHORT_SEED = 20260918


def z_full_session(x):
    """Their normalisation: z-score on the mean and SD of the complete trace."""
    x = np.asarray(x, float)
    sd = np.nanstd(x)
    return (x - np.nanmean(x)) / (sd if sd > 0 else 1.0)


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    k = np.isfinite(a) & np.isfinite(b)
    if k.sum() < 30 or np.std(a[k]) == 0 or np.std(b[k]) == 0:
        return np.nan
    return float(np.corrcoef(a[k], b[k])[0, 1])


def lag_optimised(a, b, lo=-30, hi=30):
    """Their lag sweep: -30 to +30 s in 1 s steps, keep the maximum coefficient."""
    best, bl = -np.inf, 0
    for L in range(lo, hi + 1):
        aa = np.roll(np.asarray(a, float), L).astype(float)
        if L > 0:
            aa[:L] = np.nan
        elif L < 0:
            aa[L:] = np.nan
        v = pearson(aa, b)
        if np.isfinite(v) and v > best:
            best, bl = v, L
    return best, bl


def run_subject(sub, pname, cond, w_max):
    cfg = tb.TEMPLE_PROTOCOLS[pname]
    phys = tb.physiology(sub, cfg["kind"], cfg["dur"],
                         np.random.default_rng(sub["seed"] + 101),
                         w_max=w_max)
    long_ch, short_ch = tb.sensor(sub, phys, cfg["kind"],
                                  np.random.default_rng(sub["seed"] + 4441),
                                  ablate_ic=(cond == "sham_ic"))
    ref = tb.references(sub, phys, np.random.default_rng(sub["seed"] + 7717))
    n = len(ref["tcd"])
    return dict(
        n=n, phys=phys,
        tcd=ref["tcd"], hr=ref["ecg_hr"],
        true_mcav=ref["true_mcav"],
        fc_pi=tb.flow_candidate(long_ch, mode="pi")[:n],
        fc_shape=tb.flow_candidate(long_ch, mode="shape")[:n],
        fc_ssr=tb.flow_candidate(long_ch, short_ch, mode="ssr")[:n],
    )


def main():
    master = np.random.default_rng(COHORT_SEED)
    subs = [tb.draw_subject(i, master) for i in range(N)]
    # "fast seated cycling at the participant's maximum self-selected cadence":
    # one bout, intensity self-selected, so it is drawn per subject and spans
    # the hyperventilation threshold rather than sitting on one side of it.
    w_max = {s["sid"]: float(np.clip(master.normal(0.50, 0.08), 0.34, 0.68)) for s in subs}

    CANDS = ("fc_pi", "fc_shape", "fc_ssr")
    out = {}
    for pname in tb.TEMPLE_PROTOCOLS:
        res = {c: {"real": [], "sham_ic": []} for c in CANDS}
        res["ORACLE"] = {"real": []}
        res["HR"] = {"real": []}
        lagres = {c: {"real": [], "sham_ic": []} for c in CANDS}
        for sub in subs:
            cells = {c: run_subject(sub, pname, c, w_max[sub["sid"]])
                     for c in ("real", "sham_ic")}
            tcd_z = z_full_session(cells["real"]["tcd"])
            res["ORACLE"]["real"].append(
                pearson(z_full_session(cells["real"]["true_mcav"]), tcd_z))
            res["HR"]["real"].append(pearson(z_full_session(cells["real"]["hr"]), tcd_z))
            for c in CANDS:
                for cond in ("real", "sham_ic"):
                    sig = z_full_session(cells[cond][c])
                    res[c][cond].append(pearson(sig, tcd_z))
                    lagres[c][cond].append(lag_optimised(sig, tcd_z)[0])
        out[pname] = dict(
            zero={k: {c: float(np.nanmedian(v)) for c, v in d.items()}
                  for k, d in res.items()},
            lag={k: {c: float(np.nanmedian(v)) for c, v in d.items()}
                 for k, d in lagres.items()},
            raw=res)
    return out, subs, w_max


if __name__ == "__main__":
    out, subs, w_max = main()
    REPORTED = {"tp_cycle": (0.795, 0.837), "tp_posture": (0.799, 0.817)}
    LBL = {"fc_pi": "perfusion index (AC/DC)", "fc_shape": "beat morphology",
           "fc_ssr": "short-separation regressed", "ORACLE": "oracle ceiling",
           "HR": "heart rate alone"}
    for pname, d in out.items():
        rep_z, rep_l = REPORTED[pname]
        title = "EXERCISE  240 rest / 240 cycle / 600 recovery" if pname == "tp_cycle" \
            else "POSTURE   stand / supine / stand / supine, 240 s each"
        print(f"\n=== {title} ===")
        print(f"{'':32s} {'zero-lag':>10s} {'lag-opt':>10s}")
        for k in ("ORACLE", "HR", "fc_pi", "fc_shape", "fc_ssr"):
            z = d["zero"][k]
            for cond in (["real"] if k in ("ORACLE", "HR") else ["real", "sham_ic"]):
                lg = d["lag"].get(k, {}).get(cond, float("nan"))
                tag = LBL[k] + ("" if k in ("ORACLE", "HR")
                                else ("  intact" if cond == "real" else "  ABLATED"))
                print(f"  {tag:30s} {z[cond]:10.3f} {lg:10.3f}")
        print(f"  {'-> PAPER REPORTS':30s} {rep_z:10.3f} {rep_l:10.3f}")
