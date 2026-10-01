"""
Generate the synthetic cohort.

data/streams_1hz.csv   tidy 1 Hz table: candidates, references, latent truth
data/raw64/*.npz       64 Hz sensor traces (long + short channel) and latent truth
data/subjects.csv      per-subject parameters

Conditions: real | sham_ic (cerebral coupling ablated) | sham_symp (sympathetic
complexity pathway ablated). All three share seeds, so they differ only in the
ablated pathway.
"""
import os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import tbgen as tb
import pipeline as pl

N_SUBJECTS = 24
COHORT_SEED = 20260917
OUT = os.path.join(os.path.dirname(__file__), "..", "data")


def run():
    os.makedirs(os.path.join(OUT, "raw64"), exist_ok=True)
    master = np.random.default_rng(COHORT_SEED)
    subs = [tb.draw_subject(i, master) for i in range(N_SUBJECTS)]
    rows, t0 = [], time.time()

    for sub in subs:
        for pname, pcfg in tb.PROTOCOLS.items():
            for cond in pl.CONDITIONS:
                R = pl.run_cell(sub, pname, pcfg, cond)
                n, ph, ref = R["n"], R["phys"], R["ref"]
                sl = (np.arange(n) * int(tb.FS)).astype(int)

                np.savez_compressed(
                    os.path.join(OUT, "raw64", f"s{sub['sid']:02d}_{pname}_{cond}.npz"),
                    long_ch=R["raw"].astype(np.float32),
                    short_ch=R["short"].astype(np.float32),
                    mcav_true=ph["mcav"].astype(np.float32),
                    vo2_true=ph["vo2"].astype(np.float32),
                    hr_true=ph["hr"].astype(np.float32),
                    symp_true=ph["symp"].astype(np.float32),
                    fs=tb.FS, beta=sub["beta"], sid=sub["sid"],
                    protocol=pname, cond=cond)

                d = dict(sid=sub["sid"], protocol=pname, cond=cond,
                         t=np.arange(n, dtype=float), beta=sub["beta"])
                for c in pl.FLOW_CANDIDATES + pl.CX_CANDIDATES:
                    d[c] = np.asarray(R[c], float)[:n]
                d.update(tcd=ref["tcd"][:n], mcav_true=ref["true_mcav"][:n],
                         vo2_bxb=ref["vo2_bxb"][:n], vo2_true=ref["true_vo2"][:n],
                         ecg_hr=ref["ecg_hr"][:n], hr_true=ref["true_hr"][:n],
                         symp_true=ph["symp"][sl], W=ph["W"][sl],
                         upright=ph["upright"][sl])
                rows.append(pd.DataFrame(d))
        print(f"  subject {sub['sid']:02d}  beta={sub['beta']:.2f}  "
              f"{time.time()-t0:6.1f}s", flush=True)

    df = pd.concat(rows, ignore_index=True)
    # map complexity candidates onto the 1-250 scale Temple reports, per stream
    for c in pl.CX_CANDIDATES:
        df[c + "_250"] = np.nan
        for _, g in df.groupby(["protocol", "cond"]):
            df.loc[g.index, c + "_250"] = tb.to_temple_scale(g[c].values)
    df.to_csv(os.path.join(OUT, "streams_1hz.csv"), index=False)
    pd.DataFrame(subs).to_csv(os.path.join(OUT, "subjects.csv"), index=False)
    print(f"\n{len(df):,} rows | {df.sid.nunique()} subjects | {time.time()-t0:.1f}s")
    return df


if __name__ == "__main__":
    run()
