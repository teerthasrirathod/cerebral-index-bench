"""
Two-sided bracket, computed over the candidate family.

Floor  : the same candidate, same analysis, on a stream whose cerebral coupling
         (or sympathetic-complexity pathway) has been ablated with everything else
         held bit-identical. What the candidate scores with the signal removed.
Ceiling: an oracle handed the latent truth, degraded only by the reference
         channel's own noise, dropout and rate. Nothing measured through that
         reference beats it.

Analysis rules are fixed in notes/SPEC.md §8 and are not tuned here.
"""
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import pipeline as pl

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "data")
RES = os.path.join(HERE, "..", "results")

WARMUP_S = 130          # rolling windows are full by here
TRANS = (295, 315)      # postural transition, excluded from protocol B


def epoch_mask(df, protocol, drop_recovery=False):
    m = df["t"] >= WARMUP_S
    if protocol == "posture":
        m &= ~df["t"].between(*TRANS)
    if protocol == "cycle" and drop_recovery:
        m &= df["t"] < 840
    return m


def _r(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    k = np.isfinite(a) & np.isfinite(b)
    if k.sum() < 60 or np.std(a[k]) < 1e-12 or np.std(b[k]) < 1e-12:
        return np.nan
    return float(np.corrcoef(a[k], b[k])[0, 1])


def _r_lag(a, b, lo=-30, hi=60):
    """Best correlation over an integer-second lag, and the lag that gave it."""
    best, bl = np.nan, 0
    for L in range(lo, hi + 1):
        aa = np.roll(np.asarray(a, float), L)
        aa[:max(L, 0)] = np.nan
        if L < 0:
            aa[L:] = np.nan
        v = _r(aa, b)
        if np.isfinite(v) and (not np.isfinite(best) or abs(v) > abs(best)):
            best, bl = v, L
    return best, bl


def _r2(y, Xcols):
    y = np.asarray(y, float)
    X = np.column_stack([np.asarray(c, float) for c in Xcols])
    k = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if k.sum() < 60:
        return np.nan
    X = np.column_stack([np.ones(k.sum()), X[k]])
    yy = y[k]
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    ss_res = float(np.sum((yy - X @ beta) ** 2))
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    return 1 - ss_res / ss_tot if ss_tot > 0 else np.nan


def within_subject(df, protocol, cond, xcol, ycol, drop_recovery=False, lag=False):
    """Per-subject correlation, then the cohort median. Matches what the preprint reports."""
    d = df[(df.protocol == protocol) & (df.cond == cond)]
    vals, lags = [], []
    for sid, g in d.groupby("sid"):
        g = g[epoch_mask(g, protocol, drop_recovery)]
        if len(g) < 120:
            continue
        if lag:
            v, L = _r_lag(g[xcol], g[ycol]); lags.append(L)
        else:
            v = _r(g[xcol], g[ycol])
        if np.isfinite(v):
            vals.append(v)
    if not vals:
        return dict(median=np.nan, q1=np.nan, q3=np.nan, n=0, med_lag=np.nan)
    a = np.abs(vals)
    return dict(median=float(np.median(a)), q1=float(np.percentile(a, 25)),
                q3=float(np.percentile(a, 75)), n=len(a),
                med_lag=float(np.median(lags)) if lags else np.nan)


def pooled(df, protocol, cond, xcol, ycol, drop_recovery=False):
    d = df[(df.protocol == protocol) & (df.cond == cond)]
    d = d[epoch_mask(d, protocol, drop_recovery)]
    return abs(_r(d[xcol], d[ycol]))


def smooth(df, col, win):
    """Rolling mean within each subject-protocol-condition block."""
    return df.groupby(["sid", "protocol", "cond"])[col].transform(
        lambda s: s.rolling(win, min_periods=max(2, win // 2), center=True).mean())


# ---------------------------------------------------------------------------
def main():
    os.makedirs(RES, exist_ok=True)
    df = pd.read_csv(os.path.join(DATA, "streams_1hz.csv"))
    for w in (10, 30):
        df[f"vo2_{w}s"] = smooth(df, "vo2_bxb", w)

    out = {}

    # ---- A. flow candidates vs TCD -----------------------------------------
    rows = []
    for protocol in ("cycle", "posture"):
        # ceiling: latent MCAv against the TCD trace that measures it
        ceil = within_subject(df, protocol, "real", "mcav_true", "tcd")
        rows.append(dict(protocol=protocol, candidate="ORACLE (ceiling)",
                         cond="real", **ceil))
        for cand in pl.FLOW_CANDIDATES:
            for cond in ("real", "sham_ic"):
                rows.append(dict(protocol=protocol, candidate=cand, cond=cond,
                                 **within_subject(df, protocol, cond, cand, "tcd")))
    flow = pd.DataFrame(rows)
    flow.to_csv(os.path.join(RES, "bracket_flow.csv"), index=False)

    # ---- B. complexity candidates vs VO2 -----------------------------------
    rows = []
    for protocol in ("cycle", "posture"):
        for cond in ("real", "sham_symp", "sham_ic"):
            d = df[(df.protocol == protocol) & (df.cond == cond)]
            hr_only, full, oracle = [], [], []
            for sid, g in d.groupby("sid"):
                g = g[epoch_mask(g, protocol)]
                if len(g) < 120:
                    continue
                hr_only.append(_r2(g["vo2_bxb"], [g["ecg_hr"]]))
                oracle.append(_r2(g["vo2_bxb"], [g["ecg_hr"], g["symp_true"]]))
                full.append({c: _r2(g["vo2_bxb"], [g["ecg_hr"], g[c]])
                             for c in pl.CX_CANDIDATES})
            base = float(np.nanmedian(hr_only))
            rows.append(dict(protocol=protocol, cond=cond, candidate="HR only",
                             r2=base, dr2=0.0,
                             r_vs_vo2=within_subject(df, protocol, cond,
                                                     "ecg_hr", "vo2_bxb")["median"]))
            orc = float(np.nanmedian(oracle))
            rows.append(dict(protocol=protocol, cond=cond,
                             candidate="ORACLE symp (ceiling)", r2=orc,
                             dr2=orc - base,
                             r_vs_vo2=within_subject(df, protocol, cond,
                                                     "symp_true", "vo2_bxb")["median"]))
            for c in pl.CX_CANDIDATES:
                v = float(np.nanmedian([f[c] for f in full]))
                rows.append(dict(protocol=protocol, cond=cond, candidate=c, r2=v,
                                 dr2=v - base,
                                 r_vs_vo2=within_subject(df, protocol, cond,
                                                         c, "vo2_bxb")["median"]))
    cx = pd.DataFrame(rows)
    cx.to_csv(os.path.join(RES, "bracket_cx.csv"), index=False)

    # ---- C. free-parameter sweep -------------------------------------------
    rows = []
    for protocol in ("cycle", "posture"):
        for cond in ("real", "sham_symp"):
            for vcol in ("vo2_bxb", "vo2_10s", "vo2_30s", "vo2_true"):
                for agg in ("within", "pooled"):
                    for drop_rec in (False, True):
                        for lag in (False, True):
                            if agg == "pooled" and lag:
                                continue
                            rec = dict(protocol=protocol, cond=cond, vo2=vcol,
                                       agg=agg, drop_recovery=drop_rec, lag=lag)
                            for c in ("ecg_hr",) + pl.CX_CANDIDATES:
                                if agg == "within":
                                    rec[c] = within_subject(df, protocol, cond, c,
                                                            vcol, drop_rec, lag)["median"]
                                else:
                                    rec[c] = pooled(df, protocol, cond, c, vcol, drop_rec)
                            rows.append(rec)
    sweep = pd.DataFrame(rows)
    sweep.to_csv(os.path.join(RES, "sweep.csv"), index=False)

    # ---- D. headline summary ----------------------------------------------
    out["flow"] = flow.to_dict("records")
    out["cx"] = cx.to_dict("records")
    # the combination that best reproduces Temple's reported 0.93 vs 0.55 pair
    s = sweep.copy()
    for c in pl.CX_CANDIDATES:
        s[f"gap_{c}"] = (s[c] - 0.93).abs() + (s["ecg_hr"] - 0.55).abs()
    best = {}
    for c in pl.CX_CANDIDATES:
        i = s[f"gap_{c}"].idxmin()
        if np.isfinite(s.loc[i, f"gap_{c}"]):
            best[c] = s.loc[i, ["protocol", "cond", "vo2", "agg", "drop_recovery",
                                "lag", "ecg_hr", c]].to_dict()
    out["reproduce_0.93_vs_0.55"] = best
    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)

    print("\n=== A. flow candidates vs TCD, median within-subject |r| ===")
    print(flow.pivot_table(index=["protocol", "candidate"], columns="cond",
                           values="median").round(3).to_string())
    print("\n=== B. complexity candidates vs VO2, delta R-squared over HR ===")
    print(cx.pivot_table(index=["protocol", "candidate"], columns="cond",
                         values="dr2").round(3).to_string())
    print("\n=== B2. raw |r| vs VO2 ===")
    print(cx.pivot_table(index=["protocol", "candidate"], columns="cond",
                         values="r_vs_vo2").round(3).to_string())
    print("\n=== C. analysis choices that reproduce 0.93 vs 0.55 ===")
    for k, v in best.items():
        print(f"  {k}: {v}")
    return flow, cx, sweep


if __name__ == "__main__":
    main()
