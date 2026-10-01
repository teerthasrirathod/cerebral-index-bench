"""
Bench vs Temple TD-NIRS preprint (bioRxiv 10.64898/2026.09.17.751664).
Rules: tdnirs/PREREG-01OCT-tdnirs.md. 1 Hz reduced model; nothing here is
Temple's formula. Two architectures (FLOW, VOLUME), their exact pipeline.
"""
import json, sys
import numpy as np
from scipy import signal as sps

SEED = 20261001
NSESS = 20
PROTO = {  # PAPER
    "hdt":    dict(block=180, ramp=10, win=90),
    "squat":  dict(block=90,  ramp=4,  win=45),
    "supine": dict(block=180, ramp=10, win=90),
}
NIRS_NOISE = 0.15   # ASSUM, x block amplitude
TAU_F = 6.0         # FLOW response, autoregulation speed
ALPHA = 0.35        # PAPER, forward-backward EMA


def task(p):
    B, R = PROTO[p]["block"], PROTO[p]["ramp"]
    T = 5 * B
    t = np.arange(T, dtype=float)
    u = np.zeros(T)
    for k in (1, 3):                      # task blocks 2 and 4
        on, off = k * B, (k + 1) * B
        u += np.clip((t - on) / R, 0, 1) - np.clip((t - off) / R, 0, 1)
    return t, np.clip(u, 0, 1)


def fo(x, tau):
    a = 1.0 / tau
    return sps.lfilter([a], [1, -(1 - a)], x, zi=[(1 - a) * x[0]])[0]


def band_noise(rng, n, tau):
    w = rng.normal(0, 1, n + 200)
    y = fo(w, tau)[200:]
    return y / y.std()


def causal_ma(x, N):
    if N <= 1:
        return x.copy()
    c = np.cumsum(np.concatenate([[0.0], x]))
    out = np.empty_like(x)
    for i in range(len(x)):
        j0 = max(0, i + 1 - N)
        out[i] = (c[i + 1] - c[j0]) / (i + 1 - j0)
    return out


def ema_fb(x, a=ALPHA):
    y = sps.lfilter([a], [1, -(1 - a)], x, zi=[(1 - a) * x[0]])[0]
    y = sps.lfilter([a], [1, -(1 - a)], y[::-1], zi=[(1 - a) * y[-1]])[0]
    return y[::-1]


def pearson(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def lagged(ref, dev, lo=-30, hi=30):
    """positive lag = NIRS leads the index: corr(ref[t], dev[t+L])"""
    best, bl = -2, 0
    n = len(ref)
    for L in range(lo, hi + 1):
        if L >= 0:
            r = pearson(ref[:n - L], dev[L:])
        else:
            r = pearson(ref[-L:], dev[:n + L])
        if r > best:
            best, bl = r, L
    return best, bl


def z(x):
    return (x - x.mean()) / x.std()


def transitions(p, ref, dev):
    B, W = PROTO[p]["block"], PROTO[p]["win"]
    zr, zd = z(ref), z(dev)
    rr, dd = [], []
    for k in (1, 3):
        t0 = k * B
        rr.append(zr[t0:t0 + W].mean() - zr[t0 - W:t0].mean())
        dd.append(zd[t0:t0 + W].mean() - zd[t0 - W:t0].mean())
    return float(np.mean(rr)), float(np.mean(dd))


def session(p, arch, f, N, sd, rng, tau_box=None):
    t, u = task(p)
    n = len(t)
    tauB = rng.uniform(15, 29)                         # CAL
    s_flow = band_noise(rng, n, 4.0)
    s_vol = fo(s_flow, tauB)
    blk = fo(u, tauB)
    # scale spontaneous so block share of brain variance = f
    k = np.sqrt(blk.var() * (1 - f) / f) / s_vol.std()
    s_vol, s_flow = k * s_vol, k * s_flow
    brain = blk + s_vol
    ref = brain + rng.normal(0, NIRS_NOISE, n)
    if arch == "VOLUME":
        lat = brain
    elif arch == "FLOW":
        lat = fo(u, TAU_F) + s_flow
    elif arch == "TIMING":                              # floor: protocol only
        lat = fo(u, tau_box)
    elif arch == "ORACLE":
        lat = brain
    else:
        raise ValueError(arch)
    lat = lat / blk.std()                              # common scale
    dev = causal_ma(lat, N) + rng.normal(0, sd, n)
    R, D = ema_fb(ref), ema_fb(dev)
    rl, L = lagged(R, D)
    r0 = pearson(R, D)
    tr, td = transitions(p, R, D)
    return rl, L, r0, tr, td


def cohort(p, arch, f, N, sd, seed, reps=8, tau_box=None, lagsearch=True):
    out = []
    for rep in range(reps):
        rng = np.random.default_rng(seed + 1000 * rep)
        S = [session(p, arch, f, N, sd, rng, tau_box) for _ in range(NSESS)]
        S = np.array(S)
        conc = np.mean(np.sign(S[:, 3]) == np.sign(S[:, 4]))
        out.append(dict(r_lag=np.median(S[:, 0]), lag=np.median(S[:, 1]),
                        r0=np.median(S[:, 2]), ba=np.mean(S[:, 4] - S[:, 3]),
                        conc=conc))
    return {k: float(np.median([o[k] for o in out])) for k in out[0]}


def calibrate(arch, f, seed):
    """N to hit HDT lag +14, sd to hit HDT r_lag 0.912. HDT only."""
    best = None
    for N in range(0, 91, 2):
        c = cohort("hdt", arch, f, N, 0.3, seed, reps=3)
        d = abs(c["lag"] - 14)
        if best is None or d < best[0]:
            best = (d, N)
        if c["lag"] > 18:
            break
    N = best[1]
    lo, hi = 0.0, 3.0
    c0 = cohort("hdt", arch, f, N, 0.0, seed, reps=4)
    if c0["r_lag"] < 0.912:
        return N, 0.0, c0, False                       # cannot reach 0.912 at all
    for _ in range(14):
        mid = (lo + hi) / 2
        c = cohort("hdt", arch, f, N, mid, seed, reps=4)
        if c["r_lag"] > 0.912:
            lo = mid
        else:
            hi = mid
    sd = (lo + hi) / 2
    return N, sd, cohort("hdt", arch, f, N, sd, seed, reps=8), True


if __name__ == "__main__":
    res = {}
    for f in (0.60, 0.75, 0.90):
        for arch in ("VOLUME", "FLOW"):
            N, sd, chdt, ok = calibrate(arch, f, SEED)
            row = dict(N=N, sd=round(sd, 3), reached_0912=ok, hdt=chdt)
            for p in ("squat", "supine"):
                row[p] = cohort(p, arch, f, N, sd, SEED + 7, reps=8)
            res[f"{arch}_f{f}"] = row
            print(arch, f, json.dumps(row), flush=True)
        # floor and ceiling with VOLUME's fitted noise
        sdv = res[f"VOLUME_f{f}"]["sd"]
        for p in ("hdt", "squat", "supine"):
            best = None
            for tb in (5, 10, 15, 20, 25, 30, 40):
                c = cohort(p, "TIMING", f, 0, sdv, SEED + 11, reps=6, tau_box=tb)
                if best is None or c["r_lag"] > best[1]["r_lag"]:
                    best = (tb, c)
            orc = cohort(p, "ORACLE", f, 0, sdv, SEED + 13, reps=6)
            res[f"FLOOR_{p}_f{f}"] = dict(tau=best[0], **best[1])
            res[f"CEIL_{p}_f{f}"] = orc
            print("floor/ceil", p, f, best[0], json.dumps(best[1]), json.dumps(orc), flush=True)
    json.dump(res, open("results_tdnirs.json", "w"), indent=1)
