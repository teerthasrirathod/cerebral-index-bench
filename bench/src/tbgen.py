"""
Synthetic temple-wearable cohort generator.

Latent physiology -> 64 Hz temple sensor waveform -> 1 Hz derived streams,
plus simulated TCD, metabolic-cart and ECG reference channels.

Every component is anchored to published physiology, not to Temple's data.
See notes/SPEC.md. Nothing here is a model of their device; it is a model of
the measurement situation their device is in.
"""
import math
import numpy as np
from scipy import signal as sps

OVERRIDE = {}      # robustness harness writes here; empty means published defaults


def _p(key, default):
    """Load-bearing constant. Registered in params.py with value, unit and source."""
    return OVERRIDE.get(key, default)


FS = 64.0          # raw sensor rate, Hz
FS_OUT = 1.0       # derived stream rate, Hz

PROTOCOLS = {
    "cycle": dict(dur=960.0, kind="cycle"),
    "posture": dict(dur=600.0, kind="posture"),
}

# Kumar / van Rosmalen et al., medRxiv 10.64898/2026.09.02.26362022, as described
# in their Section 2.3. Reproduced exactly: block durations, block order, and the
# single bout of self-selected-maximum cycling rather than a graded ramp.
TEMPLE_PROTOCOLS = {
    "tp_cycle": dict(dur=1080.0, kind="tp_cycle"),     # 240 rest / 240 cycle / 600 recovery
    "tp_posture": dict(dur=960.0, kind="tp_posture"),  # 4 x 240: stand, supine, stand, supine
}


# ----------------------------------------------------------------------------
# subject parameters
# ----------------------------------------------------------------------------
def draw_subject(sid, rng):
    age = rng.uniform(22, 45)
    return dict(
        sid=sid,
        age=age,
        mass=np.clip(rng.normal(68, 11), 45, 100),
        hr_rest=np.clip(rng.normal(58, 7), 44, 76),
        hr_max=208 - 0.7 * age + rng.normal(0, 5),
        vo2_rest=3.5,
        vo2_max=float(np.clip(rng.normal(38, 7), 24, 58)),
        tau_vo2=float(np.clip(rng.normal(35, 6), 22, 50)),
        tau_hr=float(np.clip(rng.normal(28, 5), 18, 42)),
        mcav_rest=float(np.clip(rng.normal(58, 9), 40, 80)),
        co2_react=float(np.clip(rng.normal(0.034, 0.006), 0.020, 0.050)),  # frac / mmHg
        # fraction of the temple-surface signal that is extracranial (scalp + STA)
        beta=float(0.25 + 0.50 * rng.beta(2.2, 2.6)),
        head_h=float(np.clip(rng.normal(28, 3), 20, 36)),   # cm above heart, standing
        cadence=float(rng.normal(1.333, 0.05)),             # Hz, 80 rpm
        probe_gain=float(np.exp(rng.normal(0, 0.12))),      # TCD insonation angle error
        snr_db=float(np.clip(rng.normal(21, 3), 14, 28)),
        rr_rest=float(np.clip(rng.normal(14, 2), 10, 19)),  # breaths/min
        seed=int(rng.integers(0, 2**31 - 1)),
    )


# ----------------------------------------------------------------------------
# workload profiles
# ----------------------------------------------------------------------------
def workload(t, kind, w_max=0.70):
    """W(t) in [0,1] of VO2max, and posture: 1 = upright, 0 = supine."""
    W = np.zeros_like(t)
    upright = np.ones_like(t)
    if kind == "tp_cycle":
        # 240 s seated baseline, 240 s cycling at self-selected maximum, 600 s recovery
        W[(t >= 240) & (t < 480)] = w_max
        W[t >= 480] = 0.05
        return W, upright
    if kind == "tp_posture":
        # standing / supine / standing / supine, 240 s each
        # sigmoid ramps between the four blocks
        up = (1.0 / (1.0 + np.exp((t - 240.0) / 1.8))
              + 1.0 / (1.0 + np.exp(-(t - 480.0) / 1.8))
              - 1.0 / (1.0 + np.exp(-(t - 720.0) / 1.8)))
        return W, np.clip(up, 0.0, 1.0)
    if kind == "cycle":
        W[(t >= 180) & (t < 360)] = 0.30
        W[(t >= 360) & (t < 600)] = 0.50
        W[(t >= 600) & (t < 840)] = 0.70
        W[t >= 840] = 0.05
    else:
        # standing -> supine, 10 s transition centred on 305 s
        upright = 1.0 / (1.0 + np.exp((t - 305.0) / 1.8))
        W[:] = 0.0
    return W, upright


def _lowpass_ode(target, tau, dt, delay_s=0.0, y0=None):
    """
    First-order response with optional pure delay.

    Vectorised as a one-pole IIR. Bit-identical to the explicit loop it replaces;
    that loop ran at 64 Hz over tens of thousands of samples and dominated
    runtime, which made joint parameter sweeps impractical.
    """
    if delay_s > 0:
        n = int(round(delay_s / dt))
        target = np.concatenate([np.full(n, target[0]), target[:-n]])
    y_init = target[0] if y0 is None else y0
    a = dt / tau
    if len(target) < 2:
        return np.asarray([y_init], dtype=float)
    tail = sps.lfilter([a], [1.0, -(1.0 - a)], target[1:], zi=[(1.0 - a) * y_init])[0]
    return np.concatenate([[y_init], tail])


# ----------------------------------------------------------------------------
# latent physiology
# ----------------------------------------------------------------------------
def physiology(sub, kind, dur, rng, ablate_symp=False, symp_fixed=0.30, w_max=0.70):
    dt = 1.0 / FS
    t = np.arange(0, dur, dt)
    W, upright = workload(t, kind, w_max)

    # --- VO2, ml/kg/min. Phase II kinetics with cardiodynamic delay.
    vo2_t = sub["vo2_rest"] + (sub["vo2_max"] - sub["vo2_rest"]) * W
    vo2 = _lowpass_ode(vo2_t, sub["tau_vo2"], dt, delay_s=15.0, y0=sub["vo2_rest"])
    vo2 = vo2 * (1.0 - 0.08 * (1.0 - upright))          # supine metabolic rate lower
    vo2 = np.maximum(vo2, 2.6)

    # --- sympathetic drive. Faster onset than HR, much slower withdrawal.
    #     This is the only route by which a complexity metric can carry information
    #     that heart rate does not: during recovery HR falls quickly while
    #     sympathetic drive, and VO2, decay slowly.
    symp = np.empty(len(t)); symp[0] = 0.0
    for i in range(1, len(t)):
        tau = 12.0 if W[i] >= symp[i - 1] else 90.0
        symp[i] = symp[i - 1] + dt / tau * (W[i] - symp[i - 1])
    if kind.endswith("posture"):
        symp = symp + 0.22 * upright          # orthostatic sympathoexcitation
    symp = np.clip(symp, 0, 1.2)

    # --- heart rate
    hr_t = sub["hr_rest"] + (sub["hr_max"] - sub["hr_rest"]) * W
    hr = _lowpass_ode(hr_t, sub["tau_hr"], dt, delay_s=2.0, y0=sub["hr_rest"])
    drift = np.cumsum(np.where(W > 0.40, 0.04 * dt, 0.0))          # cardiac drift
    hr = hr + drift
    if kind.endswith("posture"):
        hr = hr + _p("orthostatic_hr", 14.0) * upright                                   # orthostatic offset

    # beat-interval variability: a regular respiratory component and a
    # complexity-bearing 1/f component, both suppressed by sympathetic drive.
    rr = sub["rr_rest"] / 60.0 * (1.0 + 1.6 * W)
    phase_resp = 2 * np.pi * np.cumsum(rr) * dt
    sup = symp if ablate_symp is False else np.full(len(t), float(symp_fixed))
    rsa_amp = 3.6 * np.clip(1.0 - 0.78 * sup, 0.08, 1.0)
    cplx_amp = 2.4 * np.clip(1.0 - 0.80 * sup, 0.10, 1.0)
    hr_trend = hr.copy()
    pink = _pink_noise(rng, len(t))
    hr = hr + rsa_amp * np.sin(phase_resp) + cplx_amp * pink
    hr = hr + _smooth_noise(rng, len(t), sd=0.5, tau_s=2.0, dt=dt)
    hr = np.clip(hr, 38, 205)

    # --- arterial CO2. Rises with metabolic load, falls once hyperventilation
    #     kicks in above ~60% VO2max. This is what makes MCAv non-monotone.
    load = W
    dco2_t = (_p("co2_rise_slope", 7.0) * load
              - _p("hypervent_slope", 55.0)
              * np.maximum(load - _p("hypervent_threshold", 0.52), 0.0))
    dpaco2 = _lowpass_ode(dco2_t, 22.0, dt, delay_s=8.0, y0=0.0)
    dpaco2 = dpaco2 + 1.4 * (1.0 - upright)          # supine slight hypercapnia

    # --- mean arterial pressure and cerebral perfusion pressure
    map_rest = 88.0
    mapv = map_rest + _p("map_exercise_gain", 25.0) * _lowpass_ode(W, 20.0, dt, y0=0.0)
    # hydrostatic: the head-above-heart column disappears when supine
    cpp = mapv - _p("hydrostatic", 0.77) * sub["head_h"] * upright
    cpp_rest = map_rest - _p("hydrostatic", 0.77) * sub["head_h"]

    # --- cardiac output proxy (normalised)
    co_norm = (hr_trend / sub["hr_rest"]) \
        * (1.0 + _p("co_exercise_gain", 0.35) * _lowpass_ode(W, 25.0, dt, y0=0.0)) \
        * (1.0 + _p("supine_sv_rise", 0.52) * (1.0 - upright))
    co_norm = co_norm / co_norm[: int(30 * FS)].mean()

    # --- MCAv ground truth
    k = np.where(dpaco2 >= 0, sub["co2_react"],
                 sub["co2_react"] * _p("co2_react_hypo_ratio", 0.55))
    auto = _lowpass_ode((cpp - cpp_rest), 6.0, dt, y0=0.0)   # autoregulatory lag
    # CO2 reactivity dominates; perfusion pressure is strongly autoregulated;
    # cardiac output couples weakly. The supine term carries the combined
    # hydrostatic / venous-return effect that autoregulation does not remove.
    mcav = sub["mcav_rest"] * (1 + k * dpaco2) \
        * (1 + _p("autoreg_gain", 0.0015) * auto) \
        * (1 + _p("mcav_co_coupling", 0.035) * (co_norm - 1.0)) \
        * (1 + _p("mcav_supine_extra", 0.060) * (1.0 - upright))
    mcav = mcav + _smooth_noise(rng, len(t), sd=0.9, tau_s=3.0, dt=dt)
    mcav = np.clip(mcav, 18, 140)

    # --- core / skin temperature -> skin blood flow (extracranial driver)
    tc = _lowpass_ode(W, 300.0, dt, y0=0.0)
    skin_flow = 1.0 + _p("skin_thermal_gain", 2.4) * np.maximum(tc - 0.05, 0.0)
    # (a) hydrostatic + loss of orthostatic cutaneous vasoconstriction when supine
    skin_flow = skin_flow * (1.0 + _p("skin_supine_rise", 0.10) * (1.0 - upright))
    # (b) unautoregulated: tracks perfusion pressure and cardiac output directly
    skin_flow = skin_flow * (1.0 + _p("skin_co_coupling", 0.45) * (co_norm - 1.0)) \
        * (1.0 + 0.004 * (mapv - map_rest))
    skin_flow = skin_flow + _smooth_noise(rng, len(t), sd=0.04, tau_s=12.0, dt=dt)

    return dict(t=t, W=W, upright=upright, vo2=vo2, hr=hr, dpaco2=dpaco2,
                map=mapv, cpp=cpp, mcav=mcav, skin=skin_flow, co=co_norm, rr=rr,
                symp=symp, hr_trend=hr_trend)


def _pink_noise(rng, n):
    """Unit-variance 1/f noise. Carries the complexity a regular oscillator has not."""
    w = rng.normal(0, 1, n)
    f = np.fft.rfftfreq(n, d=1.0)
    f[0] = f[1] if len(f) > 1 else 1.0
    y = np.fft.irfft(np.fft.rfft(w) / np.sqrt(f), n=n)
    sd = y.std()
    return y / sd if sd > 0 else y


def _smooth_noise(rng, n, sd, tau_s, dt):
    """Band-limited noise: white noise through a one-pole lag, rescaled to sd."""
    w = rng.normal(0, 1, n)
    a = dt / tau_s
    if n < 2:
        return w
    tail = sps.lfilter([a], [1.0, -(1.0 - a)], w[1:], zi=[(1.0 - a) * w[0]])[0]
    y = np.concatenate([[w[0]], tail])
    s_ = y.std()
    return y * (sd / s_) if s_ > 0 else y


# ----------------------------------------------------------------------------
# forward sensor model
# ----------------------------------------------------------------------------
def _pulse_ic(ph, notch):
    """Intracranial-coupled pulse: sharp early systolic upstroke, distinct notch."""
    return (np.exp(-((ph - 0.22) / 0.085) ** 2)
            + notch * np.exp(-((ph - 0.48) / 0.125) ** 2)
            + 0.14 * np.exp(-((ph - 0.72) / 0.20) ** 2))


def _pulse_ec(ph, notch):
    """Extracranial (scalp + superficial temporal): broader, later, damped notch."""
    return (0.92 * np.exp(-((ph - 0.31) / 0.135) ** 2)
            + 0.62 * notch * np.exp(-((ph - 0.62) / 0.185) ** 2)
            + 0.20 * np.exp(-((ph - 0.86) / 0.24) ** 2))


def sensor(sub, phys, kind, rng, ablate_ic=False):
    """
    64 Hz temple sensor trace.

    ablate_ic=True is the floor condition: the intracranial term is replaced by an
    extracranial-shaped term of the same amplitude budget, with cardiac phase, scalp
    waveform, motion, drift and noise held bit-identical. Same hardware, no brain.
    """
    dt = 1.0 / FS
    t, hr = phys["t"], phys["hr"]
    beta = sub["beta"]

    ph = np.cumsum(hr / 60.0) * dt
    ph = np.mod(ph, 1.0)

    notch = 0.30 + 0.55 * np.clip((phys["map"] - 88.0) / 35.0, 0, 1.4)

    # The superficial temporal pulse arrives later than the intracranial one, and
    # the offset moves with peripheral vasomotor tone. Without this the two
    # components sit at identical cardiac phase and no filter can separate them,
    # which would make short-separation regression impossible by construction.
    skin_n = phys["skin"] / phys["skin"][: int(30 * FS)].mean()
    ptt = _p("ptt_offset", 0.026) + _p("ptt_vasomotor", 0.022) * np.clip(1.4 - skin_n, 0, 1.0)
    ph_ec = np.mod(ph - ptt * phys["hr"] / 60.0, 1.0)

    a_ic = (phys["mcav"] / sub["mcav_rest"]) * (1 + 0.30 * (phys["map"] - 88.0) / 88.0)
    a_ec = phys["skin"] * (1 + 0.18 * (phys["hr_trend"] / sub["hr_rest"] - 1))

    if ablate_ic:
        # same amplitude budget, extracranial content only
        ic_term = (1 - beta) * a_ec * _pulse_ec(ph_ec, notch) * (a_ic.mean() / a_ec.mean())
    else:
        ic_term = (1 - beta) * a_ic * _pulse_ic(ph, notch)
    ec_term = beta * a_ec * _pulse_ec(ph_ec, notch)
    pulsatile = ic_term + ec_term

    # --- motion
    if kind == "cycle":
        f = sub["cadence"]
        m = (0.55 * np.sin(2 * np.pi * f * t)
             + 0.22 * np.sin(2 * np.pi * 2 * f * t + 0.7)
             + 0.11 * np.sin(2 * np.pi * 3 * f * t + 1.9))
        m = m * (0.30 + 1.5 * _lowpass_ode(phys["W"], 4.0, dt, y0=0.0))
    else:
        m = 2.6 * np.exp(-0.5 * ((t - 305.0) / 2.2) ** 2) * np.sin(2 * np.pi * 1.1 * t)
    m = m * 0.35 * np.abs(pulsatile).mean()

    # --- drift: contact random walk + thermal ramp
    rw = np.cumsum(rng.normal(0, 1, len(t))) * dt
    rw = rw / (np.abs(rw).max() + 1e-9)
    thermal = _lowpass_ode(phys["W"], 240.0, dt, y0=0.0)
    d = (0.55 * rw + 0.35 * thermal) * np.abs(pulsatile).mean() * 2.0

    # --- sensor noise
    sig_p = pulsatile.std()
    eps = rng.normal(0, sig_p / (10 ** (sub["snr_db"] / 20.0)), len(t))

    # optical-style DC pedestal: the signal is a light level, strictly positive,
    # with the pulsatile component riding on it as a small AC fraction.
    dc0 = 10.0 * np.abs(pulsatile).mean()
    long_ch = np.maximum(dc0 + pulsatile + m + d + eps, 0.05 * dc0)

    # Short-separation channel. Sees scalp and superficial temporal artery only,
    # shares the motion and most of the drift, has independent sensor noise.
    # Any device that regresses this out of the long channel is doing what the
    # fNIRS literature has required since short-separation regression existed.
    ec_only = 1.30 * a_ec * _pulse_ec(ph_ec, notch)
    eps_s = rng.normal(0, ec_only.std() / (10 ** (sub["snr_db"] / 20.0)), len(t))
    dc0s = 10.0 * np.abs(ec_only).mean()
    short_ch = np.maximum(dc0s + ec_only + 0.92 * m + 0.85 * d + eps_s, 0.05 * dc0s)

    return long_ch, short_ch


# ----------------------------------------------------------------------------
# device-side metrics
# ----------------------------------------------------------------------------
def _bandpass(x, lo, hi, fs=FS, order=3):
    b, a = sps.butter(order, [lo / (fs / 2), hi / (fs / 2)], btype="band")
    return sps.filtfilt(b, a, x)


def flow_candidate(raw, short=None, mode="pi"):
    """
    A candidate cerebral-flow index at 1 Hz. NOT a reconstruction of anyone's
    proprietary metric - one member of a spanning family of plausible constructions,
    used so the bracket can be taken over the family rather than over one guess.

    mode='pi'    : perfusion index. Beat pulse amplitude of the band-passed signal
                   divided by the DC light level in the same window - the standard
                   AC/DC construction every optical wearable uses. Carries no
                   separation of intracranial from extracranial content.
    mode='shape' : early-systolic energy relative to late-diastolic energy. The
                   intracranial and extracranial pulse shapes differ in where their
                   energy sits, so a morphology-based estimator can partially
                   separate them.
    All normalised to the subject's first 60 s.
    """
    ac = _bandpass(raw, 0.5, 8.0)
    b, a = sps.butter(3, 0.20 / (FS / 2), btype="low")
    dc = np.maximum(sps.filtfilt(b, a, raw), 1e-6)

    n_sec = int(len(raw) // FS)
    out = np.full(n_sec, np.nan)
    w = int(5 * FS)                              # 5 s analysis window, 1 s hop

    if mode == "pi":
        for s_i in range(n_sec):
            i0 = max(0, int(s_i * FS) - w // 2)
            seg, segdc = ac[i0:i0 + w], dc[i0:i0 + w]
            if len(seg) < w // 2:
                continue
            pk, _ = sps.find_peaks(seg, distance=int(0.40 * FS))
            tr, _ = sps.find_peaks(-seg, distance=int(0.40 * FS))
            amp = (seg[pk].mean() - seg[tr].mean()) if (len(pk) >= 2 and len(tr) >= 2) \
                else (seg.max() - seg.min())
            out[s_i] = amp / segdc.mean()
    else:
        xe = _bandpass(raw, 2.5, 8.0)            # early-systolic band
        xl = _bandpass(raw, 0.5, 2.2)            # late / diastolic band
        for s_i in range(n_sec):
            i0 = max(0, int(s_i * FS) - w // 2)
            se, sl = xe[i0:i0 + w], xl[i0:i0 + w]
            if len(se) < w // 2:
                continue
            out[s_i] = np.sqrt(np.mean(se ** 2)) / (np.sqrt(np.mean(sl ** 2)) + 1e-12)

    if mode == "shape":
        out = _beat_shape_index(raw, n_sec, w)
    elif mode == "ssr":
        out = _ssr_index(raw, short, n_sec, w)

    out = _fill_nan(out)
    base = np.nanmean(out[:60])
    return out / base if base and np.isfinite(base) else out


def _ssr_index(raw, short, n_sec, w, taps=21):
    """
    Short-separation-regressed perfusion index.

    A multi-tap FIR is fitted from the short channel to the long channel over the
    whole record and subtracted. Fitting once over the record, rather than per
    window, is deliberate: it removes the stationary scalp contribution and leaves
    the part that changes differently from the scalp. What remains is not clean -
    scalp leakage survives wherever the scalp changed - and that leakage is the
    honest residual a real device has to live with.

    A device that skips this step is reporting scalp perfusion with a cerebral
    label on it.
    """
    if short is None:
        return np.full(n_sec, np.nan)
    L = _bandpass(raw, 0.5, 8.0)
    S = _bandpass(short, 0.5, 8.0)
    b, a = sps.butter(3, 0.20 / (FS / 2), btype="low")
    dc = np.maximum(sps.filtfilt(b, a, raw), 1e-6)

    X = np.column_stack([np.roll(S, k) for k in range(-(taps // 2), taps // 2 + 1)])
    X[: taps, :] = 0.0
    X[-taps:, :] = 0.0
    coef, *_ = np.linalg.lstsq(X, L, rcond=None)
    resid = L - X @ coef

    out = np.full(n_sec, np.nan)
    for s_i in range(n_sec):
        i0 = max(0, int(s_i * FS) - w // 2)
        seg, segdc = resid[i0:i0 + w], dc[i0:i0 + w]
        if len(seg) < w // 2:
            continue
        pk, _ = sps.find_peaks(seg, distance=int(0.40 * FS))
        tr, _ = sps.find_peaks(-seg, distance=int(0.40 * FS))
        amp = (seg[pk].mean() - seg[tr].mean()) if (len(pk) >= 2 and len(tr) >= 2) \
            else (seg.max() - seg.min())
        out[s_i] = amp / segdc.mean()
    return out


def _beat_shape_index(raw, n_sec, w):
    """
    Heart-rate-invariant morphology index.

    Beats are segmented, each resampled onto a fixed 0-1 phase grid, and averaged
    across the window. The index is early-cycle area over mid-cycle area. Because
    every beat is phase-normalised first, the index does not move simply because
    heart rate moved - which the band-ratio version does.
    """
    x = _bandpass(raw, 0.5, 10.0)
    pk, _ = sps.find_peaks(x, distance=int(0.36 * FS), prominence=0.25 * np.std(x))
    out = np.full(n_sec, np.nan)
    if len(pk) < 6:
        return out
    grid = np.linspace(0, 1, 50)
    shapes = np.full((len(pk) - 1, 50), np.nan)
    for i in range(len(pk) - 1):
        seg = x[pk[i]:pk[i + 1]]
        if len(seg) < 8:
            continue
        seg = (seg - seg.min()) / (seg.max() - seg.min() + 1e-12)
        shapes[i] = np.interp(grid, np.linspace(0, 1, len(seg)), seg)
    beat_t = pk[:-1] / FS
    early, mid = slice(0, 15), slice(20, 40)
    feat = np.nanmean(shapes[:, early], axis=1) / (np.nanmean(shapes[:, mid], axis=1) + 1e-9)
    ok = np.isfinite(feat)
    if ok.sum() < 5:
        return out
    for s_i in range(n_sec):
        m = ok & (np.abs(beat_t - s_i) <= 2.5)
        if m.sum() >= 3:
            out[s_i] = np.nanmedian(feat[m])
    return out


def _fill_nan(a):
    idx = np.arange(len(a))
    good = np.isfinite(a)
    if good.sum() < 2:
        return np.nan_to_num(a)
    return np.interp(idx, idx[good], a[good])


def perm_entropy(x, order=3, delay=1):
    """Normalised permutation entropy. O(N), cheap enough for a wearable."""
    n = len(x) - delay * (order - 1)
    if n < 10:
        return np.nan
    emb = np.empty((n, order))
    for i in range(order):
        emb[:, i] = x[i * delay: i * delay + n]
    perms = np.argsort(emb, axis=1, kind="stable")
    codes = np.zeros(n, dtype=np.int64)
    for i in range(order):
        codes = codes * order + perms[:, i]
    _, cnt = np.unique(codes, return_counts=True)
    p = cnt / cnt.sum()
    return float(-(p * np.log(p)).sum() / np.log(math.factorial(order)))


def sample_entropy(x, m=2, r_frac=0.2):
    """Sample entropy, vectorised. Expects a short (<=800 pt) window."""
    n = len(x)
    if n < m + 20:
        return np.nan
    r = r_frac * np.std(x)
    if r <= 0:
        return np.nan

    def _count(mm):
        k = n - mm + 1
        emb = np.lib.stride_tricks.sliding_window_view(x, mm)[:k]
        c = 0
        for i in range(k - 1):
            d = np.max(np.abs(emb[i + 1:] - emb[i]), axis=1)
            c += int((d <= r).sum())
        return c

    b, a = _count(m), _count(m + 1)
    if b == 0 or a == 0:
        return np.nan
    return float(-np.log(a / b))


def spectral_entropy(x, fs):
    f, p = sps.welch(x, fs=fs, nperseg=min(len(x), 256))
    p = p[f > 0]
    p = p / (p.sum() + 1e-15)
    return float(-(p * np.log(p)).sum() / np.log(len(p)))


def entropy_stream(x, fs, win_s=60, kind="perm", hop_s=1, dec=None):
    """
    Rolling entropy at 1 Hz, mapped to Temple's stated 1-250 scale at the end.
    x may be the 64 Hz waveform or the 1 Hz Brain Flow stream.
    """
    n_sec = int(len(x) / fs)
    w = int(win_s * fs)
    out = np.full(n_sec, np.nan)
    for s in range(0, n_sec, hop_s):
        i1 = int((s + 1) * fs)
        i0 = max(0, i1 - w)
        seg = x[i0:i1]
        if len(seg) < w * 0.6:
            continue
        if dec and dec > 1:
            seg = sps.decimate(seg, dec, ftype="fir", zero_phase=True)
        if kind == "perm":
            out[s] = perm_entropy(seg, order=4)
        elif kind == "samp":
            out[s] = sample_entropy(seg[-600:])
        else:
            out[s] = spectral_entropy(seg, fs / (dec or 1))
    if hop_s > 1:
        out = _fill_nan(out)
    return _fill_nan(out)


def beat_intervals(raw):
    """Detected beat times (s) and intervals (s) from the noisy sensor trace."""
    x = _bandpass(raw, 0.6, 6.0)
    pk, _ = sps.find_peaks(x, distance=int(0.34 * FS), prominence=0.30 * np.std(x))
    if len(pk) < 8:
        return np.array([]), np.array([])
    rri = np.diff(pk) / FS
    bt = pk[:-1] / FS
    keep = (rri > 0.28) & (rri < 1.6)
    return bt[keep], rri[keep]


def entropy_beats(raw, n_sec, n_beats=140, kind="samp"):
    """
    Rolling entropy of the last n_beats detected intervals, reported at 1 Hz.
    Beat detection runs on the noisy trace, so motion artifact corrupts the
    intervals exactly when workload is highest. That corruption is left in.
    """
    bt, rri = beat_intervals(raw)
    out = np.full(n_sec, np.nan)
    if len(rri) < n_beats // 2:
        return out
    j = 0
    for s_i in range(n_sec):
        while j < len(bt) and bt[j] <= s_i:
            j += 1
        if j < n_beats:
            continue
        seg = rri[j - n_beats:j]
        out[s_i] = sample_entropy(seg, m=2, r_frac=0.2) if kind == "samp" \
            else perm_entropy(seg, order=4)
    return _fill_nan(out)


def rri_series(raw, n_sec):
    """
    Beat-to-beat interval series recovered from the sensor trace, at 1 Hz.

    Beat detection runs on the noisy waveform, so motion artifact corrupts the
    intervals exactly when workload is highest. That corruption is part of the
    measurement situation and is left in.
    """
    x = _bandpass(raw, 0.6, 6.0)
    pk, _ = sps.find_peaks(x, distance=int(0.34 * FS), prominence=0.30 * np.std(x))
    if len(pk) < 8:
        return np.full(n_sec, np.nan)
    rri = np.diff(pk) / FS
    tb_ = pk[:-1] / FS
    keep = (rri > 0.28) & (rri < 1.6)
    if keep.sum() < 8:
        return np.full(n_sec, np.nan)
    return np.interp(np.arange(n_sec), tb_[keep], rri[keep])


def entropy_rri(rri_1hz, win_s=120, kind="samp"):
    """Rolling entropy of the recovered beat-interval series."""
    n = len(rri_1hz)
    out = np.full(n, np.nan)
    for s_i in range(n):
        i0 = max(0, s_i - win_s)
        seg = rri_1hz[i0:s_i + 1]
        if len(seg) < win_s * 0.5 or not np.all(np.isfinite(seg)):
            continue
        out[s_i] = sample_entropy(seg, m=2, r_frac=0.2) if kind == "samp" \
            else perm_entropy(seg, order=4)
    return _fill_nan(out)


def to_temple_scale(e):
    """Map a raw entropy series onto the reported 1-250 range, cohort-wise."""
    lo, hi = np.nanpercentile(e, 1), np.nanpercentile(e, 99)
    if not np.isfinite(lo) or hi <= lo:
        return np.full_like(e, np.nan)
    return 1 + 249 * np.clip((e - lo) / (hi - lo), 0, 1)


# ----------------------------------------------------------------------------
# reference channels
# ----------------------------------------------------------------------------
def references(sub, phys, rng):
    """1 Hz TCD MCAv, breath-by-breath VO2 resampled to 1 Hz, 1 Hz ECG HR."""
    n_sec = int(len(phys["t"]) / FS)
    sl = (np.arange(n_sec) * int(FS)).astype(int)

    # --- TCD: per-subject probe gain, slow insonation drift, noise, dropout
    true_mcav = np.array([phys["mcav"][i:i + int(FS)].mean() for i in sl])
    drift = np.cumsum(rng.normal(0, _p("tcd_drift_sd", 0.0015), n_sec))
    drift = drift - drift.mean()
    tcd = true_mcav * sub["probe_gain"] * (1 + drift) + rng.normal(0, _p("tcd_noise", 0.9), n_sec)
    drop = rng.random(n_sec) < _p("tcd_dropout", 0.02)
    tcd[drop] = np.nan

    # --- metabolic cart: breath-by-breath, then held at 1 Hz
    true_vo2 = np.array([phys["vo2"][i:i + int(FS)].mean() for i in sl])
    rr_sec = np.array([phys["rr"][i:i + int(FS)].mean() for i in sl])
    tb, vb = [], []
    tt = 0.0
    while tt < n_sec - 1:
        j = int(tt)
        tb.append(tt)
        vb.append(true_vo2[j] * (1 + rng.normal(0, 0.09)))
        tt += 1.0 / max(rr_sec[j], 0.12)
    vo2_bxb = np.interp(np.arange(n_sec), tb, vb)

    true_hr = np.array([phys["hr"][i:i + int(FS)].mean() for i in sl])
    ecg_hr = true_hr + rng.normal(0, 0.6, n_sec)

    return dict(tcd=tcd, true_mcav=true_mcav, vo2_bxb=vo2_bxb,
                true_vo2=true_vo2, ecg_hr=ecg_hr, true_hr=true_hr)
