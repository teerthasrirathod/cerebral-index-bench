"""
Cohort run: cerebral versus peripheral response, per manoeuvre, on real recordings.

Novak V, Mendez L. Cerebral Vasoregulation in Diabetes. PhysioNet 1.0.0,
doi 10.13026/m40k-4758. Open Access under CC BY 4.0.

For every subject with usable Day-2 markers, each manoeuvre is located from the
marker pulse train, then VALIDATED against the capnograph: a hyperventilation
epoch must show lower CO2 than baseline and a rebreathing epoch must show higher.
Epochs that fail that check are recorded as failures rather than averaged in,
because the marker file itself documents shifted and cancelled events for several
subjects.
"""
import os, sys, json, csv, traceback
import numpy as np
import time
import wfdb


def rdrecord_retry(rec, pn_dir, tries=5):
    """The egress tunnel drops occasionally mid-transfer and surfaces as a bare
    connection reset. Retry with backoff rather than losing the subject."""
    last = None
    for k in range(tries):
        try:
            return wfdb.rdrecord(rec, pn_dir=pn_dir)
        except Exception as e:
            last = e
            time.sleep(4 * (k + 1))
    raise last

PN = 'cerebral-vasoreg-diabetes/1.0.0/Data/Labview/Converted/Head-up-tilt_Day2'
OUT = 'cohort.jsonl'
CH = dict(marker=0, ecg=1, abp=2, mcar=3, mcal=4, radi=5, o2=8, co2=9)
TRIM = 0.20          # drop the first and last 20% of each marked interval


def marker_times(mk, fs):
    """Ordinal event times from the marker pulse train, deduplicated."""
    d = np.diff(np.round(np.nan_to_num(mk)))
    ups = np.where(d > 0)[0] + 1
    keep = [ups[0]] if len(ups) else []
    for i in ups[1:]:
        if (i - keep[-1]) / fs > 1.0:      # merge pulses within 1 s
            keep.append(i)
    return np.array(keep) / fs


def seg(sig, fs, t0, t1, ch):
    a, b = int(t0 * fs), int(t1 * fs)
    if b <= a or a < 0 or b > len(sig):
        return np.nan
    x = sig[a:b, CH[ch]]
    x = x[np.isfinite(x)]
    return float(np.mean(x)) if x.size > 100 else np.nan


def hr_from_abp(sig, fs, t0, t1):
    from scipy.signal import find_peaks
    a, b = int(t0 * fs), int(t1 * fs)
    x = sig[a:b, CH['abp']]
    x = x[np.isfinite(x)]
    if x.size < fs * 10:
        return np.nan
    pk, _ = find_peaks(x, distance=int(0.4 * fs), prominence=np.std(x) * 0.5)
    if len(pk) < 5:
        return np.nan
    return float(60.0 / np.median(np.diff(pk) / fs))


def window(times, i_start, i_end, dur_default=180.0):
    """Trimmed window between two ordinal markers (1-based indices)."""
    if i_start <= 0 or i_start > len(times):
        return None
    t0 = times[i_start - 1]
    t1 = times[i_end - 1] if (0 < i_end <= len(times)) else t0 + dur_default
    if t1 <= t0 + 20:
        return None
    pad = (t1 - t0) * TRIM
    return (t0 + pad, t1 - pad)


def run_subject(row, done):
    sid = row[0]
    rec = 's' + sid[1:] + 'DB'
    if sid in done:
        return None
    r = rdrecord_retry(rec, PN)
    sig, fs = r.p_signal, r.fs
    t = marker_times(sig[:, CH['marker']], fs)

    # positional: the header repeats "End" seven times, so names cannot be used
    # 2 Start10  3 Val1  4 End  5 Val2  6 End  7 Val3  8 End  9 Val4 10 End
    # 11 Hyper  12 End  13 Hypo  14 End  15 Tilt  16 End  17 Hypo2  18 Tilt2  19 EndTilt
    C = [int(x) if str(x).strip().isdigit() else 0 for x in row[2:20]]
    g = lambda k: C[k - 2]
    wins = {}
    b0, v1 = g(2), g(3)
    if 0 < b0 <= len(t):
        end_b = t[v1 - 1] - 10 if 0 < v1 <= len(t) else t[b0 - 1] + 280
        wins['baseline'] = (t[b0 - 1] + 40, min(t[b0 - 1] + 280, end_b))
    wins['hypervent'] = window(t, g(11), g(12))
    wins['rebreath'] = window(t, g(13), g(14))
    wins['tilt'] = window(t, g(15), g(16))

    out = dict(sid=sid, group=row[1], n_markers=len(t), fs=fs,
               dur_min=round(len(sig) / fs / 60, 1), epochs={}, flags=[])
    for name, w in wins.items():
        if w is None or w[1] <= w[0]:
            out['flags'].append(f'{name}:no_window'); continue
        e = {c: seg(sig, fs, w[0], w[1], c) for c in ('mcar','mcal','radi','abp','co2')}
        e['hr'] = hr_from_abp(sig, fs, w[0], w[1])
        e['t0'], e['t1'] = round(w[0], 1), round(w[1], 1)
        out['epochs'][name] = e

    # capnograph validation
    E = out['epochs']
    if 'baseline' in E and 'hypervent' in E:
        if not (E['hypervent']['co2'] < E['baseline']['co2'] - 1):
            out['flags'].append('hypervent:co2_not_lower')
    if 'baseline' in E and 'rebreath' in E:
        if not (E['rebreath']['co2'] > E['baseline']['co2'] + 3):
            out['flags'].append('rebreath:co2_not_higher')
    return out


def main():
    done = set()
    if os.path.exists(OUT):
        done = {json.loads(l)['sid'] for l in open(OUT) if l.strip()}
    with open('markers_day2.csv') as f:
        rows = [r for r in csv.reader(f)][1:]
    for i, row in enumerate(rows, 1):
        sid = row[0]
        try:
            res = run_subject(row, done)
            if res is None:
                print(f"[{i}/{len(rows)}] {sid} cached", flush=True); continue
            with open(OUT, 'a') as f:
                f.write(json.dumps(res) + '\n')
            print(f"[{i}/{len(rows)}] {sid} {res['group']:8s} "
                  f"{len(res['epochs'])} epochs  flags={res['flags']}", flush=True)
        except Exception as ex:
            print(f"[{i}/{len(rows)}] {sid} FAILED {type(ex).__name__}: {ex}", flush=True)
    print("DONE", flush=True)


if __name__ == '__main__':
    main()
