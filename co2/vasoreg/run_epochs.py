"""
Re-reduce every epoch through reduce.py. Resumable at (subject, epoch) grain.

Window POSITIONS are taken from cohort.jsonl, which recorded them as
header-seconds; they are converted straight back to samples (t * FS_HDR) and
never re-divided by a sampling rate. The windows are marker-anchored with a
proportional trim, so they land on the correct physiological intervals whatever
the header's fs says. Only their printed seconds were ever wrong.

Declared reduction: EPOCH_MAX_SAMPLES caps how much of a long epoch is read, to
keep each network fetch inside the shell's time limit. Recorded per row as
`capped` so its effect is visible rather than silent.
"""
import os, sys, json, time
import numpy as np
import wfdb
from reduce import reduce_epoch, FS_HDR, FS_TRUE

PN  = 'cerebral-vasoreg-diabetes/1.0.0/Data/Labview/Converted/Head-up-tilt_Day2'
CH  = [1, 2, 3, 4, 5, 9]
IDX = dict(ecg=0, abp=1, mcar=2, mcal=3, radi=4, co2=5)
EPOCH_MAX_SAMPLES = 60_000           # 120 s of true time, uniform across epochs
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'epochs2.jsonl')
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cohort.jsonl')
DEADLINE = float(sys.argv[1]) if len(sys.argv) > 1 else 150.0


def worklist():
    w = []
    for line in open(SRC):
        d = json.loads(line)
        for name, e in d['epochs'].items():
            a = int(round(e['t0'] * FS_HDR)); b = int(round(e['t1'] * FS_HDR))
            w.append((d['sid'], d['group'], name, a, b))
    return w


def main():
    done = set()
    if os.path.exists(OUT):
        done = {(json.loads(l)['sid'], json.loads(l)['epoch'])
                for l in open(OUT) if l.strip()}
    t0 = time.time()
    todo = [x for x in worklist() if (x[0], x[2]) not in done]
    print(f"{len(done)} done, {len(todo)} to go", flush=True)
    for sid, grp, name, a, b in todo:
        if time.time() - t0 > DEADLINE:
            print("deadline, stopping cleanly", flush=True); break
        capped = (b - a) > EPOCH_MAX_SAMPLES
        bb = a + EPOCH_MAX_SAMPLES if capped else b
        try:
            r = wfdb.rdrecord('s' + sid[1:] + 'DB', pn_dir=PN,
                              channels=CH, sampfrom=a, sampto=bb)
        except Exception as ex:
            print(f"  {sid} {name} FETCH-FAIL {type(ex).__name__}", flush=True); continue
        row = dict(sid=sid, group=grp, epoch=name, s0=a, s1=bb, capped=capped)
        row.update(reduce_epoch(r.p_signal.astype(float), IDX))
        with open(OUT, 'a') as f:
            f.write(json.dumps(row) + '\n')
        print(f"  {sid:6s} {name:10s} hr={row['hr']:6.1f} "
              f"co2et={row['co2_et']:6.2f} mcarPI={row['mcar_pi']:5.3f} "
              f"radiPI={row['radi_pi']:5.3f} ({time.time()-t0:.0f}s)", flush=True)
    print("EXIT", flush=True)


if __name__ == '__main__':
    main()
