"""
Sensitivity of the headline number to the reductions nobody declared.

Four parameters, swept against the within-subject CO2 fit (2.75 %/mmHg, R2 0.82):

  TRIM          fraction cut from each end of a marker-to-marker window   [0.00 .. 0.40]
  baseline      where inside the 10 min rest the baseline epoch sits      skip x duration
  MERGE         marker pulses closer than this are treated as one event   [0.5 .. 4.0 s]
  referencing   every manoeuvre against ONE fixed baseline, or against
                the rest immediately preceding it

Each record is read ONCE. R-peaks are detected once over a generous superset
span per manoeuvre and then subset per TRIM, so 360 settings cost one pass.
"""
import os, sys, json, csv, time
import numpy as np, wfdb
from reduce import rpeaks, beat_stats, capno, FS_TRUE

PN='cerebral-vasoreg-diabetes/1.0.0/Data/Labview/Converted/Head-up-tilt_Day2'
IDX=dict(marker=0,ecg=1,abp=2,mcar=3,mcal=4,radi=5,co2=9)
OUT='sweep_raw.jsonl'
DEADLINE=float(sys.argv[1]) if len(sys.argv)>1 else 450.0

TRIMS   =[0.00,0.10,0.20,0.30,0.40]
MERGES  =[0.5,1.0,2.0,4.0]
BASE_SKIP=[30,60,120]
BASE_DUR =[60,120,240]
DOC={'hypervent':180.0,'rebreath':180.0}
MINWIN=int(20*FS_TRUE)

def marker_samples(mk,merge_s):
    d=np.diff(np.round(np.nan_to_num(mk))); ups=np.where(d>0)[0]+1
    m=int(merge_s*FS_TRUE); keep=[]
    for i in ups:
        if not keep or i-keep[-1]>m: keep.append(int(i))
    return np.array(keep)

def et_mfv(block, r, a, b):
    """co2 p95 and cerebral MFV over [a,b) of the span, reusing precomputed beats."""
    seg=block[a:b]
    rr=r[(r>=a)&(r<b)]-a
    c=capno(seg[:,5])
    v=[]
    for ch in (2,3):
        s=beat_stats(seg[:,ch],rr)
        if np.isfinite(s['mfv']): v.append(s['mfv'])
    return c['et'], (float(np.mean(v)) if v else np.nan)

def main():
    done=set()
    if os.path.exists(OUT):
        done={json.loads(l)['sid'] for l in open(OUT) if l.strip()}
    rows=[r for r in csv.reader(open('markers_day2.csv'))][1:]
    t0=time.time()
    for n,row in enumerate(rows,1):
        sid=row[0]
        if sid in done: continue
        if time.time()-t0>DEADLINE: print('deadline',flush=True); break
        try: rec=wfdb.rdrecord('s'+sid[1:]+'DB',pn_dir=PN)
        except Exception as e:
            print(f'[{n}] {sid} FETCH-FAIL {type(e).__name__}',flush=True); continue
        sig=rec.p_signal
        C=[int(x) if str(x).strip().isdigit() else 0 for x in row[2:20]]
        g=lambda k:C[k-2]
        out=dict(sid=sid,group=row[1],settings=[])
        for merge in MERGES:
            t=marker_samples(sig[:,IDX['marker']],merge)
            spans={}
            for nm,(i0,i1) in (('hypervent',(11,12)),('rebreath',(13,14))):
                a_,b_=g(i0),g(i1)
                if not (0<a_<=len(t)): continue
                A=t[a_-1]
                B=t[b_-1] if 0<b_<=len(t) else A+int(DOC[nm]*FS_TRUE)
                if B-A<MINWIN: continue
                spans[nm]=(int(A),int(B))
            b0,v1=g(2),g(3)
            rest=None
            if 0<b0<=len(t):
                R0=t[b0-1]; R1=t[v1-1] if 0<v1<=len(t) else R0+int(600*FS_TRUE)
                rest=(int(R0),int(R1))
            if not spans or rest is None: continue
            # one read of channels per span, one beat detection per span
            cache={}
            for nm,(A,B) in spans.items():
                blk=sig[A:B,[IDX[c] for c in ('ecg','abp','mcar','mcal','radi','co2')]].astype(float)
                cache[nm]=(blk, rpeaks(blk[:,0]), A, B)
            rblk=sig[rest[0]:rest[1],[IDX[c] for c in ('ecg','abp','mcar','mcal','radi','co2')]].astype(float)
            rpk=rpeaks(rblk[:,0])
            # preceding-rest windows: 120 s ending 20 s before each manoeuvre
            pre={}
            for nm,(A,B) in spans.items():
                pb=A-int(20*FS_TRUE); pa=pb-int(120*FS_TRUE)
                if pa<0: continue
                blk=sig[pa:pb,[IDX[c] for c in ('ecg','abp','mcar','mcal','radi','co2')]].astype(float)
                pre[nm]=(blk, rpeaks(blk[:,0]))
            for trim in TRIMS:
                for skip in BASE_SKIP:
                    for dur in BASE_DUR:
                        ba=int(skip*FS_TRUE); bb=ba+int(dur*FS_TRUE)
                        if bb>rblk.shape[0]: continue
                        bet,bmfv=et_mfv(rblk,rpk,ba,bb)
                        if not (np.isfinite(bet) and np.isfinite(bmfv)): continue
                        for nm,(blk,pk,A,B) in cache.items():
                            L=B-A; pad=int(L*trim)
                            if L-2*pad<MINWIN: continue
                            et,mfv=et_mfv(blk,pk,pad,L-pad)
                            if not (np.isfinite(et) and np.isfinite(mfv)): continue
                            rec_={'merge':merge,'trim':trim,'skip':skip,'dur':dur,'man':nm,
                                  'span_s':round(L/FS_TRUE,1),
                                  'fixed':{'det':et-bet,'dmfv':100*(mfv-bmfv)/bmfv}}
                            if nm in pre:
                                pblk,ppk=pre[nm]
                                pet,pmfv=et_mfv(pblk,ppk,0,pblk.shape[0])
                                if np.isfinite(pet) and np.isfinite(pmfv) and pmfv>0:
                                    rec_['prev']={'det':et-pet,'dmfv':100*(mfv-pmfv)/pmfv}
                            out['settings'].append(rec_)
        with open(OUT,'a') as f: f.write(json.dumps(out)+'\n')
        print(f"[{n}/31] {sid} {len(out['settings'])} setting-rows ({time.time()-t0:.0f}s)",flush=True)
    print('EXIT',flush=True)
main()
