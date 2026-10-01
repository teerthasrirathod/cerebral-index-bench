"""Full cohort, cloud container. Resumable per subject. All work in SAMPLES."""
import os,sys,json,csv,time
import numpy as np, wfdb
from reduce import reduce_epoch, FS_TRUE, FS_HDR
PN='cerebral-vasoreg-diabetes/1.0.0/Data/Labview/Converted/Head-up-tilt_Day2'
IDX=dict(marker=0,ecg=1,abp=2,mcar=3,mcal=4,radi=5,co2=9)
RED=dict(ecg=1,abp=2,mcar=3,mcal=4,radi=5,co2=9)
OUT='epochs_all.jsonl'
TRIM=0.20
CAP=int(120*FS_TRUE)              # 120 true s per epoch, uniform
MERGE=int(2.0*FS_TRUE)            # markers closer than 2 true s are one event
MINWIN=int(20*FS_TRUE)
DOC={'hypervent':180.0,'rebreath':180.0,'tilt':300.0}   # documented protocol durations, true s
DEADLINE=float(sys.argv[1]) if len(sys.argv)>1 else 520.0

def marker_samples(mk):
    d=np.diff(np.round(np.nan_to_num(mk))); ups=np.where(d>0)[0]+1
    keep=[]
    for i in ups:
        if not keep or i-keep[-1]>MERGE: keep.append(int(i))
    return np.array(keep)

def win(t,i0,i1,default_s):
    if i0<=0 or i0>len(t): return None,'no_start'
    a=t[i0-1]
    if 0<i1<=len(t):
        b=t[i1-1]; d='marker'
    else:
        b=a+int(default_s*FS_TRUE); d='default'
    if b-a<MINWIN: return None,'too_short'
    pad=int((b-a)*TRIM)
    return (a+pad,b-pad),d

def main():
    done=set()
    if os.path.exists(OUT):
        done={json.loads(l)['sid'] for l in open(OUT) if l.strip()}
    rows=[r for r in csv.reader(open('markers_day2.csv'))][1:]
    t0=time.time()
    for n,row in enumerate(rows,1):
        sid=row[0]
        if sid in done: continue
        if time.time()-t0>DEADLINE:
            print('deadline',flush=True); break
        try:
            r=wfdb.rdrecord('s'+sid[1:]+'DB',pn_dir=PN)
        except Exception as e:
            print(f'[{n}] {sid} FETCH-FAIL {type(e).__name__}',flush=True); continue
        sig=r.p_signal
        t=marker_samples(sig[:,IDX['marker']])
        C=[int(x) if str(x).strip().isdigit() else 0 for x in row[2:20]]
        g=lambda k:C[k-2]
        specs={}
        b0,v1=g(2),g(3)
        if 0<b0<=len(t):
            a=t[b0-1]+int(60*FS_TRUE)
            b=t[v1-1] if 0<v1<=len(t) else a+CAP
            if b-a>=MINWIN: specs['baseline']=((a,min(a+CAP,b)),'marker')
        for nm,(i0,i1) in (('hypervent',(11,12)),('rebreath',(13,14)),('tilt',(15,16))):
            w,d=win(t,g(i0),g(i1),DOC[nm])
            if w: specs[nm]=(w,d)
        out=dict(sid=sid,group=row[1],n_markers=int(len(t)),max_ord=max(C),
                 rec_min=round(len(sig)/FS_TRUE/60,1),epochs={},flags=[])
        for nm,(w,d) in specs.items():
            a,b=w
            # protocol-duration ruler: full marker span against the documented length
            if d=='marker' and nm in DOC:
                span=(b-a)/(1-2*TRIM)/FS_TRUE
                out.setdefault('span_s',{})[nm]=round(span,1)
                if abs(span-DOC[nm])/DOC[nm]>0.25: out['flags'].append(f'{nm}:span_{span:.0f}s')
            if d=='default': out['flags'].append(f'{nm}:default_window')
            blk=sig[a:min(b,a+CAP),[RED[c] for c in ('ecg','abp','mcar','mcal','radi','co2')]]
            e=reduce_epoch(blk.astype(float),dict(ecg=0,abp=1,mcar=2,mcal=3,radi=4,co2=5))
            e['s0'],e['s1']=int(a),int(min(b,a+CAP)); e['src']=d
            out['epochs'][nm]=e
        E=out['epochs']
        if 'baseline' in E:
            if 'hypervent' in E and not (E['hypervent']['co2_et']<E['baseline']['co2_et']-3):
                out['flags'].append('hypervent:et_not_lower')
            if 'rebreath' in E and not (E['rebreath']['co2_et']>E['baseline']['co2_et']+3):
                out['flags'].append('rebreath:et_not_higher')
        with open(OUT,'a') as f: f.write(json.dumps(out)+'\n')
        print(f"[{n}/31] {sid} {out['group']:8s} {len(E)}ep mk={len(t)}/{max(C)} "
              f"{out['rec_min']}min spans={out.get('span_s',{})} flags={out['flags']} "
              f"({time.time()-t0:.0f}s)",flush=True)
    print('EXIT',flush=True)
main()
