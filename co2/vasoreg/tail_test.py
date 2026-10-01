"""
Where inside a window the epoch is taken is a free parameter, and for a manoeuvre
whose stimulus RISES through the window it is not a neutral one.

run_all.py capped every epoch at the first 120 true seconds. Rebreathing windows
run 167-425 s. Taking the first 120 s of a CO2 wash-in samples the ramp, not the
plateau: it understates the stimulus, which inflates any reactivity computed as
response/stimulus and makes subjects fail an end-tidal gate they actually passed.

Prediction: re-reducing the SAME epochs from the END of the window raises dET,
raises the pass rate, and moves hypercapnic reactivity down toward the published
3-4 %/mmHg band. Falsifiable and cheap.
"""
import json, numpy as np, wfdb
from reduce import reduce_epoch, FS_TRUE
PN='cerebral-vasoreg-diabetes/1.0.0/Data/Labview/Converted/Head-up-tilt_Day2'
CH=[1,2,3,4,5,9]; IDX=dict(ecg=0,abp=1,mcar=2,mcal=3,radi=4,co2=5)
CAP=int(120*FS_TRUE); TRIM=0.20
R=[json.loads(l) for l in open('epochs_all.jsonl')]
out=[]
for d in R:
    row=dict(sid=d['sid'],group=d['group'],base_et=d['epochs'].get('baseline',{}).get('co2_et'))
    for m in ('hypervent','rebreath'):
        e=d['epochs'].get(m)
        if not e or e.get('src')!='marker': continue
        span=d.get('span_s',{}).get(m)
        if span is None: continue
        a=e['s0']; b=a+int(span*(1-2*TRIM)*FS_TRUE)
        if b-a<=CAP:
            row[f'{m}_head']=e['co2_et']; row[f'{m}_tail']=e['co2_et']
            row[f'{m}_same']=True; row[f'{m}_tailmfv']=None
            continue
        try:
            r=wfdb.rdrecord('s'+d['sid'][1:]+'DB',pn_dir=PN,channels=CH,sampfrom=b-CAP,sampto=b)
        except Exception as ex:
            print(d['sid'],m,'FAIL',type(ex).__name__,flush=True); continue
        t=reduce_epoch(r.p_signal.astype(float),IDX)
        row[f'{m}_head']=e['co2_et']; row[f'{m}_tail']=t['co2_et']; row[f'{m}_same']=False
        v=[t.get('mcar_mfv'),t.get('mcal_mfv')]; v=[x for x in v if x and np.isfinite(x)]
        row[f'{m}_tailmfv']=float(np.mean(v)) if v else None
        v=[e.get('mcar_mfv'),e.get('mcal_mfv')]; v=[x for x in v if x and np.isfinite(x)]
        row[f'{m}_headmfv']=float(np.mean(v)) if v else None
    out.append(row)
    print(f"  {d['sid']} rebreath head_ET={row.get('rebreath_head')} tail_ET={row.get('rebreath_tail')}",flush=True)
json.dump(out,open('tail.json','w'))
print('DONE')
