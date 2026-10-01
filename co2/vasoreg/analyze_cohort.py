"""Cohort analysis. Every number in FINDINGS-21SEP-cohort.md comes from this file."""
import json, numpy as np, sys
R=[json.loads(l) for l in open('epochs_all.jsonl')]
DOC={'hypervent':180.0,'tilt':300.0}
def cb(e,k='mfv'):
    v=[e.get(f'mcar_{k}'),e.get(f'mcal_{k}')]
    v=[x for x in v if x is not None and np.isfinite(x)]
    return float(np.mean(v)) if v else np.nan
def pct(a,b): return 100*(a-b)/b if (np.isfinite(a) and np.isfinite(b) and b) else np.nan
def valid(d,m):
    e=d['epochs'].get(m); b=d['epochs'].get('baseline')
    if not e or not b: return False,'no_epoch'
    if e.get('src')!='marker': return False,'default_window'
    if m in DOC:
        sp=d.get('span_s',{}).get(m)
        if sp is None or abs(sp-DOC[m])/DOC[m]>0.25: return False,f'span_{sp}'
    if m=='hypervent' and not (e['co2_et']<b['co2_et']-3): return False,'et_not_lower'
    if m=='rebreath'  and not (e['co2_et']>b['co2_et']+3): return False,'et_not_higher'
    return True,'ok'
def ols(y,X,names):
    X=np.column_stack([np.ones(len(y))]+X)
    beta,_,_,_=np.linalg.lstsq(X,y,rcond=None); yh=X@beta
    rs=1-((y-yh)**2).sum()/((y-y.mean())**2).sum()
    n,k=X.shape; s2=((y-yh)**2).sum()/(n-k)
    se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return rs,beta,se
rec=[]
for d in R:
    b=d['epochs'].get('baseline')
    if not b: continue
    for m in ('hypervent','rebreath'):
        e=d['epochs'].get(m)
        if not e or e.get('src')!='marker': continue
        r=dict(sid=d['sid'],grp=d['group'],man=m,dco2=e['co2_et']-b['co2_et'],
               dmfv=pct(cb(e),cb(b)),dabp=pct(e['abp_mean'],b['abp_mean']),
               dhr=pct(e['hr'],b['hr']))
        if all(np.isfinite(list(r.values())[3:])): rec.append(r)
bysub={}
for r in rec: bysub.setdefault(r['sid'],[]).append(r)
pair={k:v for k,v in bysub.items() if len(v)>=2}
yy=[];cc=[];aa=[];hh=[]
for k,v in pair.items():
    for key,acc in (('dmfv',yy),('dco2',cc),('dabp',aa),('dhr',hh)):
        m=np.mean([x[key] for x in v]); acc.extend([x[key]-m for x in v])
yy,cc,aa,hh=map(np.array,(yy,cc,aa,hh))
r1,b1,s1=ols(yy,[cc],['dCO2'])
r2,_,_=ols(yy,[cc,aa],['dCO2','dABP'])
r3,_,_=ols(yy,[cc,aa,hh],['dCO2','dABP','dHR'])
print(f"within-subject CO2 reactivity {b1[1]:.2f} %/mmHg  95% CI "
      f"{b1[1]-1.96*s1[1]:.2f}-{b1[1]+1.96*s1[1]:.2f}  R2 {r1:.3f}  "
      f"n {len(yy)} epochs / {len(pair)} subjects")
print(f"incremental R2  dABP {r2-r1:+.4f}   dHR {r3-r2:+.4f}")
