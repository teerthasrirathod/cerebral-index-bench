# EXPLORATORY, post-hoc (not in prereg): how much extra device noise in squat/supine
# reproduces their lag-adjusted r, and does that then predict their zero-lag r and BA?
import json
from sim import cohort, SEED
out={}
for f in (0.75,0.9):
    N=30; base={0.75:1.362,0.9:1.324}[f]
    for p,target in (("squat",0.843),("supine",0.839)):
        lo,hi=base,base*4
        for _ in range(12):
            m=(lo+hi)/2
            c=cohort(p,"VOLUME",f,N,m,SEED+7,reps=4)
            lo,hi=(m,hi) if c["r_lag"]>target else (lo,m)
        sd=(lo+hi)/2; c=cohort(p,"VOLUME",f,N,sd,SEED+7,reps=8)
        out[f"{p}_f{f}"]=dict(sd=round(sd,3),noise_x=round(sd/base,2),**c)
        print(p,f,json.dumps(out[f"{p}_f{f}"]),flush=True)
json.dump(out,open("explore_noise.json","w"),indent=1)
