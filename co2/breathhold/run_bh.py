"""Download, reduce, delete. One file held at a time. Resumable."""
import os,sys,json,time,subprocess
import numpy as np, bhreduce as B
S3="https://s3.amazonaws.com/openneuro.org/ds003192/sub-{s}/ses-{e}/func/sub-{s}_ses-{e}_task-breathhold_physio.tsv.gz"
SUBS=['001','002','003','004','007','008','009']
SESS=[f"{i:02d}" for i in range(1,11)]
OUT='bh_rows.jsonl'
DEADLINE=float(sys.argv[1]) if len(sys.argv)>1 else 480.0
NSES=int(sys.argv[2]) if len(sys.argv)>2 else 10

def main():
    done=set()
    if os.path.exists(OUT):
        done={(json.loads(l)['sub'],json.loads(l)['ses']) for l in open(OUT) if l.strip()}
    t0=time.time()
    for e in SESS[:NSES]:
        for s in SUBS:
            if (s,e) in done: continue
            if time.time()-t0>DEADLINE: print('deadline',flush=True); return
            url=S3.format(s=s,e=e)
            r=subprocess.run(['curl','-sSLf','-o','w.gz',url],capture_output=True)
            if r.returncode!=0 or not os.path.exists('w.gz'):
                print(f"  {s}/{e} no file",flush=True); 
                with open(OUT,'a') as f: f.write(json.dumps(dict(sub=s,ses=e,rows=[]))+'\n')
                continue
            try:
                d=B.load('w.gz'); rows=B.features(d)
            except Exception as ex:
                print(f"  {s}/{e} FAIL {type(ex).__name__}",flush=True); rows=[]
            os.remove('w.gz')
            with open(OUT,'a') as f:
                f.write(json.dumps(dict(sub=s,ses=e,rows=rows))+'\n')
            et=[r_['et'] for r_ in rows]
            print(f"  {s}/{e} {len(rows):3d} breaths  ET {min(et):.3f}-{max(et):.3f}V"
                  if rows else f"  {s}/{e} 0 breaths", flush=True)
    print('ALL DONE',flush=True)
main()
