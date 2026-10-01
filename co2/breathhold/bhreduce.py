"""
Breath-hold hypercapnia, optical arm.

EuskalIBUR (OpenNeuro ds003192), physio only. Columns at 10 kHz:
time, respiratory_chest, trigger, cardiac, respiratory_CO2, respiratory_O2.
`cardiac` is the scanner pulse oximeter: a peripheral OPTICAL signal.

No epochs and no task timing are used. The measured capnogram is the regressor,
breath by breath, which is the same move the whole project makes: regress the
response on the stimulus that was actually delivered, never on the block
structure that was intended.
"""
import gzip, numpy as np
from scipy.signal import find_peaks, butter, filtfilt
FS_RAW=10000.0; DS=20; FS=FS_RAW/DS      # 500 Hz

def load(path):
    out=[]
    with gzip.open(path,'rt') as f:
        for i,l in enumerate(f):
            if i % DS: continue
            p=l.split('\t')
            out.append((float(p[0]),float(p[1]),float(p[3]),float(p[4]),float(p[5])))
    X=np.array(out)
    return dict(t=X[:,0],resp=X[:,1],ppg=X[:,2],co2=X[:,3],o2=X[:,4])

def breaths(co2):
    pk,_=find_peaks(co2,distance=int(1.5*FS),prominence=np.ptp(co2)*0.15)
    tr,_=find_peaks(-co2,distance=int(1.5*FS),prominence=np.ptp(co2)*0.15)
    return pk,tr

def beats(ppg):
    b,a=butter(3,[0.5/(FS/2),8/(FS/2)],'bandpass')
    y=filtfilt(b,a,np.nan_to_num(ppg))
    pk,_=find_peaks(y,distance=int(0.33*FS),height=np.percentile(y,75))
    return pk,y

def features(d):
    """One row per breath: end-tidal CO2 and PPG features over the interval
    from that breath's end-tidal point to the next."""
    pk,tr=breaths(d['co2'])
    bp,filt=beats(d['ppg'])
    ppg=d['ppg']
    rows=[]
    for i in range(len(pk)-1):
        a,b=pk[i],pk[i+1]
        if b-a < int(1.0*FS) or b-a > int(60*FS): continue
        bb=bp[(bp>=a)&(bp<b)]
        if len(bb)<3: continue
        # per-beat AC amplitude and DC level of the raw PPG
        acs=[];dcs=[];sys=[];dia=[]
        for j in range(len(bb)-1):
            s=ppg[bb[j]:bb[j+1]]
            if s.size<int(0.3*FS): continue
            acs.append(s.max()-s.min()); dcs.append(s.mean())
            sys.append(s.max()); dia.append(s.min())
        if len(acs)<2: continue
        hr=60.0*FS/np.median(np.diff(bb))
        rows.append(dict(t=float(d['t'][a]),
                         et=float(d['co2'][a]),
                         o2=float(np.percentile(d['o2'][a:b],5)),
                         ac=float(np.median(acs)), dc=float(np.median(dcs)),
                         sys=float(np.median(sys)), dia=float(np.median(dia)),
                         pi=float(np.median(acs))/abs(float(np.median(dcs))) if np.median(dcs)!=0 else np.nan,
                         hr=hr, ibi=float(b-a)/FS,
                         respsd=float(np.std(d['resp'][a:b]))))
    return rows
