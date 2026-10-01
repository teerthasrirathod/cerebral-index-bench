import numpy as np
from scipy.signal import find_peaks, butter, filtfilt, welch
FS_HDR, FS_TRUE = 1000.0, 500.0
MAX_HR_BPM, MAX_RESP_BPM = 220.0, 60.0
_MIN_RR = int(FS_TRUE*60.0/MAX_HR_BPM)
_MIN_BR = int(FS_TRUE*60.0/MAX_RESP_BPM)
def _hz(f): return f/(FS_TRUE/2.0)
def rpeaks(ecg):
    x=np.nan_to_num(np.asarray(ecg,float))
    if x.size < 5*FS_TRUE: return np.array([],int)
    b,a=butter(3,[_hz(5),_hz(25)],'bandpass'); y=np.abs(filtfilt(b,a,x))
    b,a=butter(3,_hz(8),'low'); y=filtfilt(b,a,y)
    pk,_=find_peaks(y,distance=_MIN_RR,height=np.percentile(y,98)*0.35)
    return pk
def beat_stats(sig,r):
    x=np.asarray(sig,float); ps=[];ed=[];pi=[]
    for i in range(len(r)-1):
        s=x[r[i]:r[i+1]]; s=s[np.isfinite(s)]
        if s.size<_MIN_RR//2: continue
        m=float(s.mean())
        if m<=0: continue
        ps.append(float(s.max())); ed.append(float(s.min())); pi.append((s.max()-s.min())/m)
    mean=float(np.nanmean(x)) if np.isfinite(x).any() else np.nan
    if len(pi)<5: return dict(mean=mean,mfv=np.nan,pi=np.nan,ps=np.nan,ed=np.nan,nbeat=len(pi))
    return dict(mean=mean,mfv=(np.median(ps)+2*np.median(ed))/3.0,pi=float(np.median(pi)),
                ps=float(np.median(ps)),ed=float(np.median(ed)),nbeat=len(pi))
def cardiac_fraction(sig,hr):
    x=np.asarray(sig,float); x=x[np.isfinite(x)]
    if x.size<20*FS_TRUE or not np.isfinite(hr) or hr<=0: return np.nan
    f,P=welch(x-x.mean(),fs=FS_TRUE,nperseg=int(20*FS_TRUE)); f0=hr/60.0
    m=np.zeros_like(f,bool)
    for k in range(1,6): m|=np.abs(f-k*f0)<0.12
    tot=P[(f>0.02)&(f<FS_TRUE/2-1)].sum()
    return float(P[m].sum()/tot) if tot>0 else np.nan
def capno(co2):
    x=np.asarray(co2,float); x=x[np.isfinite(x)]
    if x.size<10*FS_TRUE: return dict(et=np.nan,mean=np.nan,p90=np.nan,insp=np.nan,rr=np.nan)
    pk,_=find_peaks(x,distance=_MIN_BR,prominence=2.0)
    return dict(et=float(np.percentile(x,95)),mean=float(x.mean()),
                p90=float(np.percentile(x,90)),insp=float(np.percentile(x,5)),
                rr=float(len(pk)/(x.size/FS_TRUE/60.0)))
def reduce_epoch(block,idx):
    r=rpeaks(block[:,idx['ecg']])
    hr=float(60.0*FS_TRUE/np.median(np.diff(r))) if len(r)>5 else np.nan
    out=dict(hr=hr,nbeat=int(len(r)),dur_s=float(block.shape[0]/FS_TRUE))
    for ch in ('mcar','mcal','radi','abp'):
        if ch not in idx: continue
        for k,v in beat_stats(block[:,idx[ch]],r).items(): out[f'{ch}_{k}']=v
        out[f'{ch}_cardiac']=cardiac_fraction(block[:,idx[ch]],hr)
    for k,v in capno(block[:,idx['co2']]).items(): out[f'co2_{k}']=v
    return out
