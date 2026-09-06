"""Identity scoring on fixed detections, shared by all tracking methods.

No SVD/localizer tuning is performed here. Point-list inputs carry source IDs;
saved image localizations use the pilot's frozen 30 um one-to-one matching.
Missing ground-truth observations count in full IDF1 and are excluded in the
primary detection-conditioned IDF1, which isolates association more strictly.
"""
from pathlib import Path
import json,csv
import numpy as np
from scipy.io import loadmat
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'Nanobubble tracking assessment'

def tracks_from_mat(v):return [np.asarray(x).ravel().astype(int)-1 for x in v.ravel()]

def metric(tracks,labels,frames,gt_count):
    # gt_count is observations possible for each true trajectory, including
    # variable birth/death times. Predicted nodes are never synthetic fills.
    labels=np.asarray(labels,int);frames=np.asarray(frames,int);gt_count=np.asarray(gt_count,int)
    used=np.concatenate(tracks).astype(int) if tracks else np.empty(0,int)
    assert len(np.unique(used))==len(used), 'A detection was assigned more than once'
    assert np.all((used>=0)&(used<len(labels))), 'Invalid detection index'
    n=len(gt_count); contingency=np.zeros((n,len(tracks)),int)
    retained=0; matched=0; links=0; correct=0; switches=0
    for k,tr in enumerate(tracks):
        tr=np.asarray(tr,int); assert np.all(np.diff(frames[tr])>0)
        ids=labels[tr]; retained+=len(tr);matched+=np.sum(ids>=0)
        for ident in np.unique(ids[ids>=0]):contingency[ident,k]=np.sum(ids==ident)
        links+=max(0,len(tr)-1)
        correct+=np.sum((ids[:-1]==ids[1:])&(ids[:-1]>=0))
        switches+=np.sum((ids[:-1]!=ids[1:])&(ids[:-1]>=0)&(ids[1:]>=0))
    idtp=0; recovered=0
    if n and len(tracks):
        a,b=linear_sum_assignment(-contingency)
        idtp=int(contingency[a,b].sum())
        recovered=sum(contingency[i,k]>=.8*gt_count[i] and contingency[i,k]/len(tracks[k])>=.9 for i,k in zip(a,b))
    observed_true=np.sum(labels>=0);gt=int(gt_count.sum())
    div=lambda a,b:float(a/b) if b else np.nan
    return dict(idf1=div(2*idtp,gt+retained),association_idf1=div(2*idtp,observed_true+retained),
        retained_recall=div(matched,gt),retained_precision=div(matched,retained),
        link_precision=div(correct,links),identity_switches=int(switches),
        recovery80=div(recovered,n),retained_tracks=len(tracks),retained_points=retained,
        true_points=gt,detected_true_points=int(observed_true),false_retained_points=int(retained-matched))

def old_data(key):
    inp=loadmat(OLD/'inputs'/f'{key}.mat')
    saved=loadmat(OLD/'results'/f'{key}.mat',struct_as_record=False)
    points=list(saved['SR_Localizations'].ravel())
    frames=np.repeat(np.arange(len(points)),[len(x) for x in points])
    truth=inp['truth']; T,n,_=truth.shape
    if 'labels' in inp:labels=np.concatenate([x.ravel().astype(int) for x in inp['labels'].ravel()])
    else:
        scale=inp['scale_um'].ravel();ls=[]
        for t,p in enumerate(points):
            ids=np.full(len(p),-1,int)
            if n and len(p):
                d=cdist(p*scale,truth[t]*scale)
                a,b=linear_sum_assignment(np.where(d<=30,d*d,1e12));good=d[a,b]<=30
                ids[a[good]]=b[good]
            ls.append(ids)
        labels=np.concatenate(ls)
    legacy=tracks_from_mat(saved['output'][0,0].gap1[0,0].tracks)
    return labels,frames,np.full(n,T),legacy

def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__':
    rows=[]
    for p in sorted((ROOT/'development').glob('*.mat')):
        if p.name.startswith('._'):continue
        result=loadmat(p,struct_as_record=False)['result'][0,0]
        labels,frames,gcount,legacy=old_data(p.stem)
        for name in ['legacy_default','legacy_gap6','position_only','nb']:
            tr=legacy if name=='legacy_default' else tracks_from_mat(getattr(result,name)[0,0].tracks)
            row=dict(condition=p.stem.rsplit('_seed',1)[0],seed=int(p.stem[-2:]),method=name,**metric(tr,labels,frames,gcount))
            rows.append(row)
    if not rows:raise SystemExit('No completed development outputs')
    write_csv(ROOT/'development_metrics.csv',rows)
    for condition in sorted(set(r['condition'] for r in rows)):
        print(condition, ' | '.join(f'{method}: '+f"{np.mean([r['idf1'] for r in rows if r['condition']==condition and r['method']==method])*100:.1f}%" for method in ['legacy_default','legacy_gap6','position_only','nb']))
