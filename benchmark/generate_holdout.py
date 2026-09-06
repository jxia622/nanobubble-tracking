"""Fresh fixed-localization test scenarios; no ultrasound preprocessing.

Run only after development. The protocol and candidate hash are recorded before
the new test trajectories are evaluated. If a candidate changes afterward, this
set becomes validation data and a fresh holdout is required for a final claim.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.io import savemat

ROOT=Path(__file__).resolve().parent
CODE=ROOT.parent/'src/nbtracker.m'
SCALE=np.array([12.32,12.41446725317693]);FPS=2000;T=192
CONDITIONS=['continuous_sparse','intermittent20','intermittent40','burst_dropout',
    'crowded_dropout','curved_pulsatile','false_detections','births_deaths',
    'close_crossings','slow_stationary','high_jitter','short_trajectories',
    'null_random','persistent_false']
PRIMARY=CONDITIONS[1:11]
protocol=dict(description='Tracking only, identical localization lists for every method',
    candidate_sha256=hashlib.sha256(CODE.read_bytes()).hexdigest(),seeds=list(range(10)),
    seed_base=202609060,frames=T,fps=FPS,pixel_size_um=SCALE.tolist(),
    min_detections_all_methods=16,legacy_max_link_distance_pixels=8,
    methods=['legacy_default','legacy_gap6','nb_position_only','nb_no_evidence','nb'],
    primary_conditions=PRIMARY,control_conditions=['continuous_sparse','null_random'],
    stress_conditions=['short_trajectories','persistent_false'],
    primary_metric='detection-conditioned IDF1; global one-to-one identity mapping',
    secondary_metrics=['full IDF1','link precision','identity switches','80% recovery','false retained points','runtime'],
    success_criteria=[
        'Paired 95% bootstrap CI for mean primary association-IDF1 gain over legacy_default is entirely above zero',
        'Paired 95% bootstrap CI for mean primary association-IDF1 gain over legacy_gap6 is entirely above zero',
        'Continuous sparse mean association-IDF1 is not more than 3 percentage points below legacy_default',
        'Mean false retained observations in null_random does not exceed legacy_default by more than 1 per 192-frame sequence',
        'All invariants and MATLAB integration tests pass; fixed localizations remain unchanged'],
    limitations='NB dropout and motion stress tests, not calibrated acoustic or in-vivo validation',
    development='Three existing seeds per condition; no holdout scores used to choose this frozen candidate')
if (ROOT/'holdout_protocol.json').exists():
    old=json.loads((ROOT/'holdout_protocol.json').read_text())
    if old['candidate_sha256']!=protocol['candidate_sha256']:
        raise SystemExit('Candidate changed after holdout creation: preserve this set and create a fresh evaluation round.')
(ROOT/'holdout_protocol.json').write_text(json.dumps(protocol,indent=2))
dest=ROOT/'holdout_inputs';dest.mkdir(exist_ok=True)
manifest=[]
def cells(items):
    a=np.empty((len(items),1),object)
    for i,v in enumerate(items):a[i,0]=v
    return a

for c,condition in enumerate(CONDITIONS):
    for seed in range(10):
        rng=np.random.default_rng(protocol['seed_base']+c*100+seed)
        n=96 if condition=='crowded_dropout' else (12 if condition=='close_crossings' else 24)
        if condition=='continuous_sparse':n=16
        if condition=='null_random':n=0
        noise=2 if condition=='high_jitter' else (.5 if condition=='close_crossings' else 1)
        p=.6 if condition=='intermittent40' else (.85 if condition=='close_crossings' else .8)
        if condition=='continuous_sparse':p=1
        initial=rng.uniform(20,148,(n,2))
        angles=rng.uniform(0,2*np.pi,n)
        speeds=rng.uniform(.5,5,n)*1000
        if condition=='continuous_sparse':speeds[:]=1000
        if condition=='slow_stationary':speeds[:]=20
        t=np.arange(T)/FPS
        velocities=np.broadcast_to(np.column_stack((np.sin(angles),np.cos(angles)))*speeds[:,None],(T,n,2)).copy()
        if condition=='curved_pulsatile':
            angle=angles[None]+t[:,None]*rng.uniform(-40,40,n)[None]
            speeds=6000*(1+.4*np.sin(2*np.pi*6*t[:,None]+rng.uniform(0,2*np.pi,n)[None]))
            velocities=np.stack((np.sin(angle)*speeds,np.cos(angle)*speeds),axis=-1)
        positions=np.empty((T,n,2));positions[0]=initial
        if n:
            brown=rng.normal(0,np.sqrt(2*1.816/FPS),(T-1,n,2))
            positions[1:]=initial+np.cumsum((velocities[:-1]/FPS+brown)/SCALE,axis=0)
        if condition=='close_crossings':
            for k in range(0,n,2):
                z=25+(k//2)*22
                positions[:,k,0]=z;positions[:,k+1,0]=z+.3
                step=rng.uniform(.25,.45)
                positions[:,k,1]=80+step*(np.arange(T)-T/2)
                positions[:,k+1,1]=80-step*(np.arange(T)-T/2)
                velocities[:,k]=[0,step*SCALE[1]*FPS]
                velocities[:,k+1]=[0,-step*SCALE[1]*FPS]
        exists=np.ones((T,n),bool)
        if condition in ['births_deaths','short_trajectories']:
            exists[:]=False
            for k in range(n):
                lifetime=int(rng.integers(48,T+1) if condition=='births_deaths' else rng.integers(8,33))
                start=int(rng.integers(0,T-lifetime+1));exists[start:start+lifetime,k]=True
        detected=(rng.random((T,n))<p)&exists
        if condition=='burst_dropout':
            detected=exists.copy()
            for k in range(n):
                f=0
                while f<T:
                    if rng.random()<.08:
                        gap=int(rng.integers(1,6));detected[f:f+gap,k]=False;f+=gap
                    else:f+=1
        measured=positions+rng.normal(0,noise,positions.shape)
        lists=[];labels=[]
        persistent=rng.uniform(20,148,(3,2))
        for f in range(T):
            pts=measured[f,detected[f]];lab=np.flatnonzero(detected[f])
            false_count=int(rng.poisson(20 if condition=='false_detections' else (40 if condition=='null_random' else 0)))
            if false_count:
                pts=np.vstack((pts,rng.uniform(20,148,(false_count,2))))
                lab=np.r_[lab,np.full(false_count,-1)]
            if condition=='persistent_false':
                pts=np.vstack((pts,persistent+rng.normal(0,.3,persistent.shape)))
                lab=np.r_[lab,[-1,-1,-1]]
            order=rng.permutation(len(pts));lists.append(pts[order]);labels.append(lab[order])
        key=f'{condition}_seed{seed:02d}'
        savemat(dest/f'{key}.mat',dict(SR_Localizations=cells(lists),labels=cells(labels),
            truth=positions,exists=exists,velocities=velocities,frame_times=t,scale_um=SCALE),do_compression=True)
        manifest.append(dict(key=key,condition=condition,seed=seed,frames=T,n=n,sigma_px=noise))
(ROOT/'holdout_manifest.json').write_text(json.dumps(manifest,indent=2))
print(f'Prepared {len(manifest)} fresh fixed-localization cases; candidate frozen as {protocol["candidate_sha256"]}')
