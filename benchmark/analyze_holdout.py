"""Audit and score the frozen candidate; no parameter selection here."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from score_tracks import metric,tracks_from_mat,write_csv

ROOT=Path(__file__).resolve().parent
protocol=json.loads((ROOT/'holdout_protocol.json').read_text())
candidate=ROOT.parent/'src/nbtracker.m'
assert hashlib.sha256(candidate.read_bytes()).hexdigest()==protocol['candidate_sha256'], 'Frozen candidate changed'
jobs=json.loads((ROOT/'holdout_manifest.json').read_text());rows=[]
methods=protocol['methods'];audited=0
for job in jobs:
    inp=loadmat(ROOT/'holdout_inputs'/f"{job['key']}.mat")
    out=loadmat(ROOT/'holdout_results'/f"{job['key']}.mat",struct_as_record=False)['result'][0,0]
    points=inp['SR_Localizations'].ravel();flat=np.concatenate(points)
    labels=np.concatenate([a.ravel().astype(int) for a in inp['labels'].ravel()])
    frames=np.repeat(np.arange(len(points)),[len(a) for a in points])
    gcount=inp['exists'].sum(axis=0).astype(int)
    assert len(labels)==len(flat)
    for name in methods:
        result=getattr(out,name)[0,0];tracks=tracks_from_mat(result.tracks)
        assert all(len(tr)>=16 for tr in tracks)
        row=dict(condition=job['condition'],seed=job['seed'],method=name,
                 **metric(tracks,labels,frames,gcount),seconds=float(result.seconds.item()))
        rows.append(row)
        if name.startswith('nb'):
            records=result.details[0,0].records.ravel()
            assert len(records)==len(tracks)
            for tr,record in zip(tracks,records):
                assert np.array_equal(record.positionsPixels,flat[tr]),'Localization coordinates changed'
                assert np.array_equal(record.frameIndices.ravel().astype(int)-1,frames[tr])
                assert np.allclose(record.timesSeconds.ravel(),inp['frame_times'].ravel()[frames[tr]],rtol=0,atol=0)
                assert np.isfinite(record.smoothedPositionsPixels).all()
                assert np.isfinite(record.velocityUmPerSec).all()
        audited+=1
write_csv(ROOT/'holdout_metrics.csv',rows)
conditions=list(dict.fromkeys(j['condition'] for j in jobs))
summary={}
for c in conditions:
    summary[c]={}
    for m in methods:
        data=[r for r in rows if r['condition']==c and r['method']==m]
        summary[c][m]={k:float(np.mean([r[k] for r in data])) for k in data[0] if k not in ['condition','seed','method']}

lookup={(r['condition'],r['seed'],r['method']):r for r in rows}
primary=protocol['primary_conditions'];rng=np.random.default_rng(45203)
gain={}
for comparator in ['legacy_default','legacy_gap6','nb_position_only','nb_no_evidence']:
    diffs=np.array([[lookup[c,s,'nb']['association_idf1']-lookup[c,s,comparator]['association_idf1']
                     for s in protocol['seeds']] for c in primary])
    # Paired within-condition resampling; conditions receive equal fixed weight.
    boot=np.zeros(10000)
    for d in diffs:boot+=d[rng.integers(0,len(d),(10000,len(d)))].mean(axis=1)/len(diffs)
    lo,hi=np.quantile(boot,[.025,.975])
    gain[comparator]=dict(mean_pp=float(diffs.mean()*100),ci95_pp=[float(lo*100),float(hi*100)])
mean_primary={m:float(np.mean([summary[c][m]['association_idf1'] for c in primary])) for m in methods}
checks={
    'positive_gain_vs_original':gain['legacy_default']['ci95_pp'][0]>0,
    'positive_gain_vs_longer_gap':gain['legacy_gap6']['ci95_pp'][0]>0,
    'continuous_control':summary['continuous_sparse']['nb']['association_idf1']>=summary['continuous_sparse']['legacy_default']['association_idf1']-.03,
    'null_control':summary['null_random']['nb']['false_retained_points']<=summary['null_random']['legacy_default']['false_retained_points']+1,
    'behavioral_tests':json.loads((ROOT.parent/'results'/'unit_test_results.json').read_text())['passed'],
    'assignment_and_immutable_coordinate_audits':audited==len(jobs)*len(methods)
}
report=dict(candidate_sha256=protocol['candidate_sha256'],cases=len(jobs),audited_outputs=audited,
            primary_mean=mean_primary,paired_gain=gain,checks=checks,
            benchmark_success=all(checks.values()),conditions=summary,
            ci_method='10,000 paired bootstrap resamples within each fixed condition; equal condition weights')
# JSON null rather than non-standard NaN for undefined null-scene accuracy.
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):return [clean(v) for v in x]
    if isinstance(x,float) and not np.isfinite(x):return None
    return x
(ROOT/'holdout_summary.json').write_text(json.dumps(clean(report),indent=2))
print('PRIMARY ASSOCIATION IDF1',mean_primary)
print('PAIRED GAIN (percentage points)',json.dumps(gain,indent=2))
print('CHECKS',checks)
for c in conditions:
    print(c,' | '.join(f"{m}: {summary[c][m]['association_idf1']*100:.1f}%" for m in ['legacy_default','legacy_gap6','nb']))
