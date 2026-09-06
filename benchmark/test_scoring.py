"""Known-answer identity scoring checks independent of tracker output."""
import json
import numpy as np
from score_tracks import metric,ROOT
labels=np.array([0,1,0,1,0,1,0,1]);frames=np.repeat(np.arange(4),2)
perfect=[np.arange(0,8,2),np.arange(1,8,2)]
m=metric(perfect,labels,frames,[4,4]);assert m['association_idf1']==1 and m['idf1']==1
fragmented=[p[:2] for p in perfect]+[p[2:] for p in perfect]
assert metric(fragmented,labels,frames,[4,4])['association_idf1']==.5
assert metric([np.arange(3)],np.zeros(3,int),np.arange(3),[4])['association_idf1']==1
assert metric([np.arange(3)],np.zeros(3,int),np.arange(3),[4])['idf1']==6/7
assert metric([np.arange(4)],np.array([0,0,0,-1]),np.arange(4),[3])['association_idf1']==6/7
try:metric([perfect[0],perfect[0]],labels,frames,[4,4])
except AssertionError:pass
else:raise AssertionError('Duplicate assignment must fail')
assert np.isnan(metric([],[],[],[])['association_idf1'])
(ROOT/'scoring_test_results.json').write_text(json.dumps(dict(passed=True,
    checks=['perfect identities','fragmentation penalty','conditioning on fixed detections',
            'false retention penalty','duplicate detection rejected','empty-scene undefined accuracy']),indent=2))
print('PASS: 6 scoring test groups')
