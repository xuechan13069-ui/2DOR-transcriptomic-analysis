"""Conditional NOR aging and reserve-versus-age contrasts, using frozen scores."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
P=Path(__file__).resolve().parent
s=pd.read_csv(P/'stage12_sample_scores.csv');s=s[s.dataset.eq('GSE286918')]
a=np.where(s.group.str.endswith('old'),.5,-.5);r=np.where(s.group.str.startswith('DOR'),.5,-.5)
X=np.column_stack([np.ones(len(a)),r,a,r*a]); inv=np.linalg.inv(X.T@X);rows=[]
for metric in ['A','B','C','D','E','A-C','A-D','D-E','A-E','B-C']:
 y=s[metric].to_numpy(); b=inv@X.T@y; e=y-X@b; h=np.sum((X@inv)*X,axis=1); cov=inv@X.T@np.diag((e/(1-h))**2)@X@inv
 for name,l in [('NOR_old_minus_young',[0,0,1,-.5]),('DOR_old_minus_young',[0,0,1,.5]),('reserve_minus_age',[0,1,-1,0])]:
  l=np.array(l);est=l@b;se=np.sqrt(l@cov@l);q=stats.t.ppf(.975,len(y)-4)*se
  rows.append(dict(metric=metric,contrast=name,effect=est,low=est-q,high=est+q,p=2*stats.t.sf(abs(est/se),len(y)-4)))
pd.DataFrame(rows).to_csv(P/'stage12_age_conditional_contrasts.csv',index=False)
