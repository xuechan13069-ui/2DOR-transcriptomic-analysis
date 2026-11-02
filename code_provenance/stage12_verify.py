"""Independent reconstruction of the 2x2 mean contrasts and exact small-n tests."""
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
P=Path(__file__).resolve().parent
s=pd.read_csv(P/'stage12_sample_scores.csv'); m=pd.read_csv(P/'stage12_primary_model.csv')
z=s[s.dataset.eq('GSE286918')];out=[]
for metric in ['A','B','C','D','E','A-C','A-D','D-E','A-E','B-C']:
    means=z.groupby('group')[metric].mean()
    young=means.DOR_young-means.NOR_young;old=means.DOR_old-means.NOR_old
    age_nor=means.NOR_old-means.NOR_young;age_dor=means.DOR_old-means.DOR_young
    for term,estimate in [('reserve',(young+old)/2),('age',(age_nor+age_dor)/2),('interaction',old-young)]:
        fit=m[(m.metric==metric)&(m.term==term)].effect.iloc[0]
        out.append(dict(metric=metric,term=term,independent=estimate,model=fit,abs_error=abs(estimate-fit)))
pd.DataFrame(out).to_csv(P/'stage12_verification_means.csv',index=False)
rows=[]
for ds,g1,g0 in [('GSE191322','COS','mild'),('GSE193136','DOR','NOR'),('GSE232306','DOR','NOR'),('GSE274832','DOR','NOR')]:
    a=s[s.dataset.eq(ds)]
    for metric in ['D-E','A-E','A-C','B-C']:
        x=a.loc[a.group.eq(g1),metric].to_numpy(); y=a.loc[a.group.eq(g0),metric].to_numpy(); pool=np.r_[x,y];n=len(x)
        obs=x.mean()-y.mean(); null=[]
        for comb in itertools.combinations(range(len(pool)),n):
            ix=np.array(comb); mask=np.ones(len(pool),bool);mask[ix]=False
            null.append(pool[ix].mean()-pool[mask].mean())
        p=np.mean(np.abs(null)>=abs(obs)-1e-12)
        rows.append(dict(dataset=ds,metric=metric,contrast=g1+'-'+g0,effect=obs,permutation_p=p,permutations=len(null)))
pd.DataFrame(rows).to_csv(P/'stage12_external_exact_permutation.csv',index=False)
print('Maximum independent 2x2 reconstruction error:',max(x['abs_error'] for x in out))
