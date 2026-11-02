"""Frozen-model reproduction and the authorised 112/116 sensitivity comparison.
Run with the existing project Python environment. No discovery models are fit.
"""
from pathlib import Path
import json, hashlib, sys
import numpy as np
import pandas as pd
from scipy import stats

OUT=Path(__file__).resolve().parents[1]; P=OUT.parent; ROOT=P.parent
T=OUT/'Tables'; A=OUT/'Audit'; T.mkdir(parents=True,exist_ok=True); A.mkdir(exist_ok=True)
MET=['A','B','C','D','E','A-C','A-D','D-E','A-E','B-C']
GENES=pd.read_csv(P/'stage12_frozen_gene_sets.csv').query('node in ["A","B","C","D","E"]')
MOD={k:g.gene.tolist() for k,g in GENES.groupby('node',sort=False)}
GROUPS=['DOR_young','NOR_young','DOR_old','NOR_old']
input_paths=[ROOT/f for f in ['GSE286918_count.xlsx','GSE286918_atlas_comparison.csv','GSE286918_ensembl_symbol.csv','GSE286918_sample_audit.csv']]
input_paths += [P/f for f in ['stage12_topology.py','stage12_frozen_gene_sets.csv','stage12_primary_model.csv','stage12_primary_strata.csv','stage12_sample_scores.csv','stage12_coverage.csv','GSE286918_mito_fraction.csv','stage12_external_effects.csv','stage12_external_exact_permutation.csv']]
hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in input_paths}
(A/'source_hashes_before.json').write_text(json.dumps(hashes,indent=2),encoding='utf8')
audit=[]; coverage=[]; gene_coverage=[]

def numeric_audit(ds, q, sym, cols):
    raw=q[cols]; val=raw.apply(pd.to_numeric,errors='coerce')
    invalid=raw.notna() & val.isna()
    examples=[]
    for i,j in np.argwhere(invalid.to_numpy())[:12]:
        examples.append({'row':int(i),'column':cols[j],'value':str(raw.iloc[i,j])})
    row=dict(dataset=ds,source_rows=len(q),expression_columns=len(cols),
             original_missing_cells=int(raw.isna().sum().sum()),new_coercion_missing_cells=int(invalid.sum().sum()),
             filled_zero_cells=int(val.isna().sum().sum()),negative_cells=int((val<0).sum().sum()),
             missing_symbol_rows=int(q[sym].isna().sum()),duplicate_symbol_rows=int(q[sym].duplicated(keep=False).sum()),
             STARD1_rows=int(q[sym].eq('STARD1').sum()),STAR_rows=int(q[sym].eq('STAR').sum()),invalid_examples=json.dumps(examples))
    audit.append(row)
    return val

def score(ds, mat, kind):
    mat=mat.rename(columns={'STARD1':'STAR'})
    collided=int(mat.columns.duplicated().sum())
    if collided: mat=mat.T.groupby(level=0).sum().T
    log=np.log2(mat+.5) if kind=='FPKM' else np.log2(mat.div(mat.sum(axis=1),axis=0)*1e6+.5)
    z=(log-log.mean())/log.std(ddof=1).replace(0,np.nan)
    sc=pd.DataFrame(index=mat.index)
    for node,genes in MOD.items():
        usable=[g for g in genes if g in z and z[g].notna().all()]
        sc[node]=z[usable].mean(axis=1)
        coverage.append(dict(dataset=ds,node=node,measured=len(usable),total=len(genes),genes=';'.join(usable),alias_collision_columns=collided))
        for g in genes:
            status='usable' if g in usable else ('absent' if g not in z else 'constant_or_nonfinite')
            gene_coverage.append(dict(dataset=ds,node=node,gene=g,status=status))
    for pair in MET[5:]:
        l,r=pair.split('-');sc[pair]=sc[l]-sc[r]
    for g in ['FSHR','NR5A1','CYP19A1']:
        sc['gene_'+g]=z[g]
    return sc

def fit(sc,groups):
    age=np.where(groups.str.endswith('old'),.5,-.5);res=np.where(groups.str.startswith('DOR'),.5,-.5)
    X=np.column_stack([np.ones(len(sc)),res,age,age*res]);inv=np.linalg.inv(X.T@X)
    lev=np.sum((X@inv)*X,axis=1);df=len(sc)-4;rows=[];strata=[]
    for metric in MET:
        y=sc[metric].to_numpy(); beta=inv@X.T@y;err=y-X@beta
        cov=inv@(X.T@((err/(1-lev))[:,None]**2*X))@inv
        for name,L in [('reserve',[0,1,0,0]),('age',[0,0,1,0]),('interaction',[0,0,0,1]),('young_DOR-NOR',[0,1,0,-.5]),('old_DOR-NOR',[0,1,0,.5])]:
            L=np.array(L);effect=L@beta;se=np.sqrt(L@cov@L);ci=stats.t.ppf(.975,df)*se
            row=dict(metric=metric,term=name,effect=effect,se=se,low=effect-ci,high=effect+ci,p=2*stats.t.sf(abs(effect/se),df),n=len(sc),residual_df=df)
            (rows if name in ['reserve','age','interaction'] else strata).append(row)
    m=pd.DataFrame(rows);pv=m.p.to_numpy();order=np.argsort(pv)
    q=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1].clip(0,1)
    m.loc[order,'q30']=q
    return m,pd.DataFrame(strata)

print('Reading original 116-sample matrix',flush=True)
raw=pd.read_excel(ROOT/'GSE286918_count.xlsx')
cols=[c for c in raw if c!='gene_id']
raw_num=raw[cols].apply(pd.to_numeric,errors='coerce')
assert raw_num.notna().all().all(),'Missing primary counts: record before any analysis change'
assert (raw_num>=0).all().all()
mapping=pd.read_csv(ROOT/'GSE286918_atlas_comparison.csv',usecols=['gene_id','symbol']).drop_duplicates('gene_id')
joined=raw.set_index('gene_id').join(mapping.set_index('gene_id'))
primary_mapping={'count_rows':len(raw),'count_id_duplicates':int(raw.gene_id.duplicated().sum()),'mapping_rows_after_first_id':len(mapping),
                 'unmapped_rows':int(joined.symbol.isna().sum()),'mapped_rows':int(joined.symbol.notna().sum()),
                 'duplicate_mapped_symbol_rows':int(joined.loc[joined.symbol.notna(),'symbol'].duplicated(keep=False).sum())}
q=joined.dropna(subset=['symbol']).reset_index()
numeric_audit('GSE286918',q,'symbol',cols)
v=q.set_index('symbol')[cols].groupby(level=0).sum().T
meta=pd.read_csv(ROOT/'GSE286918_sample_audit.csv');lu=dict(zip(meta.description.str.replace('Library name: ','',regex=False),meta.group))
groups=pd.Series([lu[s] for s in v.index],index=v.index)
mt=pd.read_csv(P/'GSE286918_mito_fraction.csv').set_index('sample')
mtmap=pd.read_csv(ROOT/'GSE286918_ensembl_symbol.csv').drop_duplicates('gene_id').set_index('gene_id')
counts=raw.set_index('gene_id')
mtids=mtmap.index[mtmap.symbol.str.startswith('MT-',na=False)].intersection(counts.index)
mt_recalc=counts.loc[mtids].sum()/counts.sum()
mt_error=float(np.max(np.abs(mt_recalc.loc[mt.index]-mt.mitochondrial_fraction)))
assert mt_error<1e-12
sample_qc=mt.copy();sample_qc['recalculated_fraction']=mt_recalc;sample_qc['included_112']=sample_qc.mitochondrial_fraction<=.20
sample_qc['group_from_metadata']=groups;assert sample_qc.group.eq(sample_qc.group_from_metadata).all()
sample_qc.reset_index().to_csv(T/'Table_S5a_sample_inclusion.csv',index=False)
excluded=sample_qc.loc[~sample_qc.included_112].reset_index();excluded.to_csv(A/'excluded_samples.csv',index=False)
keep=mt.index[mt.mitochondrial_fraction<=.20]
old_scores=pd.read_csv(P/'stage12_sample_scores.csv');old_main=old_scores.query('dataset=="GSE286918"').set_index('sample')
sc112=score('GSE286918',v.loc[keep],'counts');sc116=score('GSE286918_116',v,'counts')
score_error=float(np.max(np.abs(sc112.loc[old_main.index,MET].to_numpy()-old_main[MET].to_numpy())))
assert score_error<1e-10,score_error
m112,s112=fit(sc112,groups.loc[keep]);m116,s116=fit(sc116,groups)
old_m=pd.read_csv(P/'stage12_primary_model.csv').query('metric in @MET')
joined_m=m112.merge(old_m,on=['metric','term'],suffixes=('_new','_old'))
model_errors={c:float(np.max(np.abs(joined_m[c+'_new']-joined_m[c+'_old']))) for c in ['effect','se','low','high','p','q30']}
assert max(model_errors.values())<1e-10,model_errors
# Independent four-cell calculation of HC3 contrasts and uncertainties.
ind=[]
for met in MET:
    tmp=pd.DataFrame({'y':sc112[met],'group':groups.loc[keep]})
    g=tmp.groupby('group').y;means=g.mean().reindex(GROUPS).to_numpy();vv=(g.var(ddof=1)/(g.count()-1)).reindex(GROUPS).to_numpy()
    for term,w in [('reserve',[.5,-.5,.5,-.5]),('age',[-.5,-.5,.5,.5]),('interaction',[-1,1,1,-1])]:
        w=np.array(w);r=m112.query('metric==@met and term==@term').iloc[0]
        ind.append(dict(metric=met,term=term,cell_effect=w@means,cell_se=np.sqrt(w*w@vv),effect_error=abs(w@means-r.effect),se_error=abs(np.sqrt(w*w@vv)-r.se)))
ind=pd.DataFrame(ind);assert ind[['effect_error','se_error']].to_numpy().max()<1e-10
ind.to_csv(A/'independent_HC3_cell_reconstruction.csv',index=False)
for n,m,s,sc,g in [(112,m112,s112,sc112,groups.loc[keep]),(116,m116,s116,sc116,groups)]:
    m.to_csv(T/f'Table_S5_model_{n}.csv',index=False)
    s.to_csv(T/f'Table_S5_strata_{n}.csv',index=False)
    sc.assign(group=g).reset_index(names='sample').to_csv(T/f'Table_S5_scores_{n}.csv',index=False)
comparison=m112.merge(m116,on=['metric','term'],suffixes=('_112','_116'))
comparison['direction_consistent']=np.sign(comparison.effect_112)==np.sign(comparison.effect_116)
comparison['effect_change_116_minus_112']=comparison.effect_116-comparison.effect_112
comparison['effect_change_percent_of_abs112']=100*comparison.effect_change_116_minus_112/abs(comparison.effect_112)
comparison['q05_112']=comparison.q30_112<.05;comparison['q05_116']=comparison.q30_116<.05
comparison['q05_status_changed']=comparison.q05_112!=comparison.q05_116
comparison['p05_status_changed']=(comparison.p_112<.05)!=(comparison.p_116<.05)
comparison.to_csv(T/'Table_S5b_112_vs_116_all30.csv',index=False)
comparison.query('term=="reserve"').to_csv(T/'Table_S5c_112_vs_116_reserve.csv',index=False)

print('Reproduced primary model; auditing external input processing',flush=True)
external_specs=[('GSE193136',P/'stage3_sources/GSE193136_fpkm.txt.gz','gene_name','count.','counts'),
                ('GSE232306',P/'stage3_sources/GSE232306_fpkm.xlsx','gene_name','count.','counts'),
                ('GSE274832',P/'stage5_sources/GSE274832_gene_fpkm.xls.gz','gene_name',None,'FPKM'),
                ('GSE191322',P/'stage12_sources/GSE191322_gene_exp.xls.gz','tracking_id','_COUNT','counts')]
external_errors={}
for ds,path,sym,pattern,kind in external_specs:
    hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    q=pd.read_excel(path) if path.suffix=='.xlsx' else pd.read_csv(path,sep='\t')
    cols=([c for c in q if str(c).startswith(pattern)] if pattern=='count.' else
          [c for c in q if str(c).endswith(pattern)] if pattern else [f'{g}{i}' for g in ['DOR','NOR'] for i in [1,2,3]])
    num=numeric_audit(ds,q,sym,cols)
    num.index=q[sym];mat=num.fillna(0).groupby(level=0).sum().T
    mat.index=mat.index.str.replace('count.','',regex=False).str.replace('_COUNT','',regex=False)
    sc=score(ds,mat,kind);old=old_scores.query('dataset==@ds').set_index('sample')
    err=float(np.max(np.abs(sc.loc[old.index,MET].to_numpy()-old[MET].to_numpy())))
    external_errors[ds]=err;assert err<1e-10,(ds,err)

pd.DataFrame(audit).to_csv(A/'expression_input_audit.csv',index=False)
pd.DataFrame(coverage).to_csv(A/'reproduced_node_coverage.csv',index=False)
pd.DataFrame(gene_coverage).to_csv(T/'Table_S1b_gene_coverage.csv',index=False)
GENES.to_csv(T/'Table_S1a_frozen_genes.csv',index=False)
old_cov=pd.read_csv(P/'stage12_coverage.csv').query('node in ["A","B","C","D","E"]')
new_cov=pd.DataFrame(coverage).query('dataset!="GSE286918_116"')
cov_check=old_cov.merge(new_cov,on=['dataset','node'],suffixes=('_old','_new'))
assert (cov_check.measured_old==cov_check.measured_new).all()
assert (cov_check.genes_old==cov_check.genes_new).all()
new_cov.to_csv(T/'Table_S1c_node_coverage.csv',index=False)
# Descriptive E-member values from existing saved scores, without a new single-gene model.
e=old_main.reset_index().melt(id_vars=['sample','dataset','group'],value_vars=['gene_FSHR','gene_NR5A1','gene_CYP19A1'],var_name='gene',value_name='gene_z_score')
e.to_csv(T/'Table_S6a_E_gene_existing_scores.csv',index=False)
e.groupby(['gene','group']).gene_z_score.agg(['count','mean','std']).reset_index().to_csv(T/'Table_S6b_E_gene_descriptive.csv',index=False)
pd.read_csv(P/'stage12_primary_model.csv').query('metric in ["gene_FSHR","gene_CYP19A1"]').to_csv(T/'Table_S6c_existing_single_gene_HC3.csv',index=False)

summary=dict(status='PASS: no material implementation discrepancy found',primary_score_max_error=score_error,
             model_max_errors=model_errors,external_score_max_errors=external_errors,mito_fraction_max_error=mt_error,
             mitochondrial_gene_count=len(mtids),primary_mapping=primary_mapping,excluded=excluded.to_dict('records'),
             groups_112=groups.loc[keep].value_counts().to_dict(),groups_116=groups.value_counts().to_dict(),
             residual_df_112=108,residual_df_116=112,BH_tests_per_model=30,all_reserve_directions_consistent=bool(comparison.query('term=="reserve"').direction_consistent.all()),
             q05_changes=comparison.query('q05_status_changed')[['metric','term','q30_112','q30_116']].to_dict('records'),
             processing_audit=audit,versions={'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__})
(A/'stage14_audit_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf8')
for rel,h in hashes.items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h,rel
(A/'source_hashes_verified.json').write_text(json.dumps(hashes,indent=2),encoding='utf8')
print(comparison.query('term=="reserve"')[['metric','effect_112','effect_116','q30_112','q30_116','direction_consistent','q05_status_changed']].to_string(index=False),flush=True)
print('INPUT_AND_MODEL_AUDIT_PASS',flush=True)
