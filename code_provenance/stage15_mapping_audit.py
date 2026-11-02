"""Bounded annotation/count-retention audit. No alternative disease model is fit."""
from pathlib import Path
import json,hashlib,re
import numpy as np
import pandas as pd
import requests

O=Path(__file__).resolve().parents[1];P=O.parent;R=P.parent
A=O/'Audit';T=O/'Tables';S=O/'Sources'
for d in [A,T,S]:d.mkdir(parents=True,exist_ok=True)
inputs=[R/n for n in ['GSE286918_count.xlsx','GSE286918_atlas_comparison.csv','GSE286918_ensembl_symbol.csv','GSE286918_four_cell_exploratory.csv','GSE286918_sample_audit.csv','audit_gse286918.py','compare_aging_atlases.py','map_gse286918_genes.py']]
hashes={str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs}
(A/'input_hashes.json').write_text(json.dumps(hashes,indent=2),encoding='utf8')
url='https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
ref=S/'hgnc_complete_set_2026-10-31.txt'
if not ref.exists():
    r=requests.get(url,timeout=45);r.raise_for_status();ref.write_bytes(r.content)
    (S/'hgnc_download_metadata.json').write_text(json.dumps({'url':url,'accessed':'2026-10-31','http_Last_Modified':r.headers.get('Last-Modified'),'sha256':hashlib.sha256(r.content).hexdigest(),'purpose':'audit only; never substituted into primary scoring'},indent=2),encoding='utf8')
h=pd.read_csv(ref,sep='\t',dtype=str).fillna('')
print('HGNC audit reference available; reading original matrix',flush=True)
raw=pd.read_excel(R/'GSE286918_count.xlsx');ids=raw.gene_id.astype(str);cols=raw.columns[1:]
Y=raw.iloc[:,1:].to_numpy(dtype=float);assert np.isfinite(Y).all() and (Y>=0).all()
tot=Y.sum(axis=0);cpm=Y/tot[None,:]*1e6;expressed=(cpm>=1).sum(axis=1)>=10
exp=pd.read_csv(R/'GSE286918_four_cell_exploratory.csv',usecols=['gene_id'])
assert set(ids[expressed])==set(exp.gene_id),'Historical filter reconstruction conflict'
mapping=pd.read_csv(R/'GSE286918_atlas_comparison.csv',usecols=['gene_id','symbol'])
second=pd.read_csv(R/'GSE286918_ensembl_symbol.csv')
assert set(mapping.gene_id)==set(exp.gene_id)
assert not mapping.gene_id.duplicated().any()
m=mapping.set_index('gene_id').symbol;symbol=ids.map(m);retained=symbol.notna().to_numpy()
status=np.where(~expressed,'below_historical_expression_filter',np.where(retained,'retained','queried_but_no_persisted_symbol'))
rows=pd.DataFrame({'gene_id':ids,'historical_filter_pass':expressed,'historical_symbol':symbol,'retained':retained,'exclusion_category':status,
                   'samples_CPM_ge1':(cpm>=1).sum(axis=1),'total_counts_116':Y.sum(axis=1),'mean_CPM':cpm.mean(axis=1),'max_CPM':cpm.max(axis=1)})

# Current HGNC identifies some rows; absent current annotations stay unknown.
records={}
for r in h.to_dict('records'):
    for geneid in r['ensembl_gene_id'].split('|'):
        if geneid:records.setdefault(geneid,[]).append(r)
for key in ['symbol','locus_group','locus_type','hgnc_id']:
    rows['current_HGNC_'+key]=ids.map(lambda g:'|'.join(sorted({r[key] for r in records.get(g,[])})))
rows['current_HGNC_hits']=ids.map(lambda g:len(records.get(g,[])))
rows['current_type_known']=rows.current_HGNC_locus_group.ne('')
def agreement(row):
    if pd.isna(row.historical_symbol):return 'no historical symbol'
    rs=records.get(row.gene_id,[])
    if not rs:return 'not resolved by current HGNC ID field'
    if any(row.historical_symbol==r['symbol'] for r in rs):return 'exact symbol agreement'
    if any(row.historical_symbol in (r['alias_symbol']+'|'+r['prev_symbol']).split('|') for r in rs):return 'documented alias or previous symbol'
    return 'current-reference disagreement; historical version unresolved'
rows['symbol_comparison']=rows.apply(agreement,axis=1)
rows.to_csv(T/'Mapping_S1_row_disposition.csv',index=False)
rows.groupby(['exclusion_category','current_HGNC_locus_group'],dropna=False).agg(rows=('gene_id','size'),pooled_counts=('total_counts_116','sum')).reset_index().to_csv(T/'Mapping_S2_current_type_summary.csv',index=False)
rows.query('not retained').sort_values('total_counts_116',ascending=False).head(50).to_csv(T/'Mapping_S3_largest_excluded_rows.csv',index=False)
pc=rows[(~rows.retained)&rows.current_HGNC_locus_group.eq('protein-coding gene')]
pc.to_csv(T/'Mapping_S4_excluded_current_protein_coding.csv',index=False)
rows[rows.symbol_comparison.str.startswith('current-reference disagreement')].to_csv(T/'Mapping_S5_annotation_disagreements.csv',index=False)

meta=pd.read_csv(R/'GSE286918_sample_audit.csv');lu=dict(zip(meta.description.str.replace('Library name: ','',regex=False),meta.group))
qc=pd.read_csv(P/'GSE286918_mito_fraction.csv').set_index('sample')
kept=Y[retained].sum(axis=0);unqueried=Y[~expressed].sum(axis=0);unresolved=Y[expressed&~retained].sum(axis=0)
sample=pd.DataFrame({'sample':cols,'group':[lu[c] for c in cols],'raw_total_counts':tot,'mapped_retained_counts':kept,
                     'excluded_counts':tot-kept,'count_retention_fraction':kept/tot,'below_filter_counts':unqueried,
                     'queried_unresolved_counts':unresolved,'normalization_multiplier_full_over_retained':tot/kept,
                     'log2_full_over_retained':np.log2(tot/kept),'mitochondrial_fraction':qc.loc[cols,'mitochondrial_fraction'].to_numpy()})
sample['included_112']=sample.mitochondrial_fraction<=.2
sample['reserve']=np.where(sample.group.str.startswith('DOR'),'DOR','NOR')
assert np.allclose(sample.excluded_counts,sample.below_filter_counts+sample.queried_unresolved_counts)
sample.to_csv(T/'Mapping_S6_sample_count_retention.csv',index=False)
summ=[]
for setting,x in [('all116',sample),('primary112',sample[sample.included_112])]:
    for groupcol in ['group','reserve']:
        for g,ss in x.groupby(groupcol):
            z=ss.count_retention_fraction
            summ.append({'set':setting,'grouping':groupcol,'group':g,'n':len(z),'mean':z.mean(),'sd':z.std(ddof=1),'min':z.min(),'q25':z.quantile(.25),'median':z.median(),'q75':z.quantile(.75),'max':z.max(),'pooled_count_retention':ss.mapped_retained_counts.sum()/ss.raw_total_counts.sum(),'mean_log2_full_over_retained':ss.log2_full_over_retained.mean()})
summary=pd.DataFrame(summ);summary.to_csv(T/'Mapping_S7_retention_descriptive.csv',index=False)
descdiff=[]
for setting,x in [('all116',sample),('primary112',sample[sample.included_112])]:
    for age in ['young','old']:
        y=x[x.group.eq('DOR_'+age)];z=x[x.group.eq('NOR_'+age)]
        descdiff.append({'set':setting,'age':age,'DOR_minus_NOR_mean_retention_percentage_points':100*(y.count_retention_fraction.mean()-z.count_retention_fraction.mean()),'DOR_minus_NOR_mean_log2_multiplier':y.log2_full_over_retained.mean()-z.log2_full_over_retained.mean()})
pd.DataFrame(descdiff).to_csv(T/'Mapping_S8_descriptive_group_differences.csv',index=False)

# Exact aggregation and frozen programme coverage reconstruction.
mapped=raw.set_index('gene_id').join(mapping.set_index('gene_id')).dropna(subset=['symbol'])
mat=mapped.groupby('symbol').sum().T
check_total=np.max(abs(mat.sum(axis=1).loc[cols].to_numpy()-kept));assert check_total==0
dups=rows[rows.historical_symbol.notna() & rows.historical_symbol.duplicated(keep=False)]
dups.to_csv(T/'Mapping_S9_duplicate_symbol_rows.csv',index=False)
alternate=mapping.merge(second,on='gene_id',suffixes=('_atlas','_other'))
different=alternate[alternate.symbol_atlas.fillna('NA').ne(alternate.symbol_other.fillna('NA'))]
different.to_csv(T/'Mapping_S10_two_historical_maps.csv',index=False)
panels=pd.read_csv(P/'stage12_frozen_gene_sets.csv').query('node in ["A","B","C","D","E"]')
keep=qc.index[qc.mitochondrial_fraction<=.2];v=mat.loc[keep]
log=np.log2(v.div(v.sum(axis=1),axis=0)*1e6+.5);z=(log-log.mean())/log.std(ddof=1).replace(0,np.nan)
panelrows=[];scores=pd.DataFrame(index=v.index)
for node,gg in panels.groupby('node',sort=False):
    usable=[]
    for gene in gg.gene:
        historical=rows[rows.historical_symbol.eq(gene)]
        refrows=h[h.symbol.eq(gene)];refids=set('|'.join(refrows.ensembl_gene_id).split('|'))-{''}
        rawref=rows[rows.gene_id.isin(refids)]
        available=gene in z and z[gene].notna().all()
        if available:usable.append(gene)
        panelrows.append({'node':node,'gene':gene,'usable112':available,'historical_IDs':';'.join(historical.gene_id),
                          'current_HGNC_Ensembl_IDs':';'.join(sorted(refids)),'current_IDs_present_in_raw':';'.join(rawref.gene_id),
                          'current_IDs_filter_status':';'.join(rawref.exclusion_category),
                          'current_IDs_samples_CPM_ge1':';'.join(rawref.samples_CPM_ge1.astype(str)),
                          'current_IDs_max_CPM':';'.join(rawref.max_CPM.map(lambda x:f'{x:.8g}')),
                          'historical_mapping_agreement':';'.join(historical.symbol_comparison),
                          'log_sample_SD112':float(log[gene].std(ddof=1)) if gene in log else np.nan})
    scores[node]=z[usable].mean(axis=1)
panel=pd.DataFrame(panelrows);panel.to_csv(T/'Mapping_S11_frozen_panel_identity_audit.csv',index=False)
old=pd.read_csv(P/'stage12_sample_scores.csv').query('dataset=="GSE286918"').set_index('sample')
error=float(np.max(abs(scores.loc[old.index,list('ABCDE')].to_numpy()-old[list('ABCDE')].to_numpy())))
assert error<1e-10

result={'raw_rows':len(raw),'human_stable_ID_rows':int(ids.str.fullmatch(r'ENSG\d{11}').sum()),'version_suffix_rows':int(ids.str.contains(r'\.\d+$').sum()),
        'historical_filter_rows':int(expressed.sum()),'below_filter_rows':int((~expressed).sum()),'queried_no_symbol_rows':int((expressed&~retained).sum()),'retained_rows':int(retained.sum()),
        'retained_unique_symbols':mat.shape[1],'duplicate_symbol_row_count':len(dups),'duplicate_symbol_count':dups.historical_symbol.nunique(),
        'aggregation_count_error':float(check_total),'frozen_programme_max_score_error':error,
        'historical_mapping_differences':len(different),'HGNC_rows':len(h),'HGNC_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),
        'current_reference_disagreement_retained_rows':int(rows.symbol_comparison.str.startswith('current-reference disagreement').sum()),
        'excluded_current_protein_coding_rows':len(pc),'excluded_current_protein_coding_after_filter':int(pc.historical_filter_pass.sum()),
        'retention_all116':sample.count_retention_fraction.describe().to_dict(),'retention_primary112':sample.loc[sample.included_112,'count_retention_fraction'].describe().to_dict(),
        'descriptive_group_differences':descdiff,'annotation_version':'historical MyGene backend snapshot and response multiplicity not recorded; current HGNC is audit-only',
        'gate_A':'PENDING_EVIDENCE_REVIEW; no alternative models run'}
(A/'mapping_audit_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
for rel,sha in hashes.items():assert hashlib.sha256((R/rel).read_bytes()).hexdigest()==sha
print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
