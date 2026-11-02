"""Frozen stage12 programme projection; see stage12_frozen_protocol."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

P=Path(__file__).resolve().parent; ROOT=P.parent
t=pd.read_csv(P/'stage10_frozen_three_axis_modules.csv')
modules={'A':list(t.loc[t.module.eq('sterol_22_core'),'gene']),
         'B':['LDLR','SCARB1','SCAP','INSIG1','INSIG2'],
         'C':['NPC1','NPC2','STAR'],
         'D':['CYP11A1','HSD3B1','HSD3B2'],
         'E':['FSHR','NR5A1','CYP19A1']}
for m in ['energy_glycolysis','energy_TCA','energy_OXPHOS']:
    modules[m]=list(t.loc[t.module.eq(m),'gene'])
pd.DataFrame([{'node':k,'gene':g,'source':'stage10 frozen' if k in ['A','E'] or k.startswith('energy') else 'stage12 frozen protocol'} for k,gs in modules.items() for g in gs]).to_csv(P/'stage12_frozen_gene_sets.csv',index=False)

def ext(path,excel=False):
    q=pd.read_excel(path) if excel else pd.read_csv(path,sep='\t')
    cols=[x for x in q if str(x).startswith('count.')]
    m=q.set_index('gene_name')[cols].apply(pd.to_numeric,errors='coerce').fillna(0).groupby(level=0).sum().T
    m.index=m.index.str.replace('count.','',regex=False); return m

def main():
    mapping=pd.read_csv(ROOT/'GSE286918_atlas_comparison.csv',usecols=['gene_id','symbol']).drop_duplicates('gene_id')
    v=pd.read_excel(ROOT/'GSE286918_count.xlsx').set_index('gene_id').join(mapping.set_index('gene_id')).dropna(subset=['symbol']).groupby('symbol').sum().T
    mt=pd.read_csv(P/'GSE286918_mito_fraction.csv').set_index('sample'); v=v.loc[mt.index[mt.mitochondrial_fraction<=.20]]
    audit=pd.read_csv(ROOT/'GSE286918_sample_audit.csv'); lu=dict(zip(audit.description.str.replace('Library name: ','',regex=False),audit.group))
    group=pd.Series([lu[x] for x in v.index],index=v.index)
    a=ext(P/'stage3_sources/GSE193136_fpkm.txt.gz'); b=ext(P/'stage3_sources/GSE232306_fpkm.xlsx',True)
    x=pd.read_csv(P/'stage5_sources/GSE274832_gene_fpkm.xls.gz',sep='\t')
    c=x.set_index('gene_name')[[f'{g}{i}' for g in ['DOR','NOR'] for i in [1,2,3]]].apply(pd.to_numeric,errors='coerce').fillna(0).groupby(level=0).sum().T
    h=pd.read_csv(P/'stage12_sources/GSE191322_gene_exp.xls.gz',sep='\t')
    cols=[j for j in h if j.endswith('_COUNT')]; d=h.set_index('tracking_id')[cols].apply(pd.to_numeric,errors='coerce').fillna(0).groupby(level=0).sum().T
    d.index=d.index.str.replace('_COUNT','',regex=False)
    def dg(m): return pd.Series(np.where(m.index.str.startswith('DOR'),'DOR','NOR'),index=m.index)
    dg191=pd.Series(np.where(d.index.str.startswith('granulosa_cell_2'),'COS','mild'),index=d.index)
    items=[('GSE286918',v,group,'counts'),('GSE193136',a,dg(a),'fpkm_relative'),('GSE232306',b,dg(b),'fpkm_relative'),('GSE274832',c,dg(c),'fpkm'),('GSE191322',d,dg191,'counts')]
    allscore=[]; cover=[]; allgene=[]
    for ds,mat,groups,kind in items:
        mat=mat.rename(columns={'STARD1':'STAR'})
        if mat.columns.duplicated().any(): mat=mat.T.groupby(level=0).sum().T
        log=np.log2(mat+0.5) if kind=='fpkm' else np.log2(mat.div(mat.sum(axis=1),axis=0)*1e6+0.5)
        z=(log-log.mean(axis=0))/log.std(axis=0,ddof=1).replace(0,np.nan)
        out=pd.DataFrame(index=mat.index)
        for node,genes in modules.items():
            usable=[g for g in genes if g in z and z[g].notna().all()]
            out[node]=z[usable].mean(axis=1) if usable else np.nan
            cover.append(dict(dataset=ds,node=node,measured=len(usable),total=len(genes),genes=';'.join(usable)))
        out['energy_equal3']=out[['energy_glycolysis','energy_TCA','energy_OXPHOS']].mean(axis=1)
        for g in ['FSHR','CYP19A1','NR5A1','CYP11A1','HSD3B1','HSD3B2','STAR']:
            out['gene_'+g]=z[g] if g in z else np.nan
        for pair in ['A-C','A-D','D-E','A-E','B-C']:
            k,l=pair.split('-');out[pair]=out[k]-out[l]
        out.insert(0,'group',groups); out.insert(0,'dataset',ds); out=out.reset_index(names='sample')
        for ix,row in out.iterrows():
            if ds=='GSE286918': ref='NOR_'+row['group'].split('_')[-1]
            elif ds=='GSE191322': ref='mild'
            else: ref='NOR'
            base=out.loc[out.group.eq(ref)]
            if row['group']==ref: base=base[base['sample']!=row['sample']]
            for node in ['A','B','C','D','E']:
                sd=base[node].std(ddof=1)
                out.loc[ix,'relative_'+node]=(row[node]-base[node].mean())/sd if sd>0 else np.nan
        for pair in ['A-C','A-D','D-E','A-E','B-C']:
            k,l=pair.split('-');out['relative_'+pair]=out['relative_'+k]-out['relative_'+l]
        allscore.append(out)
    s=pd.concat(allscore,ignore_index=True); s.to_csv(P/'stage12_sample_scores.csv',index=False)
    pd.DataFrame(cover).to_csv(P/'stage12_coverage.csv',index=False)
    metrics=['A','B','C','D','E','A-C','A-D','D-E','A-E','B-C']
    primary=s[s.dataset.eq('GSE286918')].copy()
    age=np.where(primary.group.str.endswith('old'),.5,-.5); reserve=np.where(primary.group.str.startswith('DOR'),.5,-.5)
    X=np.column_stack([np.ones(len(age)),reserve,age,age*reserve]); inv=np.linalg.inv(X.T@X)
    model=[]; strata=[]
    for metric in metrics+['gene_FSHR','gene_CYP19A1','energy_OXPHOS','energy_glycolysis','energy_TCA','energy_equal3']:
        y=primary[metric].to_numpy(); beta=inv@X.T@y; err=y-X@beta; lev=np.sum((X@inv)*X,axis=1)
        cov=inv@(X.T@((err/(1-lev))[:,None]**2*X))@inv
        for label,L in [('reserve',[0,1,0,0]),('age',[0,0,1,0]),('interaction',[0,0,0,1]),('young_DOR-NOR',[0,1,0,-.5]),('old_DOR-NOR',[0,1,0,.5])]:
            L=np.array(L); est=L@beta; se=np.sqrt(L@cov@L); df=len(y)-4; ci=stats.t.ppf(.975,df)*se
            row=dict(metric=metric,term=label,effect=est,se=se,low=est-ci,high=est+ci,p=2*stats.t.sf(abs(est/se),df),n=len(y))
            (model if label in ['reserve','age','interaction'] else strata).append(row)
    m=pd.DataFrame(model); sel=m.metric.isin(metrics); p=m.loc[sel,'p'].to_numpy(); order=np.argsort(p); q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1); temp=np.empty(len(p));temp[order]=q;m.loc[sel,'q30']=temp
    m.to_csv(P/'stage12_primary_model.csv',index=False); pd.DataFrame(strata).to_csv(P/'stage12_primary_strata.csv',index=False)
    extrows=[]
    for ds in ['GSE193136','GSE232306','GSE274832','GSE191322']:
        a=s[s.dataset.eq(ds)]; xname,yname=('COS','mild') if ds=='GSE191322' else ('DOR','NOR')
        for met in metrics+['gene_FSHR','gene_CYP19A1','energy_OXPHOS','energy_equal3']:
            x=a.loc[a.group.eq(xname),met].dropna().to_numpy(); y=a.loc[a.group.eq(yname),met].dropna().to_numpy()
            diff=x.mean()-y.mean(); vx=x.var(ddof=1)/len(x); vy=y.var(ddof=1)/len(y); se=np.sqrt(vx+vy); df=(vx+vy)**2/(vx**2/(len(x)-1)+vy**2/(len(y)-1)); ci=stats.t.ppf(.975,df)*se
            extrows.append(dict(dataset=ds,contrast=xname+'-'+yname,metric=met,effect=diff,low=diff-ci,high=diff+ci,p=2*stats.t.sf(abs(diff/se),df),nx=len(x),ny=len(y)))
    pd.DataFrame(extrows).to_csv(P/'stage12_external_effects.csv',index=False)
    cor=[]
    for (ds,g),a in s[s.dataset.ne('GSE191322')].groupby(['dataset','group']):
        for x,y in [('A','D'),('D-E','energy_OXPHOS')]:
            rho,p=stats.spearmanr(a[x],a[y]);cor.append(dict(dataset=ds,group=g,x=x,y=y,n=len(a),rho=rho,p=p if len(a)>=5 else np.nan))
    pd.DataFrame(cor).to_csv(P/'stage12_correlations.csv',index=False)
    print('completed',len(s),'samples')
if __name__=='__main__': main()
