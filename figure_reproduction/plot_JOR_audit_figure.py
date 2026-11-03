"""Descriptive technical audit only; existing Python backend, no tests or models."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
F=Path(__file__).resolve().parent;T=F.parent/'audit_inputs'
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
r=pd.read_csv(T/'Mapping_S1_row_disposition.csv');s=pd.read_csv(T/'Mapping_S6_sample_count_retention.csv')
cats=['retained','below_historical_expression_filter','queried_but_no_persisted_symbol']
summ=r.groupby('exclusion_category').agg(rows=('gene_id','size'),counts=('total_counts_116','sum')).reindex(cats)
summ.to_csv(F/'Audit_Figure_A1_disposition_source.csv')
s.to_csv(F/'Audit_Figure_A1_sample_source.csv',index=False)
fig,(a,b)=plt.subplots(1,2,figsize=(170/25.4,112/25.4),gridspec_kw={'width_ratios':[1,1.25]})
fig.subplots_adjust(left=.09,right=.98,bottom=.27,top=.83,wspace=.35)
colors=['#4A6E8A','#CBD3D8','#C95B5B'];left=np.zeros(2)
for cat,c,label in zip(cats,colors,['Retained','Below expression filter','Queried, no symbol']):
    vals=np.array([summ.loc[cat,'rows']/summ.rows.sum(),summ.loc[cat,'counts']/summ.counts.sum()])*100
    a.barh([1,0],vals,left=left,color=c,label=label,height=.45);left+=vals
a.set(yticks=[0,1],yticklabels=['Counts','Rows'],xlim=(0,100),xlabel='Share of the original matrix (%)',ylim=(-.6,1.6))
a.text(0,1.05,'a  Rows and counts are distinct',transform=a.transAxes,fontweight='bold')
a.text(1,1.28,f"Retained rows: {summ.loc['retained','rows']:,.0f} / {summ.rows.sum():,.0f}",fontsize=7)
a.text(1,-.42,f"Pooled retained counts: {100*summ.loc['retained','counts']/summ.counts.sum():.3f}%",fontsize=7)
a.legend(loc='upper left',bbox_to_anchor=(-.12,-.23),fontsize=7,frameon=False)
groups=['NOR_young','DOR_young','NOR_old','DOR_old'];rng=np.random.default_rng(1501);display=[]
for x,g in enumerate(groups):
    z=s[s.group.eq(g)];j=rng.uniform(-.16,.16,len(z));col='#C95B5B' if g.startswith('DOR') else '#4A6E8A'
    keep=z.included_112.to_numpy();v=z.count_retention_fraction.to_numpy()*100
    b.scatter((x+j)[keep],v[keep],s=13,c=col,alpha=.8,lw=0)
    b.scatter((x+j)[~keep],v[~keep],s=35,c=col,marker='x',lw=1.1)
    b.hlines(z.loc[z.included_112,'count_retention_fraction'].mean()*100,x-.24,x+.24,color='#263D4A',lw=1.7)
    display.extend(dict(sample=ss,display_x=float(xx),retained_pct=float(yy)) for ss,xx,yy in zip(z['sample'],x+j,v))
b.set(xticks=range(4),xticklabels=['NOR\n33 + 1','DOR\n25 + 1','NOR\n18 + 2','DOR\n36 + 0'],ylim=(99.3,99.8),ylabel='Count retention (%)')
b.text(0,1.05,'b  Sample-level count retention',transform=b.transAxes,fontweight='bold')
b.text(.18,-.22,'Younger',ha='center',transform=b.transAxes);b.text(.83,-.22,'Older',ha='center',transform=b.transAxes)
b.text(0,-.36,'Dots: primary 112; crosses: excluded 4\nLabels: primary + excluded n; bars: primary means',transform=b.transAxes,fontsize=7)
for ext in ['pdf','svg','tiff']:
    kw={'dpi':600} if ext in ['png','tiff'] else {}
    if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
    fig.savefig(F/f'Audit_Figure_A1.{ext}',**kw)
fig.savefig(F/'Audit_Figure_A1_preview.png',dpi=150)
fig.savefig(F/'Audit_Figure_A1_for_supplement.png',dpi=600)
pd.DataFrame(display).to_csv(F/'Audit_Figure_A1_display_positions.csv',index=False)
(F/'Audit_Figure_A1_legend.md').write_text('''Audit Figure A1. Historical expression filtering, symbol mapping and count retention. (a) Fractions of the 58,735 source rows and pooled counts across all 116 samples assigned to retained rows (15,417), below-filter rows not queried (42,600) and queried identifiers without a persisted symbol (718). Count fractions are pooled across samples; they are distinct from the median of sample-level retention fractions. (b) Every sample is shown on a deliberately expanded 99.3–99.8% y-axis. Filled dots mark the 112 primary samples and crosses the four samples excluded by mitochondrial fraction >20%. Horizontal bars show primary group means. X-axis n denotes primary plus excluded samples. Horizontal jitter is display-only and recorded in a separate table. The primary median retention was 99.637% (range 99.346–99.710%). This is a descriptive technical audit, with no group hypothesis test or alternative normalization fit. Source data: Audit_Figure_A1_disposition_source.csv and Audit_Figure_A1_sample_source.csv; full gene-level audit in Mapping_S1_row_disposition.csv.
''',encoding='utf8')
print('Exported descriptive mapping audit figure')
