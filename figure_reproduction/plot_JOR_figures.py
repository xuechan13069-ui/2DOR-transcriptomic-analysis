"""Publication figures from the accompanying full-precision source CSV tables.
Use the existing project Python environment. No statistical models are fitted here.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

F=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
 'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'savefig.facecolor':'white'})
W=170/25.4
C={'DOR':'#C95B5B','NOR':'#4A6E8A','ink':'#263D4A','light':'#F2F5F6','grey':'#69757C','protocol':'#78668C','line':'#B7C1C6'}
manifest=[]
superdigits=str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
def fmtprob(v):
    value=float(v)
    s=f'{value:.2e}' if value < .001 else f'{value:.3g}'
    if 'e' in s:
        a,b=s.split('e');return a+r'$\times 10^{'+str(int(b))+'}$'
    return s

def export(fig,name,description):
    fig.canvas.draw()
    for fmt in ['pdf','svg','tiff']:
        kw={'dpi':600} if fmt in ['png','tiff'] else {}
        if fmt=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(F/f'{name}.{fmt}',**kw)
    fig.savefig(F/f'{name}_preview.png',dpi=150)
    manifest.append({'figure':name,'width_mm':round(fig.get_figwidth()*25.4,1),
                     'height_mm':round(fig.get_figheight()*25.4,1),'description':description,
                     'formats':'PDF SVG TIFF600; PNG150 preview','synthetic_data':False})
    plt.close(fig)

def letter(ax,s,x=-.12,y=1.04):ax.text(x,y,s,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
def box(ax,x,y,w,h,label,fc=C['light'],size=8,weight='normal'):
    ax.add_patch(plt.Rectangle((x,y),w,h,linewidth=.55,edgecolor='#404040',facecolor='none'))
    ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=size,color=C['ink'],fontweight='normal',linespacing=1.35)
def arr(ax,p,q,dashed=False):
    ax.add_patch(FancyArrowPatch(p,q,arrowstyle='-|>',mutation_scale=9,lw=.85,color=C['grey'],linestyle='--' if dashed else '-'))

# Figure 1: all sample counts and node labels come from source tables.
co=pd.read_csv(F/'Figure_1_cohort_source.csv');nodes=pd.read_csv(F/'Figure_1_node_source.csv')
qc=pd.read_csv(F/'Figure_1_sample_QC_source.csv')
fig=plt.figure(figsize=(W,177/25.4));a=fig.add_axes([.06,.54,.90,.40]);a.set(xlim=(0,1),ylim=(0,1));a.axis('off')
letter(a,'a',0,1.03);a.text(.045,1.04,'Clinical comparisons',fontsize=10,fontweight='bold')
primary=co[co.dataset.eq('GSE286918')].set_index('group')
box(a,.01,.76,.52,.20,f"GSE286918 · {int(primary.original_n.sum())} granulosa-cell samples\nStimulated IVF/ICSI · age × reserve",size=8,weight='bold')
arr(a,(.27,.75),(.27,.65))
excluded=qc[~qc.included_112]
box(a,.01,.48,.52,.17,f"Primary: n={int(primary.analysed_n.sum())} · mitochondrial fraction ≤20%\nExcluded n={len(excluded)}; all 116 retained for sensitivity",size=8)
lines=['                     Younger       Older',
 f"DOR                {int(primary.loc['younger DOR','analysed_n'])}                 {int(primary.loc['older DOR','analysed_n'])}",
 f"NOR                {int(primary.loc['younger NOR','analysed_n'])}                 {int(primary.loc['older NOR','analysed_n'])}"]
box(a,.01,.06,.52,.32,'\n'.join(lines),size=9)
a.text(.01,.01,'Younger <35 years; older ≥35 years',fontsize=7,color=C['grey'])
extlines=['Cross-cohort directional assessment']
for ds in ['GSE193136','GSE232306','GSE274832']:
    n=co[co.dataset.eq(ds)].analysed_n.tolist();extlines.append(f'{ds} · DOR {n[0]} / NOR {n[1]}')
box(a,.60,.49,.39,.47,'\n\n'.join(extlines),size=8)
box(a,.60,.06,.39,.31,'Separate protocol comparison\nGSE191322 · poor responders\nConventional 5 / mild 5',size=8,fc='#F5F2F7')

b=fig.add_axes([.06,.07,.90,.39]);b.set(xlim=(-.02,1.02),ylim=(0,1));b.axis('off')
letter(b,'b',0,1.05);b.text(.045,1.06,'Five-node transcriptional framework',fontsize=10,fontweight='bold')
box(b,.00,.55,.23,.25,'A · De novo\nsterol synthesis\n22 frozen genes',size=8)
box(b,.00,.15,.23,.25,'B · Uptake and\nhomeostasis\n5 genes',size=8)
box(b,.32,.35,.18,.25,'C · Intracellular /\nmitochondrial\ndelivery\nNPC1 · NPC2 · STAR',size=6.4)
box(b,.57,.35,.18,.25,'D · Early enzymes\nCYP11A1\nHSD3B1 · HSD3B2',size=7.5)
box(b,.82,.30,.18,.34,'E · FSH–aromatase-\nassociated terminal\ntranscriptional\nprogramme\nFSHR · NR5A1\nCYP19A1',size=6.0)
arr(b,(.235,.67),(.315,.54));arr(b,(.235,.28),(.315,.41));arr(b,(.505,.475),(.565,.475));arr(b,(.755,.475),(.815,.475),True)
box(b,.71,.80,.29,.15,'Theca-derived androgen supply',size=7)
arr(b,(.91,.79),(.91,.65),True)
box(b,.31,.02,.68,.17,'Derived coordinates:   A−C     A−D     D−E     A−E     B−C',size=7.3)
b.text(.01,-.09,'Solid: adjacent framework steps     Dashed: indirect or external supply',fontsize=6.6,color=C['grey'])
export(fig,'Fig1','Design and five-node transcriptional framework')

# Figure 2: main reserve effects, full five-node comparison.
df=pd.read_csv(F/'Figure_2_source.csv').set_index('metric').loc[list('ABCDE')]
fig=plt.figure(figsize=(W,97/25.4));gs=fig.add_gridspec(1,3,width_ratios=[1.25,2.05,2.1],left=.035,right=.98,bottom=.23,top=.85,wspace=.08)
labels=fig.add_subplot(gs[0]);ax=fig.add_subplot(gs[1]);numbers=fig.add_subplot(gs[2])
for aa in [labels,ax,numbers]:aa.set_ylim(-.65,4.8)
labels.axis('off');numbers.axis('off')
names=['A  Sterol synthesis','B  Uptake / homeostasis','C  Delivery','D  Early enzymes','E  Terminal programme']
for j,(met,r) in enumerate(df.iterrows()):
    y=4-j;colour=C['ink'] if met in ['A','D'] else C['grey']
    labels.text(0,y,names[j],va='center',fontsize=7.3,fontweight='bold' if met in ['A','D'] else 'normal')
    ax.errorbar(r.effect,y,xerr=[[r.effect-r.low],[r.high-r.effect]],fmt='o',color=colour,capsize=3,ms=5,lw=1)
    numbers.text(.01,y,f'{r.effect:+.3f} [{r.low:+.3f}, {r.high:+.3f}]',va='center',fontsize=7.5)
    numbers.text(.85,y,fmtprob(r.q30),va='center',ha='center',fontsize=7.5)
numbers.text(.01,4.65,'Estimate [95% CI]',fontsize=7.3,fontweight='bold');numbers.text(.85,4.65,'BH q',fontsize=7.3,ha='center',fontweight='bold')
ax.axvline(0,color='#D2D8DC',ls='--',lw=.65);ax.set_yticks([]);ax.spines['left'].set_visible(False);ax.set_xlim(-.65,1.25)
ax.set_xlabel('DOR − NOR\n(mean gene-standardized score)',fontsize=8)
export(fig,'Fig2','Five node reserve coefficients and HC3 confidence intervals')

# Figure 3: actual continuous scores; display-only jitter.
sc=pd.read_csv(F/'Figure_3_patient_source.csv');effects=pd.read_csv(F/'Figure_3_effect_source.csv')
fig=plt.figure(figsize=(W,151/25.4));gs=fig.add_gridspec(2,2,height_ratios=[2.5,1],left=.10,right=.98,bottom=.10,top=.91,wspace=.30,hspace=.55)
order=['NOR_young','DOR_young','NOR_old','DOR_old'];pos=[0,1,3,4];rng=np.random.default_rng(1403)
maxv=max(sc['A-E'].max(),sc['D-E'].max());minv=min(sc['A-E'].min(),sc['D-E'].min())
jitter_rows=[]
for col,met in enumerate(['A-E','D-E']):
    ax=fig.add_subplot(gs[0,col]);letter(ax,chr(97+col),-.17,1.02);ax.set_title(met.replace('-','−'),loc='left',fontweight='bold')
    arrays=[sc.loc[sc.group.eq(g),met].to_numpy() for g in order]
    bp=ax.boxplot(arrays,positions=pos,widths=.58,patch_artist=True,showfliers=False,
                  medianprops={'color':C['ink'],'linewidth':1.2},boxprops={'linewidth':.7},whiskerprops={'linewidth':.7},capprops={'linewidth':.7})
    ticks=[]
    for g,x,ar,patch in zip(order,pos,arrays,bp['boxes']):
        colour=C['DOR'] if g.startswith('DOR') else C['NOR'];patch.set(facecolor='none',edgecolor='#505050',linewidth=.7)
        sub=sc[sc.group.eq(g)];jit=rng.uniform(-.18,.18,len(sub));ax.scatter(x+jit,ar,s=12,c=colour,marker='s' if g.startswith('DOR') else 'o',alpha=.56,lw=0,zorder=3)
        for sample,j,y in zip(sub['sample'],jit,ar):jitter_rows.append({'sample':sample,'group':g,'metric':met,'display_x':x+j,'score':y})
        ticks.append(f"{g.split('_')[0]}\nn={len(ar)}")
    ax.set_xticks(pos,ticks);ax.set_xlim(-.65,4.65);ax.set_ylim(minv-.25,maxv+.25)
    ax.axhline(0,color=C['line'],lw=.7,ls='--',zorder=0);ax.set_ylabel('Within-sample programme difference')
    ax.text(.22,-.20,'Younger',transform=ax.transAxes,ha='center',fontsize=8);ax.text(.79,-.20,'Older',transform=ax.transAxes,ha='center',fontsize=8)
    bot=fig.add_subplot(gs[1,col]);letter(bot,chr(99+col),-.17,1.02)
    for y,term,label in [(1,'young_DOR-NOR','Younger'),(0,'old_DOR-NOR','Older')]:
        r=effects[(effects.metric==met)&(effects.term==term)].iloc[0]
        bot.errorbar(r.effect,y,xerr=[[r.effect-r.low],[r.high-r.effect]],fmt='o',color=C['ink'],capsize=3,ms=4)
        bot.text(.98,y+.27,f'{r.effect:+.3f} [{r.low:+.3f}, {r.high:+.3f}]',transform=bot.get_yaxis_transform(),ha='right',fontsize=6.8)
    bot.axvline(0,color=C['line'],ls='--',lw=.7);bot.set(xlim=(-.25,1.85),ylim=(-.5,1.8),yticks=[0,1],yticklabels=['Older','Younger'])
    bot.set_xlabel('Age-stratified DOR − NOR [95% CI]');bot.spines['left'].set_visible(False)
pd.DataFrame(jitter_rows).to_csv(F/'Figure_3_display_positions.csv',index=False)
export(fig,'Fig3','Patient distributions and age-stratified contrasts')

# Figure 4: separate label, estimate and confidence interval columns.
data=pd.read_csv(F/'Figure_4_source.csv')
dslist=['GSE286918','GSE193136','GSE232306','GSE274832','GSE191322']; ys=[5,4,3,2,0]
fig=plt.figure(figsize=(W,180/25.4))
gs=fig.add_gridspec(2,3,width_ratios=[1.30,2.35,2.2],left=.035,right=.98,bottom=.15,top=.96,wspace=.12,hspace=.52)
for i,met in enumerate(['A-E','D-E']):
    lab=fig.add_subplot(gs[i,0]); ax=fig.add_subplot(gs[i,1]); num=fig.add_subplot(gs[i,2])
    for a in [lab,num]: a.axis('off'); a.set_ylim(-.7,5.8)
    lab.text(0,6.15,chr(97+i)+'  '+met.replace('-','−'),fontsize=10,fontweight='bold')
    num.text(0,6.15,'Estimate [95% CI]',fontsize=8)
    for ds,y,n in zip(dslist,ys,['DOR 61 / NOR 51','DOR 6 / NOR 6','DOR 6 / NOR 6','DOR 3 / NOR 3','COS 5 / mild 5']):
        r=data[(data.dataset==ds)&(data.metric==met)].iloc[0]
        colour=C['protocol'] if ds=='GSE191322' else C['DOR'] if ds=='GSE274832' else '#174A70' if ds=='GSE286918' else '#454F55'
        lab.text(0,y,ds,fontsize=8,va='center',fontweight='bold' if ds=='GSE286918' else 'normal')
        lab.text(0,y-.26,n,fontsize=7,color=C['grey'],va='center')
        marker='D' if ds=='GSE191322' else 's' if ds=='GSE286918' else 'o'
        ax.errorbar(r.effect,y,xerr=[[r.effect-r.low],[r.high-r.effect]],fmt=marker,color=colour,ms=4,capsize=2,lw=.75)
        num.text(.02,y,f'{r.effect:+.3f} [{r.low:+.3f}, {r.high:+.3f}]',va='center',fontsize=7.5,color=colour)
    lab.text(0,.65,'Protocol comparison',fontsize=7,color=C['protocol'])
    for a in [lab,ax,num]: a.axhline(.95,ls=':',lw=.6,color=C['line'])
    ax.axvline(0,ls='--',lw=.6,color=C['line'])
    ax.set(ylim=(-.7,5.8),yticks=[],xlim=(-6,4) if i==0 else (-3,2)); ax.spines['left'].set_visible(False)
    ax.set_xlabel('Study-specific score difference')
handles=[Line2D([],[],marker='s',color='#174A70',linestyle='None',markersize=4,label='Primary: factorial model, HC3'),
         Line2D([],[],marker='o',color='#454F55',linestyle='None',markersize=4,label='External: Welch DOR − NOR'),
         Line2D([],[],marker='D',color=C['protocol'],linestyle='None',markersize=4,label='Protocol: Welch COS − mild')]
fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.035,.042),ncol=3,fontsize=7,handletextpad=.3,columnspacing=1,frameon=False)
fig.text(.04,.018,'Scores standardized within each study; effect magnitudes are not directly comparable across studies.',fontsize=6.8,color=C['grey'])
export(fig,'Fig4','Cross-cohort directional assessment and separate protocol comparison')
