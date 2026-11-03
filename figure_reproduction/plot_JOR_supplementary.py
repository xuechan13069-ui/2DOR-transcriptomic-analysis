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

F=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
 'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'savefig.facecolor':'white'})
W=170/25.4
C={'DOR':'#C95B5B','NOR':'#4A6E8A','ink':'#263D4A','light':'#F2F5F6','grey':'#69757C','protocol':'#78668C','line':'#B7C1C6'}
manifest=[]

def export(fig,name,description):
    fig.canvas.draw()
    for fmt in ['pdf','svg','png','tiff']:
        kw={'dpi':600} if fmt in ['png','tiff'] else {}
        if fmt=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(F/f'{name}.{fmt}',**kw)
    fig.savefig(F/f'{name}_preview.png',dpi=150)
    manifest.append({'figure':name,'width_mm':round(fig.get_figwidth()*25.4,1),
                     'height_mm':round(fig.get_figheight()*25.4,1),'description':description,
                     'formats':'PDF SVG PNG600 TIFF600','synthetic_data':False})
    plt.close(fig)

def letter(ax,s,x=-.12,y=1.04):ax.text(x,y,s,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
def box(ax,x,y,w,h,label,fc=C['light'],size=8,weight='normal'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008,rounding_size=0.012',lw=.65,ec=C['line'],fc=fc))
    ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=size,color=C['ink'],fontweight=weight,linespacing=1.35)
def arr(ax,p,q,dashed=False):
    ax.add_patch(FancyArrowPatch(p,q,arrowstyle='-|>',mutation_scale=9,lw=.85,color=C['grey'],linestyle='--' if dashed else '-'))

# Supplementary Figure S1: existing correlations, no fit lines or new tests.
sc=pd.read_csv(F/'Figure_S1_source.csv');cor=pd.read_csv(F/'Figure_S1_correlations.csv')
fig,axes=plt.subplots(2,2,figsize=(W,140/25.4),sharex=True,sharey=True)
fig.subplots_adjust(left=.1,right=.97,bottom=.10,top=.89,hspace=.35,wspace=.27)
for i,(ax,g) in enumerate(zip(axes.flat,['NOR_young','DOR_young','NOR_old','DOR_old'])):
    ss=sc[sc.group.eq(g)];r=cor[cor.group.eq(g)].iloc[0];colour=C['DOR'] if g.startswith('DOR') else C['NOR']
    ax.scatter(ss.A,ss.D,s=17,c=colour,alpha=.8,lw=0);letter(ax,chr(97+i),-.13,1.02)
    ax.set_title(g.replace('_',' · '),loc='left');ax.text(.03,.97,f'n={int(r.n)}; Spearman ρ={r.rho:.3f}',transform=ax.transAxes,va='top',fontsize=8)
    ax.axhline(0,lw=.5,c=C['line']);ax.axvline(0,lw=.5,c=C['line']);ax.set_xlabel('A · Sterol-synthesis score');ax.set_ylabel('D · Early-enzyme score')
    ax.tick_params(axis='both', labelbottom=True, labelleft=True)
fig.canvas.draw()
limits=[(ax.get_xlim(),ax.get_ylim()) for ax in axes.flat]
assert all(np.allclose(limits[0][0],x) and np.allclose(limits[0][1],y) for x,y in limits)
print('S1 shared x and y limits:',limits[0])
export(fig,'Supplementary_FigS1','Descriptive A-D associations in the primary four groups')
print('Exported four main figures and one supporting correlation figure')
