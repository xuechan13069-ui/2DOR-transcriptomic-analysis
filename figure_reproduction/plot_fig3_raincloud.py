"""Three Fig3 alternatives. Input is the existing private patient-score table.
Density curves are descriptive KDEs; the age-stratified effects are unchanged.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde

OUT = Path(__file__).resolve().parent
SRC = OUT
scores = pd.read_csv(SRC / 'Figure_3_patient_source.csv')
effects = pd.read_csv(SRC / 'Figure_3_effect_source.csv')
assert len(scores) == 112 and scores['sample'].is_unique
order = ['NOR_young', 'DOR_young', 'NOR_old', 'DOR_old']
assert [sum(scores.group == g) for g in order] == [33, 25, 18, 36]
colors = {'NOR': '#4A6E8A', 'DOR': '#C95B5B'}
ink, faint = '#263D4A', '#B7C1C6'
plt.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':8, 'axes.spines.top':False, 'axes.spines.right':False,
    'axes.linewidth':.6, 'pdf.fonttype':42, 'svg.fonttype':'none'})
ylow = min(scores['A-E'].min(), scores['D-E'].min()) - .3
yhigh = max(scores['A-E'].max(), scores['D-E'].max()) + .3

def kde(vals):
    grid = np.linspace(vals.min(), vals.max(), 200)
    density = gaussian_kde(vals, bw_method='scott')(grid)
    return grid, density / density.max()

def summary(ax, vals, x, color):
    q1, med, q3 = np.quantile(vals, [.25,.5,.75])
    ax.plot([x,x], [q1,q3], color=color, lw=3.2, solid_capstyle='butt', zorder=5)
    ax.scatter([x],[med], color='white', edgecolor=color, linewidth=.8, s=22, zorder=6)

def top(ax, metric, mode):
    xs = [0,1,3,4] if mode != 'split' else [-.13,.13,1.87,2.13]
    rng = np.random.default_rng(314159 if metric == 'A-E' else 271828)
    if mode == 'beeswarm':
        tmp = scores.copy()
        tmp['group'] = pd.Categorical(tmp['group'], categories=order, ordered=True)
        sns.swarmplot(data=tmp, x='group', y=metric, order=order,
                      palette=[colors[g[:3]] for g in order], size=3.3, ax=ax, zorder=2)
        xs = [0,1,2,3]
    for i,(g,x) in enumerate(zip(order,xs)):
        vals = scores.loc[scores.group == g, metric].to_numpy()
        col = colors[g[:3]]
        if mode == 'raincloud':
            grid,d = kde(vals)
            ax.fill_betweenx(grid, x+.04, x+.04+.28*d, color=col, alpha=.30, lw=0)
            ax.scatter(x-.12-rng.uniform(.015,.13,len(vals)), vals, s=12, color=col, alpha=.64, linewidth=0)
            summary(ax, vals, x+.01, col)
        elif mode == 'beeswarm':
            mean = np.mean(vals); se = np.std(vals, ddof=1)/np.sqrt(len(vals))
            ax.plot([x,x], [mean-se,mean+se], color=ink, lw=1.2, zorder=5)
            ax.plot([x-.16,x+.16], [mean,mean], color=ink, lw=1.4, zorder=6)
            ax.plot([x-.07,x+.07], [mean-se,mean-se], color=ink, lw=.8,zorder=5)
            ax.plot([x-.07,x+.07], [mean+se,mean+se], color=ink, lw=.8,zorder=5)
        else:
            center = 0 if i<2 else 2
            side = -1 if g.startswith('NOR') else 1
            grid,d = kde(vals)
            ax.fill_betweenx(grid, center, center+side*.42*d, color=col, alpha=.35, lw=0)
            jitter = rng.uniform(.015,.14,len(vals))
            ax.scatter(center+side*(.17+jitter), vals, s=10, color=col, alpha=.62, linewidth=0)
            summary(ax, vals, center+side*.055, col)
    ax.axhline(0, color=faint, ls='--', lw=.7,zorder=0)
    ax.set_ylim(ylow,yhigh)
    ax.set_title(metric.replace('-', '−'), loc='left', weight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Within-sample programme difference')
    if mode == 'split':
        ax.set_xlim(-.63,2.63); ax.set_xticks([0,2], ['Younger\nNOR 33 · DOR 25','Older\nNOR 18 · DOR 36'])
    else:
        ax.set_xticks(xs, ['NOR\nn=33','DOR\nn=25','NOR\nn=18','DOR\nn=36'])
        ax.set_xlim(-.55, max(xs)+.55)
        ax.text(.23,-.20,'Younger',transform=ax.transAxes,ha='center')
        ax.text(.77,-.20,'Older',transform=ax.transAxes,ha='center')

def bottom(ax, metric):
    for y,term in [(1,'young_DOR-NOR'),(0,'old_DOR-NOR')]:
        row = effects[(effects.metric==metric)&(effects.term==term)].iloc[0]
        ax.errorbar(row.effect,y,xerr=[[row.effect-row.low],[row.high-row.effect]],
                    fmt='o',color=ink,capsize=3,ms=4,lw=1)
        ax.text(.98,y+.25,f'{row.effect:+.3f} [{row.low:+.3f}, {row.high:+.3f}]',
                transform=ax.get_yaxis_transform(),ha='right',fontsize=7.1)
    ax.axvline(0,color=faint,ls='--',lw=.7)
    ax.set(xlim=(-.5,2),ylim=(-.5,1.8),yticks=[0,1],yticklabels=['Older','Younger'])
    ax.set_xticks([-.5,0,.5,1,1.5,2]);ax.spines['left'].set_visible(False)
    ax.set_xlabel('Age-stratified DOR − NOR [95% CI]')

for mode in ['raincloud']:
    fig = plt.figure(figsize=(170/25.4,151/25.4))
    gs = fig.add_gridspec(2,2,height_ratios=[2.5,1],left=.10,right=.98,bottom=.10,top=.92,wspace=.30,hspace=.57)
    for col,metric in enumerate(['A-E','D-E']):
        ax=fig.add_subplot(gs[0,col]);top(ax,metric,mode)
        ax.text(-.17,1.02,chr(97+col),transform=ax.transAxes,weight='bold',fontsize=10)
        ax=fig.add_subplot(gs[1,col]);bottom(ax,metric)
        ax.text(-.17,1.02,chr(99+col),transform=ax.transAxes,weight='bold',fontsize=10)
    for ext in ['pdf','svg']:
        fig.savefig(OUT/f'Fig3.{ext}')
    fig.savefig(OUT/'Fig3_preview.png',dpi=160)
    fig.savefig(OUT/'Fig3.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)
print('Raincloud figure rendered from 112 existing patient scores; no model refitted.')
