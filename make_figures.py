"""Exact graph drawings and analytic performance curves for the research note."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from research import theta,score_fast,predicted_optimum,cycle_radius
from cores import subdivide

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.bbox':'tight','pdf.fonttype':42})
blue='#175B81';red='#BB423A';green='#237A68'

n=8
cycle=[(i,(i+1)%n) for i in range(n)]
angle=2*np.pi*np.arange(n)/n
circle=np.column_stack([np.cos(angle),np.sin(angle)])
_,te=theta(3,3,3)
tp=np.array([[-1.15,0],[1.15,0],[-.45,1],[.45,1],[-.4,0],[.4,0],[-.45,-1],[.45,-1]])
cases=[('Cycle: 8 comparisons',cycle,circle,'7'),
       ('Cycle + chord: 9 comparisons',cycle+[(0,4)],circle,'21/5'),
       ('Theta: 9 comparisons',te,tp,'7/2')]
fig,axs=plt.subplots(1,3,figsize=(11.6,3.35))
for ax,(title,edges,pos,value) in zip(axs,cases):
    score,wit,_,_=score_fast(n,edges)
    ix,i,j=wit
    for k,(u,v) in enumerate(edges):
        ax.plot(pos[[u,v],0],pos[[u,v],1],color=red if k==ix else blue,
                ls='--' if k==ix else '-',lw=2,zorder=1)
    ax.scatter(pos[:,0],pos[:,1],s=105,color='white',edgecolors=blue,lw=1.4,zorder=2)
    for k,(x,y) in enumerate(pos):
        ax.text(x,y,str(k+1),ha='center',va='center',fontsize=7.6,zorder=3)
    ax.set_title(title,fontsize=10.5,pad=15)
    ax.text(.5,-.12,rf'Worst variance after loss: ${value}\,\sigma^2$',
            transform=ax.transAxes,ha='center',fontsize=10)
    ax.set_aspect('equal');ax.axis('off');ax.set_xlim(-1.4,1.4);ax.set_ylim(-1.3,1.3)
fig.subplots_adjust(wspace=.35,bottom=.18)
fig.savefig(OUT/'designs.pdf');fig.savefig(OUT/'designs.png',dpi=180);plt.close(fig)

ns=np.arange(4,101)
opt=np.array([float(predicted_optimum(int(n)))/(n-1) for n in ns])
chord=np.array([float((int(n)+1)//2-1+cycle_radius(int(n)//2+1))/(n-1) for n in ns])
fig,ax=plt.subplots(figsize=(6.8,3.5))
ax.plot(ns,opt,color=green,lw=1.8,label='Optimal redesign: n+1 comparisons')
ax.plot(ns,chord,color=blue,lw=1.8,label='Best cycle augmentation: n+1 comparisons')
ax.axhline(1,color='#777777',ls=':',label='Cycle baseline: n comparisons')
ax.set_xlabel('Number of objects n')
ax.set_ylabel('Worst variance / cycle worst variance')
ax.set_ylim(.42,1.06)
ax.legend(frameon=False,fontsize=8.5,loc='upper right',bbox_to_anchor=(1,.86))
ax.grid(alpha=.16);fig.tight_layout()
fig.savefig(OUT/'variance_ratio.pdf');fig.savefig(OUT/'variance_ratio.png',dpi=180)
plt.close(fig)

nn, edges, paths = subdivide('K4', (2,)*6)
cp = np.array([[0,1.25],[-1.25,-1],[1.25,-1],[0,-.25]])
pos = [p for p in cp]
for u,v in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
    pos.append((cp[u]+cp[v])/2)
pos = np.array(pos)
fig,axs=plt.subplots(1,2,figsize=(8.2,3.4))
for panel,ax in enumerate(axs):
    kept = list(range(nn)) if panel==0 else [i for i in range(nn) if i!=4]
    for ix,(u,v) in enumerate(edges):
        if panel==1 and ix<2:
            continue
        ax.plot(pos[[u,v],0],pos[[u,v],1],lw=1.8,
                color=red if panel==0 and ix==0 else blue,
                ls='--' if panel==0 and ix==0 else '-',zorder=1)
    ax.scatter(pos[kept,0],pos[kept,1],s=120,color='white',edgecolors=blue,lw=1.3,zorder=2)
    for i in kept:
        ax.text(*pos[i],str(i+1),ha='center',va='center',fontsize=7.5,zorder=3)
    ax.set_aspect('equal');ax.axis('off');ax.set_xlim(-1.45,1.45);ax.set_ylim(-1.2,1.45)
    ax.set_title('10 objects, 12 comparisons' if panel==0 else 'Remove one whole path: H',fontsize=10)
    ax.text(.5,-.22,r'$\Phi_1(G_2)=3$' if panel==0 else r'$\mathrm{diam}_R(H)=R_H(1,2)=2$',
            transform=ax.transAxes,ha='center',fontsize=11)
fig.subplots_adjust(bottom=.26,wspace=.3)
fig.savefig(OUT/'k4_subdivision.pdf');fig.savefig(OUT/'k4_subdivision.png',dpi=180)
plt.close(fig)
print('Created six figure files.')
