"""Render exact subdivision certificates and criterion-comparison material."""
from pathlib import Path
from fractions import Fraction as F
import json
import math
import itertools as it
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from revision_study import symmetric_cores, exact_metrics, certificates
from make_explanatory_figures import graph, save, BLUE, TEAL, RED

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'submission_sources'; OUT.mkdir(exist_ok=True)
DATA=json.loads((ROOT/'results/criterion_comparisons.json').read_text())
CERT=certificates()


def frac(x):
    x=F(x)
    return str(x.numerator) if x.denominator==1 else str(x.numerator)+'/'+str(x.denominator)


def write_certificates():
    lines=[r'''The core vertices are numbered from zero. For $K_{3,3}$ the
parts are $\{0,1,2\}$ and $\{3,4,5\}$. For $Q_3$ the labels are binary
integers from 0 to 7 and adjacency means differing in one bit. The Petersen
vertices are the ten two-element subsets of $\{0,1,2,3,4\}$ in lexicographic
order; two vertices are adjacent if the subsets are disjoint. The deleted
edge is stated for each matrix. Let $R$ denote the displayed resistance
matrix of $F=H-uv$. For $i,j\ne0$, put
$K_{ij}=(R_{i0}+R_{j0}-R_{ij})/2$ and ground vertex 0 in the Laplacian $L$.
Direct multiplication gives $LK=I$ in each case, certifying every entry.
The endpoint tests use $A=R_{wa}$, $B=R_{wb}$ and $c=1-R_{ab}$ for every
remaining edge $ab$ and both $w=u,v$. Each quadratic is bounded by its
endpoint values and its stationary value when that point lies in $(0,1)$.
The table below reports the maximum over all these tests. This is an exact
finite verification of the hypotheses used in the all-length proof.''',
        r'\begin{center}\begin{tabular}{lrrrr}',r'\toprule',
        r'Core & Deleted edge & $a_0$ & $\max R_{ij}$ & Endpoint maximum\\',r'\midrule']
    names={'K4':r'$K_4$','K33':r'$K_{3,3}$','Q3':r'$Q_3$','Petersen':'Petersen'}
    for x in CERT:
        lines.append(f"{names[x['core']]} & $({x['edge'][0]},{x['edge'][1]})$ & ${x['alpha']}$ & ${x['diameter']}$ & ${x['endpoint_continuous_eccentricity']}$"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}']
    for x in CERT:
        lines += [r'\subsection*{'+names[x['core']]+r' resistance certificate}',
            f"For the deleted edge $({x['edge'][0]},{x['edge'][1]})$, the resistance matrix is",
            r'\[ R=\frac{1}{'+str(x['resistance_denominator'])+r'}\begin{pmatrix}']
        lines += [' & '.join(map(str,row))+r'\\' for row in x['integer_resistance_matrix']]
        lines += [r'\end{pmatrix}.\]']
    lines += [r'''The machine-readable certificate also lists every quadratic,
its maximum and an attaining point. Checking all core-edge deletions,
rather than using edge transitivity, gives 888 exact endpoint-quadratic
checks. Independently, direct numerical deletion and inversion checks the
16 uniform subdivisions obtained from the four simple cores with
$\ell=1,2,3,4$. Direct rational reinversion is also used for all four
cores at $\ell=1$ and for $K_4,K_{3,3}$ at $\ell=2$. The finite checks
verify the implementation; the interpolation and endpoint arguments prove
the results for arbitrary subdivision lengths.''']
    (OUT/'symmetric_certificates.tex').write_text('\n'.join(lines)+'\n')


def draw_certificate():
    fig,axes=plt.subplots(1,3,figsize=(7.4,2.9))
    for ax,name in zip(axes,['K33','Q3','Petersen']):
        p,es=symmetric_cores()[name]
        if name=='K33': pos=[(-.85,y) for y in [1,0,-1]]+[(.85,y) for y in [1,0,-1]]
        elif name=='Q3':
            # First and last four vertices are the two square layers.
            pos=[(-1,-1),(1,-1),(-1,1),(1,1),(-.5,-.5),(.5,-.5),(-.5,.5),(.5,.5)]
        else:
            pairs=list(it.combinations(range(5),2)); pos=[None]*10
            outer=[(0,1),(2,3),(0,4),(1,2),(3,4)]
            inner=[(2,4),(1,4),(1,3),(0,3),(0,2)]
            for k in range(5):
                a=math.pi/2+2*math.pi*k/5
                pos[pairs.index(outer[k])]=(math.cos(a),math.sin(a))
                pos[pairs.index(inner[k])]=(.45*math.cos(a),.45*math.sin(a))
        graph(ax,es,pos,size=72)
        ax.set_xlim(-1.3,1.3); ax.set_ylim(-1.3,1.3)
        ax.set_title({'K33':r'$K_{3,3}$','Q3':r'$Q_3$','Petersen':'Petersen'}[name],fontsize=11)
    fig.subplots_adjust(wspace=.25,bottom=.27,top=.82)
    for x,label in zip([.235,.515,.795],[r'$r=4,\quad n=9\ell-3$',r'$r=5,\quad n=12\ell-4$',r'$r=6,\quad n=15\ell-5$']):
        fig.text(x,.16,label,ha='center',fontsize=10)
    fig.text(.5,.015,r'Every core edge becomes $\ell$ comparisons; the optimum is $(n-1)/r$.',ha='center',fontsize=9)
    save(fig,'subdivision_certificate')


def write_comparisons():
    selected={6:[(1,3,3),(2,2,3)],9:[(2,4,4),(3,3,4)],
              10:[(3,3,5),(3,4,4)],16:[(5,5,7),(5,6,6)]}
    lines=[r'\begin{table}[!htb]\centering',
        r'\caption{Exact comparisons at $m=n+1$. All rows are theta designs.',
        r'$U$ is expected post-loss MV under a uniformly chosen edge; $B$ is',
        r'expected post-loss A. At $n=6,9$ the two rows tie in $\Phi_1$.',
        r'At $n=10,16$ the first row is nominally MV-optimal and the second',
        r'is one-loss optimal within the bridgeless class.}\label{tab:criteria}',
        r'\begin{tabular}{rlrrrr}\toprule',
        r'$n$ & Branch lengths & $M$ & $\Phi_1$ & $U$ & $B$\\\midrule']
    for n,ls in selected.items():
        for lengths in ls:
            x=next(x for x in DATA['candidates'] if x['n']==n and x['kind']=='theta' and tuple(x['lengths'])==lengths)
            lines.append(f"{n} & $({','.join(map(str,lengths))})$ & ${x['M']}$ & ${x['Phi']}$ & ${x['expected_MV']}$ & ${x['expected_A']}$"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{table}']
    (OUT/'criterion_table.tex').write_text('\n'.join(lines)+'\n')
    lines=[r'\begin{center}\small\begin{longtable}{rrrrrr}',r'\toprule',
        r'$n$ & Candidates & $\min\Phi_1$ & Worst nominal MV eff. & Worst A eff. & Worst $U$ eff.\\\midrule\endhead']
    for x in DATA['summaries']:
        vals=[100*float(F(x['robust_efficiencies'][k][0])) for k in ['M','A','expected_MV']]
        lines.append(fr"{x['n']} & {x['candidates']} & ${x['optima']['Phi']}$ & {vals[0]:.2f}\% & {vals[1]:.2f}\% & {vals[2]:.2f}\%"+r'\\')
    lines += [r'\bottomrule\end{longtable}\end{center}']
    (OUT/'criterion_full_table.tex').write_text('\n'.join(lines)+'\n')


def draw_efficiencies():
    fig,axes=plt.subplots(1,2,figsize=(7.4,2.9))
    ns=[x['n'] for x in DATA['summaries']]
    for ax,key,title in zip(axes,['M','expected_MV'],['Nominal maximum-variance efficiency','Expected-loss efficiency']):
        lo=[100*float(F(x['robust_efficiencies'][key][0])) for x in DATA['summaries']]
        hi=[100*float(F(x['robust_efficiencies'][key][1])) for x in DATA['summaries']]
        ax.fill_between(ns,lo,hi,color=TEAL,alpha=.15)
        ax.plot(ns,hi,'o-',color=TEAL,ms=3,label='Best tied robust design')
        ax.plot(ns,lo,'s--',color=BLUE,ms=3,label='Worst tied robust design')
        ax.set_ylim(82,101.5); ax.set_xticks([4,8,12,16,20]); ax.set_xlabel('Objects n; comparisons n+1')
        ax.set_title(title,fontsize=9.4); ax.grid(alpha=.2); ax.set_ylabel('Efficiency (%)')
    axes[0].legend(fontsize=7.2,frameon=False,loc='lower right')
    fig.tight_layout(w_pad=1.7); save(fig,'criterion_efficiency')


def write_case():
    from research import theta
    n=8; cycle=[(i,(i+1)%n) for i in range(n)]
    rows=[]
    for name,es in [('Cycle',cycle),('Best added chord',cycle+[(0,4)]),('Redesign',theta(3,3,3)[1])]:
        met=exact_metrics(n,es)
        rows.append(dict(design=name,comparisons=len(es),**met))
    (ROOT/'results/calibration_case.json').write_text(json.dumps(rows,indent=2,default=str))
    lines=[r'\begin{table}[!htb]\centering',r'\caption{An eight-instrument planning example.',
       r'Variances are divided by $\sigma^2$. The final column sets',
       r'$\sigma=0.10$ measurement units and reports the worst post-loss',
       r'standard error. No observations are simulated.}\label{tab:case}',
       r'\begin{tabular}{lrrrrr}\toprule',r'Plan & Measurements & $M$ & $\Phi_1$ & $U$ & Worst s.e.\\\midrule']
    for x in rows:
        lines.append(f"{x['design']} & {x['comparisons']} & ${frac(x['M'])}$ & ${frac(x['Phi'])}$ & ${frac(x['expected_MV'])}$ & {0.1*math.sqrt(float(x['Phi'])):.4f}"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{table}']
    (OUT/'calibration_case.tex').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':
    write_certificates(); draw_certificate(); write_comparisons(); draw_efficiencies(); write_case()
    print('Created exact certificate, criterion comparison and calibration-case material.')
