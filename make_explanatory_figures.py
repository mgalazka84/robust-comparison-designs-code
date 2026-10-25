"""Reproducible explanatory diagrams from explicit graphs and exact calculations.

No saved data files are required. PDF and PNG outputs are generated locally.
The graph scores shown in the figures are checked by rational reinversion.
"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle
from scipy.stats import t

from research import theta, score_fast, score_exact, exact_resistance
from cores import CORES
from statistical_study import cases, scale_mean

OUT = Path(__file__).resolve().parent / 'figures'
BLUE, TEAL, RED, AMBER = '#205F86', '#087F78', '#B8393B', '#D48C00'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
    'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'axes.spines.top': False,
    'axes.spines.right': False})


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / (name + '.pdf'), pad_inches=.07)
    fig.savefig(OUT / (name + '.png'), dpi=180, pad_inches=.07)
    plt.close(fig)


def graph(ax, edges, pos, lost=None, pair=(), labels=None, size=125):
    pos = np.asarray(pos)
    for k, (u, v) in enumerate(edges):
        ax.plot(pos[[u, v], 0], pos[[u, v], 1],
                color=RED if k == lost else BLUE,
                ls=(0, (3, 2)) if k == lost else '-', lw=1.6, zorder=1)
    for i, (x, y) in enumerate(pos):
        ax.scatter(x, y, s=size, facecolors='#FFF0C7' if i in pair else 'white',
                   edgecolors=AMBER if i in pair else BLUE, linewidths=1.2, zorder=3)
        if labels is not None:
            ax.text(x, y, labels[i], ha='center', va='center', fontsize=7, zorder=4)
    ax.set_aspect('equal'); ax.axis('off')


def witness(n, edges, expected):
    assert score_exact(n, edges) == expected
    value, (ix, i, j), _, _ = score_fast(n, edges)
    assert abs(value - float(expected)) < 1e-9
    R = exact_resistance(n, edges[:ix] + edges[ix + 1:])
    assert R[i][j] == expected
    return ix, (i, j)


def theta_positions(lengths):
    pos = [(-1.15, 0), (1.15, 0)]
    for length, h in zip(lengths, (1., 0., -1.)):
        pos += [(-1.15 + 2.3 * j / length, h * math.sin(math.pi*j/length))
                for j in range(1, length)]
    return np.array(pos)


def opening():
    n = 8
    cycle = [(i, (i+1) % n) for i in range(n)]
    angle = 2*np.pi*np.arange(n)/n
    circle = np.column_stack([np.cos(angle), np.sin(angle)])
    _, th = theta(3, 3, 3)
    designs = [('Existing cycle', cycle, circle, Fraction(7), r'$7$'),
        ('Add one comparison', cycle+[(0, 4)], circle, Fraction(21, 5), r'$21/5=4.20$'),
        ('Redesign at the same budget', th, theta_positions((3, 3, 3)),
         Fraction(7, 2), r'$7/2=3.50$')]
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.55))
    for ax, (title, edges, pos, value, label) in zip(axes, designs):
        ix, pair = witness(n, edges, value)
        graph(ax, edges, pos, ix, pair)
        ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.2, 1.3)
        ax.set_title(title+'\n'+str(len(edges))+' measurements', fontsize=9.3, pad=7)
        ax.text(.5, -.06, label, transform=ax.transAxes, ha='center', fontsize=11)
    fig.subplots_adjust(wspace=.25, bottom=.24)
    fig.text(.5, .01, r'Largest pairwise variance after one loss, divided by $\sigma^2$',
             ha='center', fontsize=9)
    save(fig, 'opening_designs')


def tied_designs():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    for ax, lengths, title in zip(axes, [(3, 3, 4), (2, 4, 4)],
            ['Balanced branches', 'A second, unbalanced optimum']):
        n, edges = theta(*lengths)
        assert n == 9 and len(edges) == 10 and score_exact(n, edges) == Fraction(9, 2)
        graph(ax, edges, theta_positions(lengths))
        for h, length in zip([1.18, .22, -1.18], lengths):
            ax.text(0, h, str(length)+' edges', ha='center', va='center',
                    fontsize=9, color=TEAL if length == 4 else BLUE)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.45)
        ax.text(.5, -.08, r'$\Phi_1=9/2$; longest branch: 4 edges',
                transform=ax.transAxes, ha='center', fontsize=9)
    fig.subplots_adjust(wspace=.38, bottom=.16)
    save(fig, 'theta_ties')


def endpoint_loss():
    H = [(0, 1), (0, 2), (2, 3), (3, 1)]
    path = [(0, 4), (4, 5), (5, 6), (6, 7), (7, 1)]
    edges = H + path
    pos = [(0, 1), (0, -1), (-1.55, 1), (-1.55, -1),
           (1.05, 1.35), (2.05, .6), (2.05, -.6), (1.05, -1.35)]
    labels = ['u', 'v', 'w', 'z', r'$x_1$', r'$x_2$', r'$x_3$', r'$x_4$']
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1))
    for ax, cut, expected, title, sub in zip(axes, [6, 4],
          [Fraction(19, 4), Fraction(5)], ['Interior loss', 'Endpoint loss'],
          [r'Two remnants: $p=q=2$', 'One remnant: 4 edges']):
        R = exact_resistance(8, edges[:cut]+edges[cut+1:])
        D = max(map(max, R)); assert D == expected
        pair = next((i, j) for i in range(8) for j in range(i+1, 8) if R[i][j] == D)
        ax.add_patch(Rectangle((-1.8, -1.2), 2.05, 2.4,
                              facecolor='#F0F4F6', edgecolor='none', zorder=0))
        graph(ax, edges, pos, cut, pair, labels, size=145)
        ax.text(-.85, 0, r'$H$', fontsize=14, color='#687C86', ha='center')
        ax.set_title(title+'\n'+sub, fontsize=9.5, pad=7)
        ax.text(.5, -.07, r'$D(G-e)='+str(expected)+r'$', transform=ax.transAxes,
                ha='center', fontsize=10.5)
        ax.set_xlim(-2, 2.35); ax.set_ylim(-1.6, 1.7)
    fig.subplots_adjust(wspace=.17, bottom=.12)
    save(fig, 'endpoint_loss')


def draw_core(ax, name):
    p, edges = CORES[name]
    if p == 1: pos = np.array([(0., 0.)])
    elif p == 2: pos = np.array([(-.65, 0), (.65, 0)])
    elif p == 3: pos = np.array([(0, .55), (-.65, -.6), (.65, -.6)])
    elif name == 'K4': pos = np.array([(0, .8), (-.85, -.65), (.85, -.65), (0, -.12)])
    else: pos = np.array([(-.65, .6), (.65, .6), (-.65, -.6), (.65, -.6)])
    for (u, v), count in Counter(edges).items():
        if u == v:
            angles = [90, 210, 330] if count == 3 else ([180] if p == 2 and u == 0 else [0] if p == 2 else [90])
            for deg in angles:
                a = math.radians(deg)
                center = pos[u] + .44*np.array([math.cos(a), math.sin(a)])
                ax.add_patch(Ellipse(center, .64, .92, angle=deg-90,
                                    fill=False, edgecolor=BLUE, lw=1.5, zorder=1))
        else:
            radii = [0] if count == 1 else np.linspace(-.65 if count > 2 else -.2,
                                                       .65 if count > 2 else .2, count)
            for rad in radii:
                ax.add_patch(FancyArrowPatch(pos[u], pos[v], arrowstyle='-',
                    connectionstyle=f'arc3,rad={rad}', shrinkA=0, shrinkB=0,
                    linewidth=1.5, color=BLUE, zorder=1))
    ax.scatter(pos[:, 0], pos[:, 1], s=115, facecolors='white', edgecolors=BLUE, zorder=3)
    for i, point in enumerate(pos):
        ax.text(*point, str(i), fontsize=7, ha='center', va='center', zorder=4)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.35, 1.55)


def core_gallery():
    labels = [r'$B_3$', r'$P_4$', r'$LP_3$', r'$L_2P_2$',
              r'$LT_2$', r'$T_{221}$', r'$Q_{2211}$', r'$K_4$']
    fig, axes = plt.subplots(2, 4, figsize=(7.4, 3.9))
    for ax, name, label in zip(axes.flat, CORES, labels):
        draw_core(ax, name)
        p, edges = CORES[name]
        ax.set_title(label+'\n'+rf'$p={p},\ q={len(edges)}$', fontsize=9.2, pad=2)
    fig.subplots_adjust(wspace=.16, hspace=.22)
    save(fig, 'core_gallery')


def precision_tradeoff():
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.75),
                             gridspec_kw={'width_ratios': [1, 1, 1.25]})
    angle = np.pi/2+2*np.pi*np.arange(5)/5
    circle = np.column_stack([np.cos(angle), np.sin(angle)])
    robust_pos = [(-1, 0), (0, 1.1), (0, .48), (0, -1.05), (1, 0)]
    widths = []
    for ax, case, pos, title in zip(axes[:2], cases()[-2:], [circle, robust_pos],
            ['Nominal optimum', 'One-loss optimum']):
        name, n, edges, phi = case
        ix, pair = witness(n, edges, phi)
        graph(ax, edges, pos, ix, pair, [str(i) for i in range(5)], size=145)
        ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
        ax.set_title(title, fontsize=9.5)
        nominal = max(map(max, exact_resistance(n, edges)))
        ax.text(.5, -.06, r'$M='+str(nominal)+r'\quad\Phi_1='+str(phi)+r'$',
                transform=ax.transAxes, ha='center', fontsize=9)
        widths.append(t.isf(.05/(2*10), 20)*scale_mean(20)*math.sqrt(float(phi)))
    ax = axes[2]
    ax.barh([1, 0], widths, height=.43, color=[BLUE, TEAL])
    ax.set_yticks([1, 0], [r'$G_N$', r'$G_R$'])
    for y, w in zip([1, 0], widths): ax.text(w+.08, y, f'{w:.3f}', va='center', fontsize=9)
    ax.set_xlim(0, 5.4); ax.set_ylim(-.6, 1.65)
    ax.set_xticks([0, 2, 4]); ax.grid(axis='x', alpha=.15); ax.set_axisbelow(True)
    ax.set_title('Expected maximum\ninterval half-width', fontsize=9.5)
    ax.set_xlabel(r'In units of $\sigma$', fontsize=9)
    ax.text(.5, -.31, '9.86% smaller', color=TEAL, transform=ax.transAxes,
            ha='center', fontsize=10, fontweight='bold')
    assert abs(100*(1-widths[1]/widths[0])-9.861218) < .00001
    fig.subplots_adjust(wspace=.48, bottom=.23)
    save(fig, 'precision_tradeoff')


def main():
    opening(); tied_designs(); endpoint_loss(); core_gallery(); precision_tradeoff()
    print('PASS: exact figure values; created five explanatory figures (PDF and PNG).')


if __name__ == '__main__':
    main()
