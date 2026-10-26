"""Render the exact search table from its committed, machine-readable results."""
import csv
from pathlib import Path
from fractions import Fraction
import json

ROOT = Path(__file__).resolve().parent
names = {'three_loops': 'B_3', 'four_parallel': 'P_4', 'loop_triple': 'LP_3',
         'two_loops_double': 'L_2P_2', 'loop_double_triangle': 'LT_2',
         'double_triangle': 'T_{221}', 'double_square': 'Q_{2211}', 'K4': 'K_4'}
rows = json.loads((ROOT / 'results/core_optimizers.json').read_text())
lines = [r'\begin{table}[ht]', r'\centering',
         r'\caption{Exact global optima for $m=n+2$. Lengths follow the ordered',
         r'core edges in Table~\ref{tab:cores}. These are rational exhaustive',
         r'search results; Theorem~\ref{thm:k4} also proves the indicated',
         r'uniform $K_4$ cases analytically. Ties count nonisomorphic optima.}\label{tab:exactcores}',
         r'\begin{tabular}{rrrclr}', r'\toprule',
         r'$n$ & Bound $(n-1)/3$ & Minimum $\Phi_1$ & Core & Path lengths & Ties\\',
         r'\midrule']
for row in rows:
    w = row['optimizers'][0]
    value = Fraction(row['minimum'])
    f = str(value.numerator) if value.denominator == 1 else f'{value.numerator}/{value.denominator}'
    lengths = ','.join(map(str, w['lengths']))
    bound = Fraction(row['n'] - 1, 3)
    lines.append(f"{row['n']} & ${bound}$ & ${f}$ & ${names[w['core']]}$ & $({lengths})$ & {len(row['optimizers'])}" + r'\\')
lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}', '']
(ROOT / 'tables').mkdir(exist_ok=True)
(ROOT / 'tables/core_results.tex').write_text('\n'.join(lines))
print('Created tables/core_results.tex from exact search output.')
