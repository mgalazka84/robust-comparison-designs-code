"""Regenerate the manuscript's simulation table from recorded JSON results."""
import json
from pathlib import Path
root = Path(__file__).resolve().parent
rows = json.loads((root/'results/calibration_study.json').read_text())
names = ['Cycle', 'Best chord', 'Theta', r'$K_4$ subdivision', r'Nominal $G_N$', r'Robust $G_R$']
text = r'''\begin{table}[ht]
\centering\small
\caption{Synthetic Gaussian study with 50,000 replications per row and 20
independent calibration degrees of freedom. $m$ excludes calibration.
$h/\sigma$ is the analytic expected maximum half-width of the simultaneous
intervals. Parentheses give Monte Carlo standard errors for coverage.}
\label{tab:simulation}
\begin{tabular}{lrrrrr}
\toprule
Design & $n$ & $m$ & $\ph$ & $h/\sigma$ & Joint coverage\\
\midrule
'''
for name, r in zip(names, rows):
    s = r['exact_post_loss_variance']
    value = s
    text += f"{name} & {r['n']} & {r['graph_comparisons']} & ${value}$ & {r['expected_max_halfwidth']:.3f} & {r['simultaneous_t_coverage']:.4f} ({r['coverage_mcse']:.4f}) " + r'\\' + '\n'
text += '\\bottomrule\n\\end{tabular}\n\\end{table}\n'
(root/'submission_sources').mkdir(exist_ok=True)
(root/'submission_sources/study_table.tex').write_text(text)
