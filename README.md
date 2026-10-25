# Robust comparison designs

Python programs for exact optimization of pairwise comparison designs after
one lost observation, with independent checks and a Gaussian calibration study.
This repository contains source code and software documentation. Numerical
outputs, tables and figures are generated locally by the commands below.

## Install and check

Use Python 3.12 and create a virtual environment:

```bash
python -m venv .venv
```

Activate it with `source .venv/bin/activate` on Linux/macOS, or
`.venv\Scripts\Activate.ps1` in Windows PowerShell. Then run:

```bash
python -m pip install -r requirements.txt
python verify.py
python statistical_study.py --verify --reps 5000
python make_explanatory_figures.py
```

The checks work in a fresh checkout with no saved outputs. They independently
verify the eight-core catalogue, all 104 nonisomorphic subdivisions through
eight objects by full rational deletion and inversion, four labelled orbit
counts, the complete five-object benchmark, and analytic graph families.
The statistical check also verifies the least-squares covariance identity and
Gaussian interval calibration. Figure generation checks the displayed graph
values in rational arithmetic.

If a complete `results/` archive is present, `verify.py` additionally checks
the 20,133 saved candidate scores and directly recomputes every saved optimizer
through 20 objects. It reports explicitly when these archive checks are skipped.

## Generate the full outputs

Run these commands from the repository directory. The complete exact search
is more expensive than the checks above.

```bash
python research.py --nmax 7 --out results
python cores.py --nmax 20 --out results
python statistical_study.py --out results
python verify.py
python make_tables.py
python make_figures.py
python make_study_table.py
python make_explanatory_figures.py
```

Setting `OPENBLAS_NUM_THREADS=1` can reduce numerical-library overhead.
The simulation defaults to six scenarios with 50,000 replications each,
master seed 20261025, and 20 independent variance-calibration degrees of
freedom. Those 20 reference observations are additional to the graph budget.
Every observation is synthetic, and the lost edge is chosen before responses
are generated. No external dataset or credentials are required.

## Programs

| File | Purpose |
|---|---|
| `research.py` | Graph constructors, resistance calculations, labelled search and exact benchmarks |
| `cores.py` | Complete eight-core enumeration and rational endpoint scoring |
| `statistical_study.py` | Equal-budget A/D/maximum-variance benchmark and Gaussian calibration |
| `verify.py` | Independent checks, with optional verification of regenerated full results |
| `make_figures.py` | Analytic performance plots and graph diagrams |
| `make_explanatory_figures.py` | Construction, failure-localization, core and precision diagrams |
| `make_tables.py`, `make_study_table.py` | Tables from regenerated outputs |

Outputs are written to `results/`, `figures/`, `tables/` and
`submission_sources/`; these directories are ignored by Git. Regeneration
overwrites the requested output files, so use a separate `--out` directory
for partial searches.

## License

The software is distributed under the MIT license in `LICENSE`.
OpenAI ChatGPT (Codex) assisted code development and debugging.
