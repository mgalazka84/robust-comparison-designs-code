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
python revision_study.py --verify
```

The checks work in a fresh checkout with no saved outputs. They independently
verify the eight-core catalogue, all 104 nonisomorphic subdivisions through
eight objects by full rational deletion and inversion, four labelled orbit
counts, the complete five-object benchmark, and analytic graph families.
The statistical check also verifies the least-squares covariance identity and
Gaussian interval calibration. The additional subdivision checks verify four
core certificates, 888 endpoint quadratics, 152 resistance-interpolation identities,
and direct deletion profiles for uniform subdivisions of lengths 1–4. Exact
checks of random-loss and post-loss A calculations cover every candidate through
eight objects. Figure generation checks the displayed graph
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
python revision_study.py --nmax 20
python make_revision_material.py
```

Setting `OPENBLAS_NUM_THREADS=1` can reduce numerical-library overhead.
The simulation defaults to six scenarios with 50,000 replications each,
master seed 20261025, and 20 independent variance-calibration degrees of
freedom. Those 20 reference observations are additional to the graph budget.
The criterion study evaluates all 327 nonisomorphic simple bridgeless designs
with n+1 comparisons, for 4–20 objects, using rational arithmetic. Expected-loss
criteria assume a uniformly chosen missing edge. Nominal benchmarks in this
study use the same bridgeless competitor class. The numerical study covers a
finite range; the subdivision theorem requires its explicit certificate.
Every observation is synthetic, and the lost edge is chosen before responses
are generated. No external dataset or credentials are required.

## Programs

| File | Purpose |
|---|---|
| `research.py` | Graph constructors, resistance calculations, labelled search and exact benchmarks |
| `cores.py` | Complete eight-core enumeration and rational endpoint scoring |
| `statistical_study.py` | Equal-budget A/D/maximum-variance benchmark and Gaussian calibration |
| `revision_study.py` | Exact symmetric-core certificates and complete comparisons of seven criteria at budget n+1 |
| `make_revision_material.py` | Certificate matrices, criterion tables and figures, and a calibration planning example |
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
