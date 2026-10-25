"""Fast, independent checks for the mathematical research companion.

This does not repeat the complete 20-object core search or seven-object
labelled search. Run those documented commands for complete reproduction.
"""
import csv
from collections import defaultdict
from pathlib import Path
from fractions import Fraction
import json
import tempfile

import cores
import research

ROOT = Path(__file__).resolve().parent


def main():
    # A fresh code-only checkout has no generated result files. The core
    # verification and small exhaustive search below always run independently.
    result_paths = [ROOT / 'results' / name for name in
                   ['core_scores.csv', 'core_search.csv', 'core_optimizers.json']]
    has_results = all(path.is_file() for path in result_paths)
    if any(path.is_file() for path in result_paths) and not has_results:
        raise RuntimeError('Incomplete saved core results; regenerate all three files with cores.py.')
    if has_results:
        stored = list(csv.DictReader(result_paths[0].open()))
        groups = defaultdict(list)
        for row in stored:
            groups[int(row['n'])].append(Fraction(row['robust_value']))
        for row in csv.DictReader(result_paths[1].open()):
            n = int(row['n'])
            assert len(groups[n]) == int(row['nonisomorphic_candidates'])
            assert min(groups[n]) == Fraction(row['minimum'])
            assert min(groups[n]) >= Fraction(n-1, 3)
        assert len(stored) == 20133
    with tempfile.TemporaryDirectory(prefix='comparison-design-checks-') as d:
        out = Path(d)
        checks = cores.verify(out)
        assert checks['nonisomorphic_cores'] == 8
        assert checks['direct_rational_graph_checks'] == 104
        assert len(checks['independent_labelled_count_checks']) == 4
        research.exact_tradeoff(out)
        actual = cores.search(4, 8, out)
        expected = [(4, 1, '1'), (5, 3, '13/8'), (6, 11, '5/3'),
                    (7, 29, '35/16'), (8, 60, '17/7')]
        assert len(actual) == len(expected)
        for a, (n, count, minimum) in zip(actual, expected):
            assert a['n'] == n and a['minimum'] == minimum
            assert a['nonisomorphic_candidates'] == count
        winner_count = 0
        saved_winners = json.loads(result_paths[2].read_text()) if has_results else []
        for row in saved_winners:
            for winner in row['optimizers']:
                n, edges, paths = cores.subdivide(winner['core'], winner['lengths'])
                # Full direct rational inversions, independent of endpoint scoring.
                assert research.score_exact(n, edges) == Fraction(row['minimum'])
                winner_count += 1
        for n in range(4, 31):
            for lengths in research.optimal_theta(n):
                nn, edges = research.theta(*lengths)
                assert nn == n
                assert abs(research.score_direct(n, edges)
                           - float(research.predicted_optimum(n))) < 1e-8
        print(f'PASS: core catalogue, 104 complete rational cross-checks, four labelled '
              f'counts, the 120-graph tradeoff, {winner_count} exact core winners, '
              f'uniform K4 and theta families.')
        if not has_results:
            print('Saved-result checks: not run (no generated archive). Run cores.py '
                  '--nmax 20 --out results, then verify.py, to enable them.')


if __name__ == '__main__':
    main()
