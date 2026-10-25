"""Complete core search for n+2 comparisons, with rational certificates.

Every feasible simple graph is a subdivision of one of eight bridgeless
multigraphs. The reduction and endpoint-deletion lemma are proved in the note.
No numerical pruning is used by the exact search.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from fractions import Fraction
from functools import lru_cache
import itertools as it
import json
import math
from pathlib import Path
import time

from research import connected, exact_resistance, score_exact, score_fast

CORES = {
    'three_loops': (1, ((0, 0),) * 3),
    'four_parallel': (2, ((0, 1),) * 4),
    'loop_triple': (2, ((0, 0),) + ((0, 1),) * 3),
    'two_loops_double': (2, ((0, 0), (0, 1), (0, 1), (1, 1))),
    'loop_double_triangle': (3, ((0, 0), (0, 1), (0, 2), (1, 2), (1, 2))),
    'double_triangle': (3, ((0, 1), (0, 1), (0, 2), (0, 2), (1, 2))),
    'double_square': (4, ((0, 1), (0, 1), (0, 2), (1, 3), (2, 3), (2, 3))),
    'K4': (4, tuple(it.combinations(range(4), 2))),
}


def mapped_edges(edges, perm):
    return tuple(sorted(tuple(sorted((perm[u], perm[v]))) for u, v in edges))


def core_key(p, edges):
    return min(mapped_edges(edges, perm) for perm in it.permutations(range(p)))


@lru_cache(None)
def automorphisms(name):
    p, edges = CORES[name]
    return tuple(perm for perm in it.permutations(range(p))
                 if mapped_edges(edges, perm) == edges)


def length_key(name, lengths):
    edges = CORES[name][1]
    return min(tuple(z[2] for z in sorted(
        (*sorted((perm[u], perm[v])), length)
        for (u, v), length in zip(edges, lengths)))
        for perm in automorphisms(name))


def automorphism_count(name, lengths):
    """Count full subdivision automorphisms, including loop reflections."""
    edges = CORES[name][1]
    original = tuple((*edge, length) for edge, length in zip(edges, lengths))
    vertex_maps = sum(tuple(sorted((*sorted((perm[u], perm[v])), length)
                     for (u, v), length in zip(edges, lengths))) == original
                     for perm in automorphisms(name))
    equal_branches = math.prod(math.factorial(count) for count in Counter(original).values())
    loop_flips = 2 ** sum(u == v for u, v in edges)
    return vertex_maps * equal_branches * loop_flips


def length_vectors(name, total):
    """One representative per isomorphism class of simple subdivisions."""
    edges = CORES[name][1]
    lower = [3 if u == v else 1 for u, v in edges]
    suffix = [0] * (len(edges) + 1)
    for i in range(len(edges) - 1, -1, -1):
        suffix[i] = suffix[i + 1] + lower[i]

    def visit(prefix, remaining):
        i = len(prefix)
        if i == len(edges):
            value = tuple(prefix)
            if remaining == 0 and length_key(name, value) == value:
                yield value
            return
        lo = lower[i]
        if i and edges[i] == edges[i - 1]:
            lo = max(lo, prefix[-1], 2 if edges[i][0] != edges[i][1] else 3)
        for length in range(lo, remaining - suffix[i + 1] + 1):
            yield from visit(prefix + [length], remaining - length)
    yield from visit([], total)


def subdivide(name, lengths):
    p, core_edges = CORES[name]
    edges, paths = [], []
    nxt = p
    for (u, v), length in zip(core_edges, lengths):
        path = [u] + list(range(nxt, nxt + length - 1)) + [v]
        nxt += length - 1
        paths.append(path)
        edges.extend(zip(path[:-1], path[1:]))
    assert len(set(tuple(sorted(e)) for e in edges)) == len(edges)
    assert all(u != v for u, v in edges)
    return nxt, edges, paths


def endpoint_indices(paths):
    selected, offset = set(), 0
    for path in paths:
        length = len(path) - 1
        selected.update((offset, offset + length - 1))
        offset += length
    return sorted(selected)


def exact_endpoint_score(n, edges, paths):
    """Exact rank-one updates; only branch-end deletions need be considered."""
    R = exact_resistance(n, edges)
    nominal = max(map(max, R))
    worst, witness = Fraction(0), None
    for ix in endpoint_indices(paths):
        u, v = edges[ix]
        denominator = 4 * (1 - R[u][v])
        assert denominator > 0
        for i in range(n):
            for j in range(i + 1, n):
                numerator = R[i][v] + R[j][u] - R[i][u] - R[j][v]
                value = R[i][j] + numerator * numerator / denominator
                if value > worst:
                    worst, witness = value, (ix, i, j)
    return worst, nominal, witness


def verify_catalogue():
    """Independent generation of all possible core edge multisets."""
    found = {}
    raw_count = 0
    for p in range(1, 5):
        pairs = list(it.combinations_with_replacement(range(p), 2))
        for edges in it.combinations_with_replacement(pairs, p + 2):
            raw_count += 1
            deg = [0] * p
            for u, v in edges:
                deg[u] += 1
                deg[v] += 1
            if min(deg) < 3 or not connected(p, edges):
                continue
            if any(not connected(p, edges[:i] + edges[i + 1:])
                   for i in range(len(edges))):
                continue
            found[(p, core_key(p, edges))] = True
    expected = {(p, core_key(p, edges)) for p, edges in CORES.values()}
    assert set(found) == expected
    return {'raw_core_multisets_checked': raw_count, 'nonisomorphic_cores': len(found)}


def verify(out):
    result = verify_catalogue()
    direct_checks = 0
    floating_checks = 0
    # Direct deletion/reinversion is independent of both the endpoint reduction
    # and the rank-one update used in the complete core search.
    for n in range(4, 9):
        for name in CORES:
            for lengths in length_vectors(name, n + 2):
                nn, edges, paths = subdivide(name, lengths)
                assert nn == n
                actual, _, _ = exact_endpoint_score(n, edges, paths)
                assert actual == score_exact(n, edges)
                direct_checks += 1
    # Longer chains probe endpoint reduction beyond the small exact cross-check.
    for n in (10, 15, 20):
        for name in CORES:
            for lengths in it.islice(length_vectors(name, n + 2), 3):
                nn, edges, paths = subdivide(name, lengths)
                actual, _, _ = exact_endpoint_score(nn, edges, paths)
                assert abs(float(actual) - score_fast(nn, edges)[0]) < 1e-8
                floating_checks += 1
    result.update(direct_rational_graph_checks=direct_checks,
                  longer_chain_floating_checks=floating_checks)
    labelled_checks = []
    for n in range(4, 8):
        m = n + 2
        labelled = sum(math.factorial(n) // automorphism_count(name, lengths)
                       for name in CORES for lengths in length_vectors(name, m))
        # Independent labelled edge-set enumeration: no stored files and no
        # core construction enter this connectivity-only count.
        direct = 0
        for candidate in it.combinations(it.combinations(range(n), 2), m):
            degrees = [0] * n
            for u, v in candidate:
                degrees[u] += 1; degrees[v] += 1
            if min(degrees) < 2:
                continue
            if all(connected(n, candidate[:ix] + candidate[ix+1:])
                   for ix in range(m)):
                direct += 1
        assert labelled == direct
        labelled_checks.append({'n': n, 'labelled_graphs': labelled})
    result['independent_labelled_count_checks'] = labelled_checks
    uniform_checks = []
    for length in range(1, 11):
        n, edges, paths = subdivide('K4', (length,) * 6)
        # Full floating deletion search, followed by exact small cases.
        assert abs(score_fast(n, edges)[0] - (2 * length - 1)) < 1e-8
        if length <= 3:
            assert score_exact(n, edges) == 2 * length - 1
        uniform_checks.append({'path_length': length, 'n': n, 'value': 2 * length - 1})
    result['uniform_K4_checks'] = uniform_checks
    (out / 'core_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
    return result


def search(nmin, nmax, out):
    records, witnesses, all_scores = [], [], []
    for n in range(nmin, nmax + 1):
        start = time.monotonic()
        best = None
        winners = []
        counts = Counter()
        nominal_best, nominal_worst_loss = None, None
        for name in CORES:
            for lengths in length_vectors(name, n + 2):
                nn, edges, paths = subdivide(name, lengths)
                assert nn == n and len(edges) == n + 2
                value, nominal, witness = exact_endpoint_score(n, edges, paths)
                counts[name] += 1
                all_scores.append(dict(n=n, core=name, lengths='-'.join(map(str, lengths)),
                                       robust_value=str(value), nominal_value=str(nominal)))
                item = dict(core=name, lengths=lengths, edges=edges,
                            nominal_value=str(nominal), worst_deletion_and_pair=witness)
                if best is None or value < best:
                    best, winners = value, [item]
                elif value == best:
                    winners.append(item)
                if nominal_best is None or nominal < nominal_best:
                    nominal_best, nominal_worst_loss = nominal, value
                elif nominal == nominal_best:
                    nominal_worst_loss = min(nominal_worst_loss, value)
        assert best is not None
        assert best >= Fraction(n - 1, 3)
        if n % 6 == 4:
            assert best == Fraction(n - 1, 3)
        # Direct inverse and all deletions certify a representative independently.
        assert score_exact(n, winners[0]['edges']) == best
        row = dict(n=n, m=n+2, nonisomorphic_candidates=sum(counts.values()),
                   minimum=str(best), optimizer_count=len(winners),
                   optimizer_cores=';'.join(sorted({x['core'] for x in winners})),
                   nominal_optimum_within_bridgeless_class=str(nominal_best),
                   best_robust_value_among_those_nominal_optimizers=str(nominal_worst_loss))
        records.append(row)
        witnesses.append(dict(n=n, minimum=str(best), counts_by_core=dict(counts),
                              optimizers=winners))
        # Timings are logged, not included in the deterministic result files.
        print(json.dumps({**row, 'seconds': round(time.monotonic()-start, 3)}), flush=True)
        write_results(records, witnesses, all_scores, out)
    return records


def write_results(records, witnesses, all_scores, out):
    for filename, rows in [('core_search.csv', records), ('core_scores.csv', all_scores)]:
        with (out / filename).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    (out / 'core_optimizers.json').write_text(json.dumps(witnesses, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nmin', type=int, default=4)
    parser.add_argument('--nmax', type=int, default=12)
    parser.add_argument('--mode', choices=['all', 'search', 'verify'], default='all')
    parser.add_argument('--out', type=Path, default=Path('results'))
    args = parser.parse_args()
    if args.nmin < 4 or args.nmax < args.nmin:
        parser.error('require 4 <= nmin <= nmax')
    args.out.mkdir(parents=True, exist_ok=True)
    if args.mode in ('all', 'verify'):
        verify(args.out)
    if args.mode in ('all', 'search'):
        search(args.nmin, args.nmax, args.out)
