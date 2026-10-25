"""Exact equal-budget benchmarks and a Gaussian calibration experiment.

All observations in the simulation are synthetic. The failure is chosen from
the graph before responses are generated. No response-dependent deletion is
simulated. Calibration consists of q independent zero-mean Gaussian residuals
with the same variance as a comparison, additional to the stated graph budget.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import itertools as it
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln
from scipy.stats import norm, t

from cores import subdivide
from research import (connected, exact_resistance, optimal_theta,
                      predicted_optimum, resistance, score_exact, score_fast, theta)


def spanning_trees(n, edges):
    # Independent combinatorial count, avoiding floating determinant comparisons.
    return sum(connected(n, list(s)) for s in it.combinations(edges, n - 1))


def benchmark_five():
    records = []
    for es in it.combinations(list(it.combinations(range(5), 2)), 7):
        edges = list(es)
        assert connected(5, edges)
        R = exact_resistance(5, edges)
        records.append(dict(edges=edges,
                            average=sum(R[i][j] for i in range(5)
                                        for j in range(i + 1, 5)) / 10,
                            maximum=max(map(max, R)),
                            trees=spanning_trees(5, edges),
                            robust=score_exact(5, edges)))
    results = []
    winner_sets = []
    for label, key, sign in [('A', 'average', 1), ('D', 'trees', -1),
                             ('MV', 'maximum', 1), ('one_loss_MV', 'robust', 1)]:
        best = min(sign * r[key] for r in records)
        winners = [r for r in records if sign * r[key] == best]
        winner_sets.append({tuple(r['edges']) for r in winners})
        results.append(dict(criterion=label, optimum=str(sign * best),
            labelled_optimizers=len(winners),
            average_range=sorted({str(r['average']) for r in winners}),
            nominal_max_range=sorted({str(r['maximum']) for r in winners}),
            spanning_tree_range=sorted({str(r['trees']) for r in winners}),
            post_loss_range=sorted({str(r['robust']) for r in winners}),
            representative_edges=winners[0]['edges']))
    assert len(records) == 120
    assert winner_sets[0] == winner_sets[1] == winner_sets[2]
    assert winner_sets[0].isdisjoint(winner_sets[3])
    return dict(n=5, m=7, candidates=120,
                connected_candidates=120,
                bridgeless_candidates=sum(math.isfinite(r['robust']) for r in records),
                criteria=results)


def scale_mean(df):
    return math.sqrt(2 / df) * math.exp(gammaln((df + 1) / 2) - gammaln(df / 2))


def cases():
    n = 28
    cycle = [(i, (i + 1) % n) for i in range(n)]
    nn, th = theta(*optimal_theta(n)[0])
    nk, k4, _ = subdivide('K4', (5,) * 6)
    assert nn == nk == n
    return [
        ('cycle_28', 28, cycle, Fraction(27)),
        ('chord_28', 28, cycle + [(0, 14)], Fraction(251, 15)),
        ('theta_28', 28, th, predicted_optimum(28)),
        ('K4_28', 28, k4, Fraction(9)),
        ('nominal_MV_5', 5, [(0,1),(0,3),(0,4),(1,2),(2,3),(2,4),(3,4)], Fraction(2)),
        ('robust_MV_5', 5, [(0,1),(0,2),(0,3),(0,4),(1,4),(2,4),(3,4)], Fraction(13,8)),
    ]


def simulate_case(name, n, edges, exact_phi, reps, calibration_df, seed):
    phi, witness, _, _ = score_fast(n, edges)
    assert abs(phi - float(exact_phi)) < 1e-8
    edge_index, wi, wj = witness
    kept = edges[:edge_index] + edges[edge_index + 1:]
    R, _ = resistance(n, kept)
    assert abs(R[wi, wj] - float(exact_phi)) < 1e-8
    X = np.zeros((len(kept), n - 1))
    for row, (u, v) in enumerate(kept):
        if u < n - 1: X[row, u] = 1
        if v < n - 1: X[row, v] = -1
    operator = np.linalg.solve(X.T @ X, X.T)
    # Independent numerical covariance identity from the observation map.
    C = np.zeros((n, n)); C[:-1, :-1] = operator @ operator.T
    Rd = C.diagonal()[:, None] + C.diagonal()[None, :] - 2 * C
    assert np.max(np.abs(Rd - R)) < 1e-8
    pairs = np.array(list(it.combinations(range(n), 2)))
    sr = np.sqrt(R[pairs[:, 0], pairs[:, 1]])
    alpha = .05; K = len(pairs)
    cq = t.isf(alpha / (2 * K), calibration_df)
    cz = norm.isf(alpha / (2 * K))
    marginal = t.isf(alpha / 2, calibration_df)
    rng = np.random.default_rng(seed)
    count_t = count_z = count_marginal = 0
    sum_w = sum_w2 = sum_scale = sum_scale2 = 0.
    # theta=0 and sigma=1 lose no generality for coverage under the linear model.
    for start in range(0, reps, 2500):
        b = min(2500, reps - start)
        noise = rng.standard_normal((b, len(kept)))
        estimates = np.zeros((b, n))
        estimates[:, :-1] = noise @ operator.T
        calibration = np.sqrt(rng.chisquare(calibration_df, b) / calibration_df)
        delta = estimates[:, pairs[:, 0]] - estimates[:, pairs[:, 1]]
        max_standardized = np.max(np.abs(delta) / sr, axis=1)
        count_t += int(np.sum(max_standardized <= cq * calibration))
        count_z += int(np.sum(max_standardized <= cz))
        worst = estimates[:, wi] - estimates[:, wj]
        count_marginal += int(np.sum(np.abs(worst) <= marginal * calibration * math.sqrt(phi)))
        sum_w += float(worst.sum()); sum_w2 += float(worst @ worst)
        sum_scale += float(calibration.sum()); sum_scale2 += float(calibration @ calibration)
    empirical_var = (sum_w2 - sum_w * sum_w / reps) / (reps - 1)
    p = count_t / reps
    scale_var = (sum_scale2 - sum_scale ** 2 / reps) / (reps - 1)
    residual_df = len(edges) - n
    residual_halfwidth = (t.isf(alpha / (2 * K), residual_df) *
                          scale_mean(residual_df) * math.sqrt(phi)) if residual_df else None
    return dict(design=name, n=n, graph_comparisons=len(edges),
        calibration_df=calibration_df, total_measurements=len(edges)+calibration_df,
        replications=reps, seed=int(seed), deleted_edge=list(edges[edge_index]),
        worst_pair=[wi,wj], exact_post_loss_variance=str(exact_phi),
        empirical_worst_variance=empirical_var,
        variance_mcse=phi * math.sqrt(2 / (reps - 1)),
        simultaneous_t_coverage=p, coverage_mcse=math.sqrt(p * (1-p) / reps),
        simultaneous_known_sigma_coverage=count_z/reps,
        marginal_t_coverage=count_marginal/reps,
        expected_max_halfwidth=cq * scale_mean(calibration_df) * math.sqrt(phi),
        empirical_mean_max_halfwidth=cq * math.sqrt(phi) * sum_scale/reps,
        mean_halfwidth_mcse=cq * math.sqrt(phi) * math.sqrt(scale_var/reps),
        known_sigma_max_halfwidth=cz * math.sqrt(phi),
        residual_only_df=residual_df,
        residual_only_expected_max_halfwidth=residual_halfwidth)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, default=Path('results'))
    p.add_argument('--reps', type=int, default=50000)
    p.add_argument('--calibration-df', type=int, default=20)
    p.add_argument('--seed', type=int, default=20261025)
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    if args.reps < 2 or args.calibration_df < 1:
        p.error('reps must be at least two; calibration-df must be positive')
    benchmark = benchmark_five()
    seeds = np.random.SeedSequence(args.seed).generate_state(len(cases()))
    rows = [simulate_case(*case, args.reps, args.calibration_df, seed)
            for case, seed in zip(cases(), seeds)]
    if args.verify:
        # Distributional smoke checks: broad thresholds, fixed data-generating model.
        for row in rows:
            assert abs(row['empirical_worst_variance'] - float(Fraction(row['exact_post_loss_variance']))) < 6 * row['variance_mcse']
            assert abs(row['empirical_mean_max_halfwidth'] - row['expected_max_halfwidth']) < 6 * row['mean_halfwidth_mcse']
            assert row['simultaneous_t_coverage'] > .94
            assert abs(row['marginal_t_coverage'] - .95) < .025
        print('PASS: exact benchmark, covariance identity and Gaussian calibration checks')
    else:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out/'equal_budget_benchmark.json').write_text(json.dumps(benchmark, indent=2)+'\n')
        (args.out/'calibration_study.json').write_text(json.dumps(rows, indent=2)+'\n')
        with (args.out/'calibration_study.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        print(json.dumps(benchmark, indent=2))
        for row in rows:
            print(row['design'], 'coverage=', round(row['simultaneous_t_coverage'], 5),
                  'expected half-width=', round(row['expected_max_halfwidth'], 4))


if __name__ == '__main__':
    main()
