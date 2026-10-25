"""Independent computational checks for one-observation-loss comparison designs.

Unit conductances; simple undirected graphs; each edge is one comparison.
Python >= 3.10 and numpy. No graph-library dependency.
"""
from __future__ import annotations
import argparse
import csv
import itertools as it
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np


def cycle_radius(s: int) -> Fraction:
    return Fraction(s*s//4, s)


def theta(a: int, b: int, c: int):
    edges = []
    nxt = 2
    for length in (a, b, c):
        path = [0] + list(range(nxt, nxt+length-1)) + [1]
        nxt += length-1
        edges.extend(zip(path[:-1], path[1:]))
    return nxt, edges


def figure_eight(a: int, b: int):
    edges = []
    nxt = 1
    for length in (a, b):
        path = [0] + list(range(nxt, nxt+length-1)) + [0]
        nxt += length-1
        edges.extend(zip(path[:-1], path[1:]))
    return nxt, edges


def connected(n: int, edges) -> bool:
    adj = [[] for _ in range(n)]
    for u,v in edges:
        adj[u].append(v); adj[v].append(u)
    seen = {0}; stack = [0]
    while stack:
        for v in adj[stack.pop()]:
            if v not in seen:
                seen.add(v); stack.append(v)
    return len(seen)==n


def resistance(n: int, edges):
    L = np.zeros((n,n))
    for u,v in edges:
        L[u,u] += 1; L[v,v] += 1
        L[u,v] -= 1; L[v,u] -= 1
    K = np.zeros((n,n))
    K[:-1,:-1] = np.linalg.inv(L[:-1,:-1])
    d = K.diagonal()
    return d[:,None]+d[None,:]-2*K, K


def score_fast(n: int, edges):
    R,K = resistance(n,edges)
    worst = -1.; witness = None
    for ix,(u,v) in enumerate(edges):
        denom = 1.-R[u,v]
        if denom < 1e-9:
            return math.inf, None, float(R.max()), float(R.sum()/(n*(n-1)))
        z = K[:,u]-K[:,v]
        updated = R+(z[:,None]-z[None,:])**2/denom
        ij = np.unravel_index(updated.argmax(),updated.shape)
        value = float(updated[ij])
        if value>worst:
            worst=value; witness=(ix,int(ij[0]),int(ij[1]))
    return worst,witness,float(R.max()),float(R.sum()/(n*(n-1)))


def score_direct(n: int, edges):
    ans=0.
    for ix in range(len(edges)):
        kept = edges[:ix]+edges[ix+1:]
        if not connected(n,kept): return math.inf
        R,_ = resistance(n,kept)
        ans=max(ans,float(R.max()))
    return ans


def exact_resistance(n: int, edges):
    """Fraction Gauss-Jordan inverse of a grounded Laplacian."""
    p=n-1
    A=[[Fraction(0) for _ in range(2*p)] for _ in range(p)]
    for i in range(p): A[i][p+i]=Fraction(1)
    for u,v in edges:
        if u<p: A[u][u]+=1
        if v<p: A[v][v]+=1
        if u<p and v<p: A[u][v]-=1; A[v][u]-=1
    for j in range(p):
        pivot=next(i for i in range(j,p) if A[i][j])
        A[j],A[pivot]=A[pivot],A[j]
        c=A[j][j]; A[j]=[x/c for x in A[j]]
        for i in range(p):
            if i!=j:
                c=A[i][j]
                A[i]=[x-c*y for x,y in zip(A[i],A[j])]
    K=[[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(p):
        K[i][:p]=A[i][p:]
    return [[K[i][i]+K[j][j]-2*K[i][j] for j in range(n)] for i in range(n)]


def score_exact(n,edges):
    values=[]
    for ix in range(len(edges)):
        kept=edges[:ix]+edges[ix+1:]
        if not connected(n,kept): return math.inf
        values.append(max(map(max,exact_resistance(n,kept))))
    return max(values)


def optimal_theta(n):
    S=n+1; c=(S+2)//3
    return [(a,S-c-a,c) for a in range(1,c+1)
            if a<=S-c-a<=c and S-c-a>=2]


def predicted_optimum(n):
    c=(n+3)//3
    return Fraction(c-1)+cycle_radius(n+1-c)


def enumerations(nmax, out):
    records=[]; representatives=[]
    for n in range(4,nmax+1):
        possible=list(it.combinations(range(n),2))
        for m in range(n,min(n+2,len(possible))+1):
            best=math.inf; best_count=0; robust_count=0; min_degree_count=0
            nominal_best=math.inf; nominal_robust=math.inf
            for edges in it.combinations(possible,m):
                degree=[0]*n
                for u,v in edges: degree[u]+=1; degree[v]+=1
                if min(degree)<2: continue
                min_degree_count+=1
                if not connected(n,edges): continue
                value,witness,nominal,_=score_fast(n,edges)
                if not math.isfinite(value): continue
                robust_count+=1
                if n<=5:
                    assert abs(value-score_direct(n,list(edges)))<1e-8
                if nominal<nominal_best-1e-8:
                    nominal_best=nominal; nominal_robust=value
                elif abs(nominal-nominal_best)<1e-8:
                    nominal_robust=min(nominal_robust,value)
                if value<best-1e-8:
                    best=value; best_count=1; representative=list(edges)
                    best_nominal=nominal
                elif abs(value-best)<1e-8:
                    best_count+=1
                    if nominal<best_nominal:
                        best_nominal=nominal; representative=list(edges)
            predicted=predicted_optimum(n) if m==n+1 else (Fraction(n-1) if m==n else None)
            if predicted is not None: assert abs(best-float(predicted))<1e-8
            exact=score_exact(n,representative)
            assert abs(best-float(exact))<1e-8
            row=dict(n=n,m=m,labelled_graphs=math.comb(len(possible),m),
                min_degree_two_graphs=min_degree_count,robust_labelled_graphs=robust_count,
                optimal_labelled_graphs=best_count,minimax_value=str(exact),
                minimax_decimal=round(best,12),nominal_max_variance=round(best_nominal,12),
                nominal_optimal_robust_value=round(nominal_robust,12))
            print(json.dumps(row),flush=True); records.append(row)
            representatives.append(dict(n=n,m=m,edges=representative,value=str(exact)))
    with (out/'enumeration.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=records[0].keys());w.writeheader();w.writerows(records)
    (out/'representatives.json').write_text(json.dumps(representatives,indent=2))
    return records


def family_checks(out):
    counts=dict(theta_direct_checks=0,figure_eight_direct_checks=0,
                rational_checks=0,integer_minimization_checks=0)
    for n in range(4,26):
        S=n+1
        for a in range(1,S):
            for b in range(max(a,2),S):
                c=S-a-b
                if c<b: continue
                nn,edges=theta(a,b,c)
                value=Fraction(c-1)+cycle_radius(a+b)
                assert nn==n and len(edges)==n+1
                assert abs(score_direct(n,edges)-float(value))<1e-8
                counts['theta_direct_checks']+=1
                if n<=9:
                    assert score_exact(n,edges)==value
                    counts['rational_checks']+=1
        for a in range(3,S//2+1):
            b=S-a
            nn,edges=figure_eight(a,b)
            value=Fraction(b-1)+cycle_radius(a)
            assert nn==n
            assert abs(score_direct(n,edges)-float(value))<1e-8
            counts['figure_eight_direct_checks']+=1
            if n<=9:
                assert score_exact(n,edges)==value
                counts['rational_checks']+=1
    for n in range(4,201):
        S=n+1; vals=[]
        for a in range(1,S):
            for b in range(max(a,2),S):
                c=S-a-b
                if c>=b: vals.append((Fraction(c-1)+cycle_radius(a+b),(a,b,c)))
        low=min(v for v,_ in vals)
        mins=sorted(abc for v,abc in vals if v==low)
        assert low==predicted_optimum(n)
        assert mins==optimal_theta(n)
        for a in range(3,S//2+1):
            assert Fraction(S-a-1)+cycle_radius(a)>low
        counts['integer_minimization_checks']+=1
    (out/'verification_counts.json').write_text(json.dumps(counts,indent=2))
    print(json.dumps(counts),flush=True)


def examples(out):
    rows=[]
    for n in [6,8,10,12,20,30,50,100]:
        a=n//2; b=n-a
        chord=Fraction(b-1)+cycle_radius(a+1)
        optimum=predicted_optimum(n)
        rows.append(dict(n=n,cycle_budget=n,augmented_budget=n+1,
                         cycle_robust_variance=str(n-1),
                         antipodal_chord_variance=str(chord),
                         optimal_redesign_variance=str(optimum),
                         optimal_theta_paths=str(optimal_theta(n)),
                         chord_reduction_pct=100*(1-float(chord)/(n-1)),
                         redesign_reduction_pct=100*(1-float(optimum)/(n-1))))
    with (out/'examples.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    print(json.dumps(rows),flush=True)


def exact_tradeoff(out):
    """Complete rational check: 5 objects, 7 comparisons (120 labelled graphs)."""
    n=5; m=7; rows=[]
    for edges in it.combinations(list(it.combinations(range(n),2)),m):
        edges=list(edges)
        assert connected(n,edges)
        nominal=max(map(max,exact_resistance(n,edges)))
        robust=score_exact(n,edges)
        rows.append((nominal,robust,edges))
    nominal_best=min(x[0] for x in rows)
    robust_best=min(x[1] for x in rows)
    nominal_choices=[x for x in rows if x[0]==nominal_best]
    robust_choices=[x for x in rows if x[1]==robust_best]
    assert nominal_best==Fraction(7,8)
    assert robust_best==Fraction(13,8)
    assert set(x[1] for x in nominal_choices)=={Fraction(2)}
    assert set(x[0] for x in robust_choices)=={Fraction(1)}
    result=dict(labelled_graphs=len(rows),nominal_optimum=str(nominal_best),
                robust_optimum=str(robust_best),
                nominal_optimizer=dict(edges=nominal_choices[0][2],
                                       nominal='7/8',robust='2'),
                robust_optimizer=dict(edges=robust_choices[0][2],
                                      nominal='1',robust='13/8'),
                nominal_optimizer_count=len(nominal_choices),
                robust_optimizer_count=len(robust_choices))
    (out/'exact_tradeoff.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--nmax',type=int,default=7)
    parser.add_argument('--out',type=Path,default=Path('results'))
    parser.add_argument('--mode',choices=['all','enumerate','families','examples','tradeoff'],default='all')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if args.mode in ['all','enumerate']: enumerations(args.nmax,args.out)
    if args.mode in ['all','families']: family_checks(args.out)
    if args.mode in ['all','examples']: examples(args.out)
    if args.mode in ['all','tradeoff']: exact_tradeoff(args.out)
