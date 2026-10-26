"""Exact subdivision certificates and comparisons of one-loss design criteria.

All reported optima in the n+1 study are over simple, connected, bridgeless
graphs. The theta/figure-eight classification makes the enumeration complete.
No empirical observations or saved output files are needed by --verify.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import itertools as it
import json
import math
from pathlib import Path

from research import (theta, figure_eight, exact_resistance, score_direct,
                      score_exact, cycle_radius, predicted_optimum)

ROOT = Path(__file__).resolve().parent


def symmetric_cores():
    pairs = list(it.combinations(range(5), 2))
    return {
        'K4': (4, list(it.combinations(range(4), 2))),
        'K33': (6, [(u, v) for u in range(3) for v in range(3, 6)]),
        'Q3': (8, [(u, u ^ (1 << k)) for u in range(8)
                   for k in range(3) if u < (u ^ (1 << k))]),
        'Petersen': (10, [(i, j) for i in range(10) for j in range(i+1, 10)
                         if not set(pairs[i]) & set(pairs[j])]),
    }


def subdivision(p, es, lengths):
    edges, nxt = [], p
    for (u, v), length in zip(es, lengths):
        path = [u] + list(range(nxt, nxt + length - 1)) + [v]
        nxt += length - 1
        edges.extend(zip(path[:-1], path[1:]))
    return nxt, edges


def quadratic_max(A, B, c):
    """Maximum of A+(B-A+c)t-c*t*t on [0,1], exactly."""
    assert c >= 0
    vals = [(A, F(0)), (B, F(1))]
    if c:
        t = (B-A+c)/(2*c)
        if 0 < t < 1:
            vals.append((A+(B-A+c)**2/(4*c), t))
    return max(vals)


def check_inverse_certificate(R, es):
    """Check grounded L*K=I from the displayed resistance matrix."""
    p = len(R)
    K = [[(R[i][0]+R[j][0]-R[i][j])/2 for j in range(1,p)]
         for i in range(1,p)]
    L = [[0]*(p-1) for _ in range(p-1)]
    for u,v in es:
        if u: L[u-1][u-1] += 1
        if v: L[v-1][v-1] += 1
        if u and v: L[u-1][v-1] -= 1; L[v-1][u-1] -= 1
    assert all(sum(L[i][k]*K[k][j] for k in range(p-1)) == int(i==j)
               for i in range(p-1) for j in range(p-1))


def certificates():
    records = []
    for name,(p,es) in symmetric_cores().items():
        q=len(es); r=q-p+1; alpha=F(p-1,r)
        R0=exact_resistance(p,es)
        assert {R0[u][v] for u,v in es} == {F(p-1,q)}
        counts=0
        for ix,e in enumerate(es):
            fs=es[:ix]+es[ix+1:]; R=exact_resistance(p,fs)
            check_inverse_certificate(R,fs)
            assert max(map(max,R)) == alpha == R[e[0]][e[1]]
            quadratics=[]
            for w in e:
                for a,b in fs:
                    A,B,c=R[w][a],R[w][b],1-R[a][b]
                    maximum,t=quadratic_max(A,B,c)
                    assert maximum <= alpha
                    quadratics.append(dict(root=w,edge=[a,b],A=str(A),B=str(B),
                        c=str(c),maximum=str(maximum),argmax=str(t)))
            counts += len(quadratics)
            if ix==0:
                denominator=math.lcm(*(x.denominator for row in R for x in row))
                record=dict(core=name,p=p,q=q,r=r,edge=list(e),alpha=str(alpha),
                    diameter=str(max(map(max,R))),
                    endpoint_continuous_eccentricity=str(max(F(x['maximum']) for x in quadratics)),
                    resistance_denominator=denominator,
                    integer_resistance_matrix=[[int(x*denominator) for x in row] for row in R],
                    edges=es,quadratics=quadratics)
        record['deleted_edges_checked']=q
        record['endpoint_quadratics_checked']=counts
        records.append(record)
    return records


def exact_metrics(n, edges):
    R=exact_resistance(n,edges); pairs=list(it.combinations(range(n),2)); K=len(pairs)
    nominal_A=sum(R[i][j] for i,j in pairs)/K
    diameters=[]; averages=[]
    for u,v in edges:
        denom=4*(1-R[u][v]); assert denom>0
        vals=[]
        for i,j in pairs:
            d=R[i][v]+R[j][u]-R[i][u]-R[j][v]
            vals.append(R[i][j]+d*d/denom)
        diameters.append(max(vals)); averages.append(sum(vals)/K)
    return dict(M=max(map(max,R)),A=nominal_A,Phi=max(diameters),
                expected_MV=sum(diameters)/len(edges),
                expected_A=sum(averages)/len(edges),worst_A=max(averages),
                loss_diameters=diameters)


def designs(n):
    S=n+1
    for a in range(1,S):
        for b in range(max(a,2),S):
            c=S-a-b
            if c<b: continue
            nn,es=theta(a,b,c); assert nn==n
            yield 'theta',(a,b,c),es,a*b+a*c+b*c
    for a in range(3,S//2+1):
        b=S-a
        nn,es=figure_eight(a,b); assert nn==n
        yield 'figure_eight',(a,b),es,a*b


def comparison_study(nmax):
    records=[]; summaries=[]
    criteria=['Phi','M','A','T','expected_MV','expected_A','worst_A']
    for n in range(4,nmax+1):
        candidates=[]
        for kind,lengths,edges,trees in designs(n):
            metrics=exact_metrics(n,edges); metrics['T']=F(trees)
            predicted=lengths[-1]-1+cycle_radius(sum(lengths[:-1]))
            assert metrics['Phi']==predicted
            candidates.append(dict(n=n,kind=kind,lengths=lengths,edges=edges,**metrics))
        opt={k:(max if k=='T' else min)(x[k] for x in candidates) for k in criteria}
        assert opt['Phi']==predicted_optimum(n)
        winners={k:[x for x in candidates if x[k]==opt[k]] for k in criteria}
        rob=winners['Phi']
        summary=dict(n=n,candidates=len(candidates),optima=opt,
            winners={k:[dict(kind=x['kind'],lengths=x['lengths']) for x in ws]
                     for k,ws in winners.items()},
            robust_efficiencies={k:[min(opt[k]/x[k] for x in rob),max(opt[k]/x[k] for x in rob)]
                                 for k in ['A','M','expected_MV','expected_A','worst_A']},
            robust_D_efficiency=[min(float(x['T']/opt['T'])**(1/(n-1)) for x in rob),
                                 max(float(x['T']/opt['T'])**(1/(n-1)) for x in rob)],
            robust_penalty_of_nominal={k:[min(x['Phi']/opt['Phi'] for x in winners[k]),
                                          max(x['Phi']/opt['Phi'] for x in winners[k])]
                                      for k in ['A','M','T','expected_MV','expected_A','worst_A']})
        constrained=[x for x in candidates if x['Phi']<=F(11,10)*opt['Phi']]
        copt=min(x['A'] for x in constrained)
        summary['A_opt_under_10pct_Phi_relaxation']=dict(A=copt,
            A_improvement_over_best_robust=1-copt/min(x['A'] for x in rob),
            winners=[dict(kind=x['kind'],lengths=x['lengths'],Phi=x['Phi']) for x in constrained if x['A']==copt])
        records.extend(candidates); summaries.append(summary)
        print(f'n={n}: {len(candidates)} exact candidates; Phi={opt["Phi"]}',flush=True)
    return dict(scope='all simple connected bridgeless n-vertex graphs with m=n+1',
                loss_model='one uniformly chosen edge for expected criteria',
                summaries=summaries,candidates=records)


def verify():
    cert=certificates(); direct=0; interpolated=0
    for name,(p,es) in symmetric_cores().items():
        r=len(es)-p+1
        # Independent check of the continuous interpolation identity at
        # all third-grid points, using a directly inverted subdivided graph.
        fs=es[1:]; R=exact_resistance(p,fs)
        nn,ee=subdivision(p,fs,[3]*len(fs)); RR=exact_resistance(nn,ee)
        for j,(a,b) in enumerate(fs):
            for k in (1,2):
                t=F(k,3); x=p+2*j+k-1
                for w in es[0]:
                    value=(1-t)*R[w][a]+t*R[w][b]+t*(1-t)*(1-R[a][b])
                    assert RR[w][x]==3*value
                    interpolated+=1
        for ell in range(1,5):
            n,edges=subdivision(p,es,[ell]*len(es)); bound=F(n-1,r)
            assert abs(score_direct(n,edges)-float(bound))<1e-8
            if ell==1 or (ell==2 and name in ['K4','K33']):
                assert score_exact(n,edges)==bound
            direct+=1
    for n in range(4,9):
        for kind,lengths,edges,trees in designs(n):
            met=exact_metrics(n,edges)
            assert met['Phi']==score_exact(n,edges)
            # Independent deletion/reinversion checks every expected-loss term.
            direct_averages=[]
            for ix in range(len(edges)):
                R=exact_resistance(n,edges[:ix]+edges[ix+1:])
                assert max(map(max,R))==met['loss_diameters'][ix]
                direct_averages.append(sum(R[i][j] for i in range(n) for j in range(i+1,n))/F(n*(n-1),2))
            assert sum(direct_averages)/len(edges)==met['expected_A']
            assert max(direct_averages)==met['worst_A']
    print(json.dumps(dict(core_certificates=len(cert),
        exact_endpoint_quadratics=sum(x['endpoint_quadratics_checked'] for x in cert),
        exact_interpolation_checks=interpolated,
        full_numerical_subdivisions=direct,verification='PASS')))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--verify',action='store_true')
    ap.add_argument('--nmax',type=int,default=20); ap.add_argument('--certificates-only',action='store_true')
    ap.add_argument('--out',type=Path,default=ROOT/'results'); args=ap.parse_args()
    if args.verify: verify(); return
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'subdivision_certificates.json').write_text(json.dumps(certificates(),indent=2))
    if not args.certificates_only:
        data=comparison_study(args.nmax)
        (args.out/'criterion_comparisons.json').write_text(json.dumps(data,indent=2,default=str))


if __name__=='__main__': main()
