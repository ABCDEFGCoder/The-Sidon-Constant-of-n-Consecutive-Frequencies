"""Finite sign quotients of the symmetric QSP critical charts.

Odd degree: simultaneous conjugation gives a two-sheeted quotient.
Even degree: variable rotation by pi gives a second independent sign.
Every use records nonzero pivots; no zero-pivot stratum is discarded
from an assertion about the original global problem.
"""
from functools import lru_cache
from pathlib import Path
import argparse, json
import sympy as s
from qsp_chart import chart, critical


def primitive(f, variables):
    p = s.Poly(f, *variables).clear_denoms()[1].primitive()[1]
    if p.LC() < 0:
        p = -p
    return p.as_expr()


@lru_cache(None)
def model(n):
    t, P, Q, E = chart(n)
    m = len(t)
    if n % 2 or m == 1:
        names = ['B'] + (['C', 'D'] if m == 3 else [f'C{j}' for j in range(2, m+1)])
        blocks = [list(range(m))]
    else:
        odd = list(range(0, m, 2)); even = list(range(1, m, 2))
        names = ['B'] + (['C'] if m == 4 else [f'C{j+1}' for j in odd[1:]])
        names += ['U'] + (['V'] if m == 4 else [f'C{j+1}' for j in even[1:]])
        blocks = [odd, even]
    q = list(s.symbols(' '.join(names), seq=True))
    coordinate = {}; images = []; pivots = []; offset = 0
    for block in blocks:
        pivot = q[offset]; pivots.append(pivot)
        coordinate[block[0]] = pivot
        images.append(t[block[0]]**2)
        for j, var in zip(block[1:], q[offset+1:offset+len(block)]):
            coordinate[j] = var
            images.append(t[block[0]]*t[j])
        offset += len(block)

    def quotient(f):
        out = s.Integer(0)
        for exponents, c in s.Poly(s.expand(f), *t).terms():
            term = c
            for block in blocks:
                a = exponents[block[0]] - sum(exponents[j] for j in block[1:])
                assert a % 2 == 0, ('Noninvariant monomial', n, exponents)
                term *= coordinate[block[0]]**(a//2)
                for j in block[1:]:
                    term *= coordinate[j]**exponents[j]
            out += term
        return s.cancel(out)

    R = [quotient(s.expand(P[j]*s.conjugate(P[j]))) for j in range(1, m+1)]
    factored_R = []
    for j in range(1, m+1):
        pieces = []
        for component in [s.re(P[j]), s.im(P[j])]:
            component = s.expand(component)
            if component == 0: continue
            exponents = s.Poly(component, *t).monoms()
            character = [sum(exponents[0][k] for k in block) % 2 for block in blocks]
            assert all([sum(e[k] for k in block) % 2 for block in blocks] == character for e in exponents)
            multiplier = s.prod(t[block[0]]**bit for block, bit in zip(blocks, character))
            denominator = s.prod(pivot**bit for pivot, bit in zip(pivots, character))
            pieces.append(quotient(component*multiplier)**2/denominator)
        factored_R.append(s.Add(*pieces, evaluate=False))
    energy = quotient(E)
    r = list(s.symbols(f'r1:{m+1}')); variables = q+r
    weights = [2]*m
    if n % 2 == 0:
        weights[-1] = 1
    expressions = [rj*rj-v for rj, v in zip(r, R)]
    gradients = []
    for x in q:
        logder = s.cancel(s.diff(energy, x)/energy)
        f = sum(w*(s.diff(v, x)-v*logder)*s.prod(r[k] for k in range(m) if k != i)
                for i, (w, v) in enumerate(zip(weights, R))) - 2*logder*s.prod(r)
        gradients.append(s.cancel(f))
    expressions += gradients
    equations = [primitive(s.fraction(s.cancel(f))[0], variables) for f in expressions]
    saturation = s.prod(pivots)*s.prod(r)
    for x in t:
        saturation *= s.fraction(quotient(1+x*x))[0]
    saturation = s.factor(saturation)
    L = 2+sum(w*x for w, x in zip(weights, r))
    return {'n': n, 't': list(t), 'variables': variables, 'q': q, 'r': r,
            'images': images, 'R': R, 'R_factored': factored_R, 'energy': energy, 'gradients': gradients,
            'equations': equations, 'saturation': saturation, 'L': L,
            'pivots': pivots, 'sheets': 2**len(blocks)}


def audit(n, chain_rule=True):
    d = model(n); t = d['t']; q = d['q']; images = d['images']
    sub = dict(zip(q, images)); _, P, _, E = chart(n)
    assert s.cancel(d['energy'].subs(sub)-E) == 0
    for j, v in enumerate(d['R'], 1):
        assert s.cancel(v.subs(sub)-s.expand(P[j]*s.conjugate(P[j]))) == 0
    jac = s.Matrix(images).jacobian(t)
    determinant = s.factor(jac.det())
    if chain_rule:
        vs, fs, *_ = critical(n)
        for j, x in enumerate(t):
            lhs = sum(g.subs(sub)*s.diff(y, x) for g, y in zip(d['gradients'], images))
            assert s.cancel(lhs-fs[len(t)+j]/(1+x*x)) == 0
    return {'n': n, 'sheets': d['sheets'], 'jacobian_determinant': str(determinant),
            'degrees': [int(s.total_degree(f)) for f in d['equations']],
            'terms': [len(s.Poly(f, *d['variables']).terms()) for f in d['equations']]}


