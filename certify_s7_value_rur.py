"""Accept a proposed S7 target RUR only after exact identities and bridging.

The full global cover is run by verify_all.py. This module rechecks
the exact RUR identities, minimality, and the local global-to-value bridge.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time
import sympy as s
from flint import fmpq, fmpq_poly, fmpz_poly, nmod_poly, arb, arb_poly, ctx
from exact_values import load_rur, primitive, simple_root_taylor
from irreducibility import prime, irreducible_mod, verify_sieve
from qsp_symmetry_quotient import model, audit
from qsp_critical_certificate import verify as verify_old_bridge
from integer_circuit_verify import verify as verify_integer_identities

sys.set_int_max_str_digits(0)
INTERVAL = [[2551658061766431151, 10**18], [2551658061766431152, 10**18]]


def irreducibility_witness(h, limit=16):
    D = h.degree(); possible = set(range(D+1)); rows = []; p = 1000003
    for _ in range(limit):
        p = int(s.nextprime(p))
        f = nmod_poly(h.coeffs(), p)
        if f.degree() != D or f.gcd(f.derivative()).degree() != 0: continue
        unit, factors = f.factor(); product = nmod_poly([int(unit)], p); subsets = {0}; entries = []
        for g, e in factors:
            assert e == 1; irreducible_mod(g, p); product *= g
            subsets |= {a+g.degree() for a in list(subsets)}
            entries.append([list(map(int, g.coeffs())), int(e)])
        assert product == f
        rows.append({'prime': p, 'unit': int(unit), 'factors': entries}); possible &= subsets
        print('MINPOLY_DEGREE_SIEVE', p, 'remaining', sorted(possible), flush=True)
        if possible == {0, D}:
            c = {'format': 'modular_degree_intersection_v1', 'degree': D, 'rows': rows}
            verify_sieve(h, c); return c
    raise ValueError('The modular degree sieve has not proved this polynomial irreducible')


def shifted(poly, centre):
    return arb_poly([arb(v) for v in poly.coeffs()])(arb_poly([arb(centre), arb(1)]))


def real_branch(raw):
    old = json.loads(Path('S7_QSP_critical_root.json').read_text())
    inherited = verify_old_bridge(old, bridge=True)
    assert inherited['global_bridge_checked']
    H = fmpz_poly(raw['parameter_polynomial']); a, b = [fmpq(*v) for v in INTERVAL]
    height = max(abs(int(v)).bit_length() for v in H.coeffs())
    for coefficients, scale in raw['numerators'].values():
        height = max(height, int(scale).bit_length(), max(abs(int(v)).bit_length() for v in coefficients))
    previous = ctx.prec; ctx.prec = height+4*(H.degree()+1)+2048
    try:
        centre = (a+b)/2; radius = (b-a)/2
        P = shifted(H, centre); derivative = P.derivative(); X = arb(0, radius)
        assert not derivative(X).contains(0)
        for _ in range(24):
            lower, upper = X.lower().fmpq(), X.upper().fmpq(); midpoint = (lower+upper)/2
            newton = arb(midpoint)-P(arb(midpoint))/derivative(X)
            lower = max(lower, newton.lower().fmpq()); upper = min(upper, newton.upper().fmpq())
            assert lower < upper
            X = arb((lower+upper)/2, (upper-lower)/2)
            if upper-lower < fmpq(1,10**80): break
        else: raise AssertionError('Interval Newton did not reach the requested branch precision')
        den = shifted(fmpz_poly(raw['common_denominator']), centre)(X); assert not den.contains(0)
        values = {name: shifted(fmpz_poly(coefficients), centre)(X)/(den*arb(scale))
                  for name, (coefficients, scale) in raw['numerators'].items()}
        B = values['B']; assert B.lower() > 0
        t1 = -B.sqrt(); point = [t1, values['C']/t1, values['D']/t1]+[values[f'r{j}'] for j in range(1,4)]
        assert all(v.lower() > 0 for v in point[3:])
        large_radius = fmpq(old['bridge_radius'])
        for value, c in zip(point, old['centre']):
            c = fmpq(*c)
            assert c-large_radius < value.lower() and value.upper() < c+large_radius, 'RUR point is outside the proved global QSP uniqueness box'
        return {'accepted': True, 'precision_bits': ctx.prec,
                'fine_root_interval': [str(centre+X.lower().fmpq()), str(centre+X.upper().fmpq())],
                'QSP_box': 'S7_QSP_critical_root.json', 'global_source': 'Sidon6_7_Exact_Closure_Certificate.py complete complex-coefficient cover',
                'full_global_cover_rerun_this_round': False}
    finally: ctx.prec = previous


def certify(path, write=True):
    start = time.monotonic(); raw = load_rur(path); d = model(7); variables = d['variables']+[s.Symbol('value')]
    assert raw['characteristic'] == 0 and raw['names'] == list(map(str, variables))
    H = fmpz_poly(raw['parameter_polynomial']); assert primitive(H) == H and H.degree() == raw['degree'] == 855
    assert fmpz_poly(raw['common_denominator']) == H.derivative()
    a, b = [fmpq(*v) for v in INTERVAL]; uniqueness = simple_root_taylor(H, a, b)
    try:
        minimal = H; irreducibility = irreducibility_witness(H)
    except ValueError:
        unit, factors = H.factor(); candidates = [f for f, e in factors if e == 1 and fmpq_poly(f)(a)*fmpq_poly(f)(b) < 0]
        assert len(candidates) == 1, 'No unique rational factor was selected by the proved root interval'
        minimal = primitive(candidates[0]); irreducibility = irreducibility_witness(minimal)
    B, C, D = d['q']; r = d['r']; target = variables[-1]*(1+B)*(B+C*C)*(B+D*D)-B*B*(1+sum(r))
    print('EXACT_RUR_IDENTITIES_START', H.degree(), flush=True)
    verify_integer_identities(raw, variables, [target]+d['equations'])
    print('EXACT_RUR_IDENTITIES_COMPLETE', flush=True)
    audit(7); bridge = real_branch(raw)
    exact = dict(raw); exact['status'] = 'EXACT_INTEGER_IDENTITIES_AND_GLOBAL_REAL_BRANCH_VERIFIED'
    exact_text = json.dumps(exact, separators=(',', ':'))
    proof = {'format': 'Sn_minimal_polynomial_and_global_bridge_v1', 'n': 7,
             'coefficient_order': 'constant_upwards', 'polynomial_coefficients': list(map(int, minimal.coeffs())),
             'isolating_interval': INTERVAL, 'root_uniqueness': uniqueness,
             'irreducibility': irreducibility, 'rur_source_file': 'S7_value_RUR_EXACT.json',
             'rur_source_sha256': hashlib.sha256(exact_text.encode()).hexdigest(),
             'global_bridge': bridge, 'seconds': time.monotonic()-start}
    if write:
        Path('S7_minimal_polynomial_and_interval.json').write_text(json.dumps(proof, separators=(',', ':')))
        Path('S7_value_RUR_EXACT.json').write_text(exact_text)
    print('S7_FULL_CLOSURE_VERIFIED', minimal.degree(), 'seconds', proof['seconds'], flush=True)
    return proof


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('rur', type=Path, nargs='?', default=Path('S7_value_RUR_EXACT.json'))
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args(); certify(args.rur, write=not args.verify_only)
