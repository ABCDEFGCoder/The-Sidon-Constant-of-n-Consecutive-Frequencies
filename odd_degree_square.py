"""An explicit squared-value minimal polynomial from odd field degree.

If h(alpha)=0 is irreducible of odd degree d, then
Q(alpha)=Q(alpha^2). Writing h(X)=A(X^2)+X B(X^2), the primitive
polynomial Y B(Y)^2-A(Y)^2 is the degree-d minimal polynomial of alpha^2.
"""
from pathlib import Path
import argparse, hashlib, json
from flint import fmpz_poly, fmpq
from exact_values import primitive, simple_root_taylor
from irreducibility import verify_sieve


def produce(source, output=None):
    source = Path(source); c = json.loads(source.read_text()); h = fmpz_poly(c['polynomial_coefficients'])
    d = h.degree(); assert d % 2 == 1 and primitive(h) == h
    verify_sieve(h, c['irreducibility'])
    a, b = [fmpq(*v) for v in c['isolating_interval']]; assert 0 < a < b
    simple_root_taylor(h, a, b)
    A, B = fmpz_poly(h.coeffs()[::2]), fmpz_poly(h.coeffs()[1::2])
    g = fmpz_poly([0,1])*B*B-A*A
    assert g.degree() == d and primitive(g) == g and B.degree() < d
    reflected = fmpz_poly([int(v)*(-1)**j for j, v in enumerate(h.coeffs())])
    substituted = [0]*(2*d+1)
    for j, v in enumerate(g.coeffs()): substituted[2*j] = int(v)
    assert h*reflected == -fmpz_poly(substituted)
    uniqueness = simple_root_taylor(g, a*a, b*b)
    result = {'format': 'Sn_odd_degree_square_corollary_v1', 'n': c['n'], 'degree': d,
              'coefficient_order': 'constant_upwards', 'polynomial_coefficients': list(map(int, g.coeffs())),
              'isolating_interval': [str(a*a), str(b*b)], 'root_uniqueness': uniqueness,
              'same_number_field_as_unsquared_value': True,
              'rational_recovery': 'alpha=-A(alpha^2)/B(alpha^2)',
              'A_coefficients': list(map(int, A.coeffs())), 'B_coefficients': list(map(int, B.coeffs())),
              'source_file': source.name, 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'proof': '[Q(alpha):Q(alpha^2)] is 1 or 2 and divides the odd degree d; hence it is 1. The displayed degree-d annihilator is therefore minimal.'}
    if output: Path(output).write_text(json.dumps(result, separators=(',', ':')))
    print('ODD_DEGREE_SQUARE_MINPOLY_VERIFIED', d, flush=True)
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('source', type=Path, nargs='?', default=Path('S7_minimal_polynomial_and_interval.json'))
    ap.add_argument('--output', type=Path, default=Path('S7_square_minimal_polynomial.json')); args = ap.parse_args()
    produce(args.source, args.output)
