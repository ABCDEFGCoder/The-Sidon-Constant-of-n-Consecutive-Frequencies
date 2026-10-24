"""Full symmetric Galois group from a verified large prime Frobenius cycle.

The external group-theoretic input is Jordan's classical prime-cycle
theorem: a primitive group containing a prime cycle fixing at least three
points contains A_n. See Jones, arXiv:1209.5169, Theorem 1.1.
"""
from pathlib import Path
from math import gcd, lcm
import argparse, hashlib, json
from flint import fmpz_poly
from irreducibility import prime, verify_sieve


def produce(source, output=None):
    source = Path(source); c = json.loads(source.read_text()); h = fmpz_poly(c['polynomial_coefficients']); n = h.degree()
    verify_sieve(h, c['irreducibility'])
    chosen = None
    for row in c['irreducibility']['rows']:
        lengths = [len(coefficients)-1 for coefficients, exponent in row['factors']]
        assert all(exponent == 1 for _, exponent in row['factors']) and sum(lengths) == n
        for p in lengths:
            if prime(p) and 2*p > n and n-p >= 3 and all(gcd(p, q) == 1 for q in lengths if q != p):
                if lengths.count(p) == 1 and (n-len(lengths)) % 2 == 1:
                    chosen = row['prime'], lengths, p; break
        if chosen: break
    assert chosen is not None, 'No sufficient symmetric-group witness was found'
    modulus_prime, lengths, p = chosen
    power = lcm(*(q for q in lengths if q != p)); assert gcd(power, p) == 1
    powered = []
    for q in lengths: powered.extend([q//gcd(q, power)]*gcd(q, power))
    assert [q for q in powered if q > 1] == [p] and powered.count(1) == n-p
    blocks = [b for b in range(2, n) if n % b == 0]
    assert all(p > b and p > n//b for b in blocks)
    # Irreducibility gives transitivity. A block system is impossible:
    # an element of prime order p acts trivially on fewer than p blocks,
    # and then trivially inside every block of size smaller than p.
    # Thus G is primitive, Jordan gives A_n <= G, and the odd original
    # Frobenius permutation gives G=S_n.
    result = {'format': 'Sn_full_symmetric_galois_certificate_v1', 'n': c['n'],
              'polynomial_degree': n, 'galois_group': f'SymmetricGroup({n})',
              'good_prime': modulus_prime, 'frobenius_cycle_type': lengths,
              'isolated_prime_cycle': p, 'power_isolating_cycle': power, 'fixed_points': n-p,
              'odd_frobenius_permutation': True, 'possible_nontrivial_block_sizes': blocks,
              'jordan_theorem': 'Primitive group + prime cycle fixing >=3 points implies A_n is contained',
              'theorem_reference': 'Gareth A. Jones, arXiv:1209.5169, Theorem 1.1',
              'theorem_url': 'https://arxiv.org/abs/1209.5169',
              'source_file': source.name, 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'corollaries': {'not_expressible_by_radicals_over_Q': True,
                              'root_field_has_no_proper_intermediate_fields': True,
                              'normal_closure_degree': f'{n}!'}}
    if output: Path(output).write_text(json.dumps(result, separators=(',', ':')))
    print('FULL_SYMMETRIC_GALOIS_GROUP_VERIFIED', n, 'prime cycle', p, 'Frobenius', modulus_prime, flush=True)
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('source', type=Path, nargs='?', default=Path('S7_minimal_polynomial_and_interval.json'))
    ap.add_argument('--output', type=Path, default=Path('S7_galois_group_certificate.json')); args = ap.parse_args()
    produce(args.source, args.output)
