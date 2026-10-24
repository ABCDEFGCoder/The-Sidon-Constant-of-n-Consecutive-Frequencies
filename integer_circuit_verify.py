"""Denominator-cleared integer RUR identities with homogeneous circuits.

Each intermediate represents D^degree times its polynomial evaluated at
N_i/D. Only integer polynomial arithmetic is used. Exact divisibility is
checked again by multiplication. No modular or numerical stopping rule is
accepted in place of this characteristic-zero identity.
"""
from collections import Counter
from functools import lru_cache
import time
import sympy as s
from flint import fmpz_poly, fmpz_mod_poly_ctx
from math import lcm
from exact_values import primitive


def integers(raw, names):
    scale = lcm(*(int(raw['numerators'][name][1]) for name in names if name != raw['names'][-1]))
    D = fmpz_poly(raw['common_denominator'])*scale
    N = {name: fmpz_poly(raw['numerators'][name][0])*(scale//int(raw['numerators'][name][1]))
         for name in names if name != raw['names'][-1]}
    if raw['names'][-1] in names: N[raw['names'][-1]] = fmpz_poly([0, 1])*D
    return D, N


def numerator(poly, variables, N, D):
    poly = s.Poly(poly, *variables).clear_denoms()[1].as_expr()
    replacements, outputs = s.cse([poly], symbols=s.numbered_symbols('_c'), order='none')
    values = {v: (N[str(v)], 1) for v in variables if str(v) in N}
    definitions = dict(replacements); uses = Counter()
    for expression in [f for _, f in replacements]+outputs:
        uses.update(expression.free_symbols & definitions.keys())

    @lru_cache(maxsize=6)
    def denominator_power(k): return D**k

    @lru_cache(maxsize=None)
    def degree(f):
        if f.is_Integer: return 0
        if f.is_Symbol:
            if f in values: return values[f][1]
            return degree(definitions[f])
        if f.is_Add: return max(map(degree, f.args), default=0)
        if f.is_Mul: return sum(map(degree, f.args))
        if f.is_Pow and f.exp.is_Integer and f.exp >= 0: return degree(f.base)*int(f.exp)
        raise ValueError('The verification circuit must be an integer polynomial')

    def evaluate(f):
        if f.is_Integer: return fmpz_poly([int(f)]), 0
        if f.is_Symbol:
            if f not in values:
                expression = definitions[f]; values[f] = evaluate(expression)
                release(expression)
            return values[f]
        if f.is_Add:
            d = degree(f); result = fmpz_poly([])
            for term in f.args:
                value, k = evaluate(term)
                if d != k: value = value*denominator_power(d-k)
                result = result+value
            return result, d
        if f.is_Mul:
            result = fmpz_poly([1]); d = 0
            for factor in sorted(f.args, key=degree):
                value, k = evaluate(factor); result = result*value; d += k
            return result, d
        if f.is_Pow and f.exp.is_Integer and f.exp >= 0:
            value, k = evaluate(f.base); return value**int(f.exp), k*int(f.exp)
        raise ValueError('Unsupported verification circuit')

    def release(expression):
        for symbol in expression.free_symbols & definitions.keys():
            uses[symbol] -= 1
            if uses[symbol] == 0: values.pop(symbol, None)

    result, d = evaluate(outputs[0]); release(outputs[0])
    values.clear(); denominator_power.cache_clear(); degree.cache_clear()
    return result, d, len(replacements)


def quotient_proposal(num, divisor):
    """Propose an integer quotient modulo one large prime power.

    Only the subsequent exact multiplication is evidence of divisibility.
    This avoids a long sequence of large integer coefficient divisions.
    """
    if num.is_zero(): return fmpz_poly([])
    bits = max(abs(int(v)).bit_length() for v in num.coeffs())+num.degree()+64
    p = 2
    while int(divisor[divisor.degree()]) % p == 0: p = int(s.nextprime(p))
    if p == 2: modulus = 1 << bits
    else:
        exponent = (bits+p.bit_length()-2)//(p.bit_length()-1)
        modulus = p**exponent
    context = fmpz_mod_poly_ctx(modulus)
    modular = context(num.coeffs())//context(divisor.coeffs())
    half = modulus//2
    return fmpz_poly([int(c)-modulus if int(c)>half else int(c) for c in modular.coeffs()])


def verify(raw, variables, polynomials):
    assert raw['characteristic'] == 0
    h = primitive(fmpz_poly(raw['parameter_polynomial']))
    assert h.gcd(fmpz_poly(raw['common_denominator'])).degree() == 0
    start = time.monotonic(); results = []
    for index, f in enumerate(polynomials):
        used = [v for v in variables if v in f.free_symbols]
        D, N = integers(raw, list(map(str, used)))
        num, degree, gates = numerator(f, used, N, D)
        print('INTEGER_CIRCUIT_NUMERATOR_CONSTRUCTED', index, 'degree', num.degree(), 'seconds', time.monotonic()-start, flush=True)
        if num.degree() > 2000:
            try: quotient = quotient_proposal(num, h)
            except (ValueError, ZeroDivisionError): quotient = num//h
            if quotient*h != num: quotient = num//h
        else: quotient = num//h
        assert quotient*h == num, ('Exact integer circuit RUR identity failed', index)
        row = {'index': index, 'homogeneous_degree': degree, 'numerator_degree': num.degree(),
               'CSE_assignments': gates, 'seconds_from_start': time.monotonic()-start}
        results.append(row); print('INTEGER_CIRCUIT_RUR_IDENTITY', row, flush=True)
    return results
