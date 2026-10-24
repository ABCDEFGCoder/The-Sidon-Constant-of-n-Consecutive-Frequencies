"""Degree-372 minimal polynomial, full integer identities, global bridge."""
from pathlib import Path
import json,time
import sympy as s
from flint import fmpz_poly,fmpq,arb,arb_poly,ctx
from exact_values import primitive,simple_root_rational,simple_root_taylor
from irreducibility import verify_square_substitution
from invariant_critical6 import model
from integer_circuit_verify import verify as integer_identities
from critical_bridge import check_printed_box
from chart6_geometry import identities

ROOT=Path(__file__).resolve().parent

def verify():
    start=time.monotonic()
    raw=json.loads((ROOT/'target6_rur_exact.json').read_text())
    interval=[[2362563342647160829,10**18],[2362563342647160830,10**18]]
    irred=json.loads((ROOT/'S6_minimal_polynomial_irreducibility.json').read_text())
    h=fmpz_poly(raw['parameter_polynomial'])
    assert primitive(h)==h and h.degree()==186 and raw['characteristic']==0
    assert raw['names']==['B','C','U','r1','r2','r3','V']
    assert fmpz_poly(raw['common_denominator'])==h.derivative()
    hp=fmpz_poly([int(h[k//2]) if k%2==0 else 0 for k in range(373)])
    patterns=verify_square_substitution(h,irred)
    def valuation5(value):
        value=abs(int(value)); assert value
        exponent=0
        while value%5==0: value//=5; exponent+=1
        return exponent
    assert valuation5(h[0])-valuation5(h[h.degree()])==-75
    a,b=[fmpq(*v) for v in interval]
    root_proof=simple_root_rational(hp,a,b)
    v0,v1=a*a,b*b
    simple_root_taylor(h,v0,v1,bits=8192)
    variables,equations,E,R,sat=model(); B,C,U,r1,r2,r3=variables; V=s.Symbol('V')
    L=2+2*r1+2*r2+r3
    target=4*B*L*L-V*(B+4)**2*(1+U)**2*(B+C*C)
    rows=integer_identities(raw,variables+[V],equations+[target])
    assert identities()
    box,bounds=check_printed_box(6)
    ctx.prec=8192; centre=(v0+v1)/2; radius=(v1-v0)/2
    shifted=arb_poly([arb(v) for v in h.coeffs()])(arb_poly([arb(centre),arb(1)]))
    derivative=shifted.derivative(); parameter_box=arb(0,radius)
    assert not derivative(parameter_box).contains(0)
    for _ in range(24):
        lower,upper=parameter_box.lower().fmpq(),parameter_box.upper().fmpq()
        midpoint=(lower+upper)/2
        newton=arb(midpoint)-shifted(arb(midpoint))/derivative(parameter_box)
        lower=max(lower,newton.lower().fmpq()); upper=min(upper,newton.upper().fmpq())
        assert lower<upper
        parameter_box=arb((lower+upper)/2,(upper-lower)/2)
        if upper-lower<fmpq(1,10**90):break
    else:raise AssertionError('Interval Newton did not isolate the squared-value branch')
    def ev(coefficients):
        return arb_poly([arb(v) for v in coefficients])(arb_poly([arb(centre),arb(1)]))(parameter_box)
    den=ev(raw['common_denominator']); assert not den.contains(0)
    point=[]
    for name,X in zip(raw['names'][:-1],box):
        coefficients,scale=raw['numerators'][name]
        assert type(scale) is int and scale>0
        z=ev(coefficients)/(arb(scale)*den)
        assert X.lower()<z.lower() and z.upper()<X.upper(),('RUR branch outside printed box',name)
        point.append(z)
    assert all(point[k].lower()>0 for k in (0,2,3,4,5))
    out={'accepted':True,'n':6,'degree':372,'squared_value_degree':186,
         'interval':interval,'root_proof':root_proof,
         'factor_patterns':patterns,'critical_box':bounds,'integer_identities':rows,
         'seconds':time.monotonic()-start}
    print('S6_ALGEBRAIC_CLOSURE_VERIFIED',out,flush=True)
    return out

if __name__=='__main__':
    verify()
