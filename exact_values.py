"""Exact primitive polynomials and local univariate root tests."""
from fractions import Fraction
from math import gcd,lcm
from functools import reduce
from flint import fmpq,fmpq_poly,fmpz_poly,arb,arb_poly,ctx
import json,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
if not __debug__:raise SystemExit('Proof verification requires assertions; do not use -O.')

def primitive(poly,positive_leading=True):
 cs=[Fraction(str(c)) for c in poly.coeffs()];den=lcm(*(c.denominator for c in cs)) if cs else 1
 a=[c.numerator*(den//c.denominator) for c in cs];g=reduce(gcd,a,0)
 if not g:return fmpz_poly([])
 if positive_leading and a[-1]<0:g=-g
 return fmpz_poly([v//g for v in a])

def simple_root_taylor(poly,a,b,bits=None):
 """Strict endpoint signs and a nonzero derivative on a rational interval.

 A local Taylor shift avoids the coefficient growth of a full high-degree
 Sturm chain. All interval arithmetic is outward rounded.
 """
 a=fmpq(str(a));b=fmpq(str(b));assert a<b;c=(a+b)/2;r=(b-a)/2
 cs=list(map(int,poly.coeffs())) if isinstance(poly,fmpz_poly) else poly.coeffs()
 height=max((abs(int(v)).bit_length() for v in cs if v),default=1)
 suggested=height+int(r.denominator).bit_length()+4*len(cs)+512
 previous=ctx.prec;ctx.prec=max(previous,suggested if bits is None else bits)
 try:
  shifted=arb_poly([arb(v) for v in cs])(arb_poly([arb(c),arb(1)]));left=shifted(arb(-r));right=shifted(arb(r));derivative=shifted.derivative()(arb(0,r))
  assert (left.upper()<0 and right.lower()>0) or (left.lower()>0 and right.upper()<0),'Endpoint signs not certified'
  assert not derivative.contains(0),'Derivative is not separated from zero'
  return {'method':'strict_Taylor_monotonicity','precision_bits':ctx.prec,'derivative_sign':1 if derivative.lower()>0 else -1}
 finally:ctx.prec=previous

def simple_root_rational(poly,a,b):
 """The same local criterion entirely in rational polynomial arithmetic."""
 a=fmpq(str(a));b=fmpq(str(b));assert a<b;c=(a+b)/2;r=(b-a)/2;p=fmpq_poly(poly);assert p(a)*p(b)<0
 shifted=p(fmpq_poly([c,1]));tail=fmpq(0)
 for j in range(shifted.degree(),1,-1):tail=(tail+j*abs(shifted[j]))*r
 assert abs(shifted[1])>tail,'Exact Taylor derivative test failed'
 return {'method':'rational_Taylor_monotonicity','derivative_sign':1 if shifted[1]>0 else -1}

def load_rur(path):
 from pathlib import Path
 return json.loads(Path(path).read_text())
