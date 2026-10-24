"""Six-variable stationary system after fourfold sign quotient.

Three positive amplitudes are retained. B=b^2, C=b*c, U=t^2.
This is a candidate chart, not a global symmetry assumption.
"""
import sympy as s
from pathlib import Path
import json

def model():
 B,C,U,r1,r2,r3=s.symbols('B C U r1 r2 r3');vs=[B,C,U,r1,r2,r3];rs=[r1,r2,r3];L=2+2*r1+2*r2+r3
 R=[B+U*(2*C-B)**2/B,((C-B/4)*(1-U)-U)**2+U*(2+2*C-B/2)**2,B*U*(C+2)**2+(C*(B-4)*(1-U)-4*B*U)**2/(4*B)]
 E=(B+4)**2*(1+U)**2*(B+C*C)/(4*B)
 fs=[s.fraction(s.cancel(r*r-v))[0] for r,v in zip(rs,R)]
 for x in [B,C,U]:
  dlog=s.cancel(s.diff(E,x)/E)
  # Use r_i^2=R_i to remove the L*rho products before expansion.
  f=sum(w*(s.diff(v,x)-v*dlog)*s.prod(rs[j] for j in range(3) if j!=i) for i,(w,v) in enumerate(zip([2,2,1],R)))-2*dlog*s.prod(rs)
  f=s.factor(s.fraction(s.cancel(f))[0]);fs.append(f)
 sat=B*(B+4)*(1+U)*(B+C*C)*r1*r2*r3
 fs=[s.Poly(f,*vs).clear_denoms()[1].as_expr() for f in fs]
 return vs,fs,E,R,sat

