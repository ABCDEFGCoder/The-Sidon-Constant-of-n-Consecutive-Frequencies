"""Rational symmetric complementary charts, with no trigonometric coefficients.

The structural identity is proved by the two-by-two matrix product. The
chart is for discovering and certifying candidates; it is never an assumed
symmetry of all Sidon extremizers.
"""
from pathlib import Path
from functools import lru_cache
import sympy as s
import json,argparse,time

if not __debug__:raise SystemExit('Proof checks require assertions; do not use Python -O')

@lru_cache(None)
def chart(n):
    if n < 2: raise ValueError('n >= 2 required')
    m=n//2;t=s.symbols('t1:'+str(m+1),real=True);z=s.Symbol('z')
    W=s.Matrix([[1+z,1-z],[1-z,1+z]])/2
    U=W
    for j in range(1,n):
        x=t[min(j,n-j)-1]
        U=(U*s.diag(1+s.I*x,1-s.I*x)*W).applyfunc(s.expand)
    U=2*U
    P=[s.expand(U[0,0]).coeff(z,k) for k in range(n+1)]
    Q=[s.expand(U[0,1]).coeff(z,k) for k in range(n+1)]
    E=4*s.prod((1+t[j-1]**2)**(1 if n%2==0 and j==m else 2) for j in range(1,m+1))
    assert P[0]==P[-1]==Q[0]==1 and Q[-1]==-1
    assert P==P[::-1] and Q==[-v for v in Q[::-1]]
    assert all(s.expand(s.im(v))==0 for v in Q)
    return t,P,Q,E

def audit(n):
    t,P,Q,E=chart(n)
    for lag in range(n+1):
        corr=s.expand(sum(P[k]*s.conjugate(P[k+lag])+Q[k]*Q[k+lag] for k in range(n+1-lag)))
        assert s.expand(corr-(E if lag==0 else 0))==0,lag
    return {'n':n,'parameters':len(t),'coefficient_terms':sum(len(s.Add.make_args(v)) for v in P),'energy':str(E)}

def inverse(P,Q):
    """Peel two degrees at a time, entirely by rational operations.

    Inputs have P=P^rev, Q=-Q^rev real, P0=Q0=1 and satisfy the
    complementary identity. The caller is responsible for those premises.
    For exact inputs all dropped coefficients and final base cases are
    checked. The sole division at a step is by (1+t*t)^2, positive over R.
    """
    n=len(P)-1
    if n<1 or len(Q)!=len(P):raise ValueError('Incompatible degrees')
    A=[s.expand(s.re(v)) for v in P];B=[s.expand(s.im(v)) for v in P]
    assert A==A[::-1] and B==B[::-1] and Q==[-v for v in Q[::-1]]
    assert A[0]==Q[0]==1 and B[0]==0
    if n==1:
        assert P==[1,1] and Q==[1,-1]
        return []
    t=s.cancel(B[1]/2)
    if n==2:
        assert s.expand(A[1])==0 and Q==[1,0,-1]
        return [t]
    assert s.expand(A[1]-Q[1])==0
    a=[s.expand(A[k]+Q[k]) for k in range(n-1)]
    b=B[1:n];c=[s.expand(A[k+2]-Q[k+2]) for k in range(n-1)]
    den=(1+t*t)**2
    aa=[s.cancel((u+2*t*v-t*t*w)/den) for u,v,w in zip(a,b,c)]
    cc=[s.cancel((w+2*t*v-t*t*u)/den) for u,v,w in zip(a,b,c)]
    bb=[s.cancel(((1-t*t)*v-t*(u+w))/den) for u,v,w in zip(a,b,c)]
    PP=[s.cancel((u+w)/2+s.I*v) for u,w,v in zip(aa,cc,bb)]
    QQ=[s.cancel((u-w)/2) for u,w in zip(aa,cc)]
    return [t]+inverse(PP,QQ)

def critical(n):
    t,P,Q,E=chart(n);m=len(t);r=s.symbols('r1:'+str(m+1));weights=[2]*m
    if n%2==0:weights[-1]=1
    R=[s.expand(P[j]*s.conjugate(P[j])) for j in range(1,m+1)]
    fs=[s.expand(rj*rj-v) for rj,v in zip(r,R)]
    for j,x in enumerate(t):
        a=2 if n%2==0 and j==m-1 else 4
        f=sum(w*((1+x*x)*s.diff(v,x)-a*x*v)*s.prod(r[k] for k in range(m) if k!=i) for i,(w,v) in enumerate(zip(weights,R)))-2*a*x*s.prod(r)
        fs.append(s.expand(f))
    sat=s.prod(1+x*x for x in t)*s.prod(r)
    L=2+sum(w*v for w,v in zip(weights,r))
    return list(t)+list(r),fs,E,R,sat,L

