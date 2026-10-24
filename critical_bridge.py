"""Exact bridges from the unrestricted global certificates to six variables.

The equalities used to set symmetric endpoint coefficients are consequences
of reversal uniqueness and full contact, proved in the accompanying paper.
Interval overlap by itself is never used as a proof of an equality.
"""
from pathlib import Path
import importlib.util
import json
import sympy as s
from flint import arb, arb_mat, fmpq, ctx

ROOT = Path(__file__).resolve().parent

def eval_poly(poly, variables, x):
    out = arb(0)
    for monomial, coefficient in s.Poly(poly, *variables).terms():
        value = arb(fmpq(int(s.numer(coefficient)), int(s.denom(coefficient))))
        for z, exponent in zip(x, monomial):
            if exponent:
                value *= z**exponent
        out += value
    return out

def norminf(A):
    return max(sum((A[i,j].abs_upper() for j in range(A.ncols())), arb(0)).upper().fmpq()
               for i in range(A.nrows()))

def old_module():
    spec = importlib.util.spec_from_file_location('sidon_global', ROOT/'Sidon6_7_Exact_Closure_Certificate.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def cm(a,b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]

def conj(a):
    return a[0], -a[1]

def magnitude(a):
    return (a[0]**2+a[1]**2).sqrt()

def reciprocal_source(n):
    assert n in (6,7)
    g = old_module(); d = g.ROOT_DATA[str(n)]
    root = g.krawczyk(n, d['point'], d['R'])
    g.verify_reversal_symmetry(n, root, d['point'], d['beta'])
    def cv(v):
        a,b = fmpq(v.lo,g.SCALE), fmpq(v.hi,g.SCALE)
        return arb((a+b)/2,(b-a)/2)
    p = [(cv(root[k]), arb(0) if k<2 else cv(root[n+1+k-2])) for k in range(n+1)]
    gamma = tuple(v/magnitude(p[-1]) for v in conj(p[-1]))
    eta = tuple(v/magnitude(p[-2]) for v in conj(cm(gamma,p[-2])))
    re = ((1+eta[0])/2).sqrt(); assert re.lower()>0
    zeta = (re,-eta[1]/(2*re)); zn = (arb(1),arb(0))
    for _ in range(n): zn=cm(zn,zeta)
    endpoint = cm(p[-1],zn)
    if n%2 and endpoint[0].upper()<0:
        zeta=(-zeta[0],-zeta[1])
    else:
        assert endpoint[0].lower()>0
    power=(arb(1),arb(0)); P=[]
    for a in p:
        z=cm(a,power); P.append((z[0]/p[0][0],z[1]/p[0][0])); power=cm(power,zeta)
    assert all((P[k][j]-P[n-k][j]).contains(0) for k in range(n+1) for j in (0,1))
    assert (P[-1][0]-1).contains(0) and P[-1][1].contains(0)
    return g,root,P,zeta,cv

def peel_interval(A,B,Q):
    n=len(A)-1
    if n==1:return []
    t=B[1]/2
    if n==2:return [t]
    a=[A[k]+Q[k] for k in range(n-1)]; b=B[1:n]
    c=[A[k+2]-Q[k+2] for k in range(n-1)]
    den=(1+t*t)**2; assert den.lower()>0
    aa=[(u+2*t*v-t*t*w)/den for u,v,w in zip(a,b,c)]
    cc=[(w+2*t*v-t*t*u)/den for u,v,w in zip(a,b,c)]
    bb=[((1-t*t)*v-t*(u+w))/den for u,v,w in zip(a,b,c)]
    AP=[(u+w)/2 for u,w in zip(aa,cc)]; QP=[(u-w)/2 for u,w in zip(aa,cc)]
    AP[0]=AP[-1]=arb(1); bb[0]=bb[-1]=arb(0); QP[0]=arb(1); QP[-1]=arb(-1)
    return [t]+peel_interval(AP,bb,QP)

def global_coordinates(n):
    g,root,P,zeta,cv=reciprocal_source(n)
    amplitudes=[magnitude(P[k]) for k in range(1,4)]
    if n==6:
        a,b=P[1]; e=P[2][1]
        t=2*(e-a*b)/(b*b+4)
        assert not b.contains(0) and not t.contains(0)
        return [b*b, b*b/2+a*b/(2*t),t*t]+amplitudes
    qq=[(arb(1),arb(0))]
    for theta in root[2*n:3*n]:
        sn,cs=g.itrign(theta); z=cm((cv(cs),cv(sn)),conj(zeta))
        out=[(arb(0),arb(0)) for _ in range(len(qq)+1)]
        for k,a in enumerate(qq):
            v=cm(a,z)
            out[k]=(out[k][0]-v[0],out[k][1]-v[1])
            out[k+1]=(out[k+1][0]+a[0],out[k+1][1]+a[1])
        qq=out
    q0=qq[0]; den=q0[0]**2+q0[1]**2
    QQ=[cm(a,conj(q0)) for a in qq]; QQ=[(a/den,b/den) for a,b in QQ]
    assert all(b.contains(0) for a,b in QQ)
    A=[a for a,b in P]; B=[b for a,b in P]; Q=[a for a,b in QQ]
    A[0]=A[-1]=arb(1); B[0]=B[-1]=arb(0); Q[0]=arb(1); Q[-1]=arb(-1)
    return peel_interval(A,B,Q)+amplitudes

def check_printed_box(n, global_bridge=True):
    """Preconditioner is the exact rational inverse at the printed centre."""
    if n==6:
        from invariant_critical6 import model
        variables,equations,*_=model()
    else:
        from qsp_chart import critical
        variables,equations,*_=critical(n)
    data=json.loads((ROOT/'critical_boxes.json').read_text())[str(n)]
    centre=[s.Rational(z,10**12) for z in data['centre_numerators']]
    radius=fmpq(1,10**8); ctx.prec=512
    J=s.Matrix(equations).jacobian(variables)
    J0=J.subs(dict(zip(variables,centre)))
    inverse=J0.inv()
    B=arb_mat(6,6,[arb(fmpq(int(z.p),int(z.q))) for z in inverse])
    xc=[arb(fmpq(int(z.p),int(z.q))) for z in centre]
    box=[arb(fmpq(int(z.p),int(z.q)),radius) for z in centre]
    residual=B*arb_mat(6,1,[eval_poly(f,variables,xc) for f in equations])
    A=arb_mat(6,6,[eval_poly(J[i,j],variables,box) for i in range(6) for j in range(6)])
    M=B*A
    for i in range(6):M[i,i]-=1
    eta,q=norminf(residual),norminf(M)
    assert q<1 and eta+q*radius<radius
    assert eta<fmpq(5,10**13)
    assert q<(fmpq(6,10**6) if n==6 else fmpq(15,10**6))
    if global_bridge:
        gx=global_coordinates(n)
        assert all(X.lower()<a.lower() and a.upper()<X.upper() for X,a in zip(box,gx))
        assert all(eval_poly(f,variables,gx).contains(0) for f in equations)
    print('PRINTED_CRITICAL_BOX_VERIFIED',n,'eta',float(eta),'q',float(q),flush=True)
    return box,{'eta':str(eta),'q':str(q),'radius':str(radius)}
