"""Exact degree-six chart identities and rational inverse."""
from invariant_critical6 import model

def identities():
 import sympy as s
 b,c,t=s.symbols('b c t',nonzero=True);B=b*b
 P=[(s.Integer(1),s.Integer(0)),(t*(2*c-b),b),((b*c-B/4)*(1-t*t)-t*t,t*(2+2*b*c-B/2)),(-b*t*(b*c+2),c*(B-4)*(1-t*t)/2-2*b*t*t)]
 P+=P[-2::-1];q1=P[1][0];q2=P[2][0]+B/2;Q=[1,q1,q2,0,-q2,-q1,-1]
 for lag in range(1,7):assert s.expand(sum(P[k][0]*P[k+lag][0]+P[k][1]*P[k+lag][1]+Q[k]*Q[k+lag] for k in range(7-lag)))==0
 E=s.factor(sum(u*u+v*v+q*q for (u,v),q in zip(P,Q)));assert s.factor(E-(B+4)**2*(1+t*t)**2*(1+c*c)/4)==0
 vs,fs,Ev,R,sat=model()
 for i,r in enumerate(R):assert s.factor(r.subs({vs[0]:B,vs[1]:b*c,vs[2]:t*t})-(P[i+1][0]**2+P[i+1][1]**2))==0
 a,e=s.symbols('a e');tt=2*(e-a*b)/(B+4);cc=(a+b*tt)/(2*tt)
 assert s.factor(P[1][0].subs({c:cc,t:tt})-a)==0
 assert s.factor(P[2][1].subs({c:cc,t:tt})-e)==0
 f=a*B/2-b*e;det=(4-B)*e-2*b*(f+2*a)
 assert s.expand(det-(B+4)*(e-a*b))==0
 # In an arbitrary normalized palindromic point, the three highest
 # nontrivial lags force q1=a, q2=d+B/2, and f=a*B/2-b*e.
 # The two remaining lags are linear in d,g, with this nonzero
 # determinant. Thus the displayed family covers a neighborhood.
 d,g=s.symbols('d g');PP=[(1,0),(a,b),(d,e),(f,g),(d,e),(a,b),(1,0)]
 QQ=[1,a,d+B/2,0,-d-B/2,-a,-1]
 lag=lambda h:s.expand(sum(PP[k][0]*PP[k+h][0]+PP[k][1]*PP[k+h][1]+QQ[k]*QQ[k+h] for k in range(7-h)))
 assert all(lag(h)==0 for h in [3,4,5,6])
 assert s.expand(lag(2)-((4-B)*d+2*b*g+2*a*f+e*e+B-B*B/4))==0
 assert s.expand(lag(1)/2-((2*a+f)*d+e*g+2*a+b*e+a*B/2))==0
 return True
