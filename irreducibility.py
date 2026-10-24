"""Short exact irreducibility proofs from several modular degree patterns.

A rational factor of a primitive integer polynomial keeps its degree at
every prime that preserves the leading coefficient. Its degree must be a
subset sum of the irreducible modular factor degrees. If the intersection
is {0,D}, no proper rational factor exists. Modular factors are checked by
the Rabin criterion, without trusting a factorization pass flag.
"""
from fractions import Fraction
from math import isqrt
from flint import fmpz_poly,nmod_poly
if not __debug__:raise SystemExit('Irreducibility verification requires assertions; do not use python -O.')

def prime(n):
 return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def prime_divisors(n):
 out=[];d=2
 while d*d<=n:
  if n%d==0:
   out.append(d)
   while n%d==0:n//=d
  d+=1
 if n>1:out.append(n)
 return out

def irreducible_mod(f,p):
 d=f.degree();assert d>=1 and int(f[d])==1
 x=nmod_poly([0,1],p)%f;r=x;checks={d//q for q in prime_divisors(d)}
 for k in range(1,d+1):
  r=r.pow_mod(p,f)
  if k in checks:assert (r-x).gcd(f).degree()==0
 assert r==x
 return True

def degree_sieve(polynomial,primes=(1000033,1000037,1000039)):
 h=fmpz_poly(polynomial);rows=[]
 for p in primes:
  f=nmod_poly(h.coeffs(),p);unit,factors=f.factor()
  rows.append({'prime':p,'unit':int(unit),'factors':[[list(map(int,g.coeffs())),int(e)] for g,e in factors]})
 c={'format':'modular_degree_intersection_v1','degree':h.degree(),'rows':rows}
 verify_sieve(h,c)
 return c

def verify_sieve(polynomial,certificate):
 h=fmpz_poly(polynomial);D=h.degree();c=certificate
 assert c['format']=='modular_degree_intersection_v1' and c['degree']==D
 possible=set(range(D+1));patterns=[]
 for row in c['rows']:
  p=row['prime'];assert type(p) is int and prime(p)
  f=nmod_poly(h.coeffs(),p);assert f.degree()==D and f.gcd(f.derivative()).degree()==0
  product=nmod_poly([row['unit']],p);subsets={0};degrees=[]
  for coefficients,multiplicity in row['factors']:
   assert multiplicity==1;g=nmod_poly(coefficients,p);irreducible_mod(g,p);product*=g
   d=g.degree();degrees.append(d);subsets|={v+d for v in list(subsets)}
  assert product==f
  possible&=subsets;patterns.append(degrees)
 assert possible=={0,D},('Degree sieve did not prove irreducibility',sorted(possible))
 return patterns

def nonsquare_norm(polynomial):
 h=fmpz_poly(polynomial);v=Fraction((-1)**h.degree()*int(h[0]),int(h[h.degree()]))
 square=v>=0 and isqrt(v.numerator)**2==v.numerator and isqrt(v.denominator)**2==v.denominator
 assert not square,'A square norm does not decide whether the element is a square'
 return [v.numerator,v.denominator]

def square_substitution_proof(polynomial):
 return {'argument':'irreducible_parameter_and_nonsquare_field_norm',
  'parameter_sieve':degree_sieve(polynomial),'norm':nonsquare_norm(polynomial)}

def verify_square_substitution(polynomial,certificate):
 assert certificate['argument']=='irreducible_parameter_and_nonsquare_field_norm'
 patterns=verify_sieve(polynomial,certificate['parameter_sieve'])
 assert nonsquare_norm(polynomial)==certificate['norm']
 # If alpha were a square in Q(alpha), its field norm would be a rational
 # square. Thus [Q(sqrt(alpha)):Q(alpha)]=2 and H(X^2) is irreducible.
 return patterns

if __name__=='__main__':
 import json,time
 from pathlib import Path
 path=Path(__file__).with_name('target6_rur_exact.json');h=json.loads(path.read_text())['parameter_polynomial'];start=time.monotonic();c=square_substitution_proof(h)
 print('IRREDUCIBILITY_VERIFIED',2*(len(h)-1),'patterns',verify_square_substitution(h,c),'seconds',time.monotonic()-start)
 path.with_name('S6_minimal_polynomial_irreducibility.json').write_text(json.dumps(c,separators=(',',':')))
