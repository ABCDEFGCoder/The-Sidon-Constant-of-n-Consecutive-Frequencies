"""Canonical integer-coefficient export; never emit coefficient divisions."""
from pathlib import Path
import sympy as s
import hashlib

def export(path,variables,polynomials,prime=0):
 canonical=[s.Poly(f,*variables).clear_denoms()[1].as_expr() for f in polynomials]
 assert all(all(c.q==1 for c in s.Poly(f,*variables).coeffs()) for f in canonical)
 rows=[str(s.expand(f)).replace('**','^').replace(' ','') for f in canonical]
 assert all('/' not in row and '.' not in row for row in rows)
 body=','.join(map(str,variables))+'\n'+str(prime)+'\n'+',\n'.join(rows);Path(path).write_text(body)
 return hashlib.sha256(body.encode()).hexdigest()
