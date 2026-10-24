"""Independent verifier of the three exact closures in the accompanying paper.

The default runs every full global cover and every integer RUR identity.
No cached success flag can skip a proof obligation. Selective modes say
explicitly which layers they have rechecked.
"""
from pathlib import Path
from math import gcd
from functools import reduce
import argparse,contextlib,hashlib,json,os,subprocess,sys,time
from flint import fmpz_poly,nmod_poly,fmpq

ROOT=Path(__file__).resolve().parent
if not __debug__:raise SystemExit('Do not disable assertions with Python -O.')
sys.set_int_max_str_digits(0)

def check_polynomial_text():
    from exact_values import primitive,simple_root_rational,simple_root_taylor
    from irreducibility import irreducible_mod
    raw=(ROOT/'S6_S7_Minimal_Polynomials.txt').read_text()
    proofs={}
    for n in (6,7):
        section=raw.split(f'BEGIN S{n}\n')[1].split(f'END S{n}')[0]
        head,coefficients=section.split('COEFFICIENTS\n')
        data={line.split(' ',1)[0]:line.split(' ',1)[1] for line in head.strip().splitlines()}
        coefficients=json.loads(coefficients)
        assert all(type(v) is int for v in coefficients)
        p=fmpz_poly(coefficients)
        assert p.degree()==int(data['DEGREE'])==(372 if n==6 else 855)
        assert primitive(p)==p
        interval=json.loads(data['INTERVAL'])
        a,b=[fmpq(*v) for v in interval]
        assert (a,b)==((fmpq(2362563342647160829,10**18),fmpq(2362563342647160830,10**18))
                      if n==6 else (fmpq(2551658061766431151,10**18),fmpq(2551658061766431152,10**18)))
        if n==6:
            rur=json.loads((ROOT/'target6_rur_exact.json').read_text())
            source=[int(v) for v in rur['parameter_polynomial']]
            assert coefficients[::2]==source and all(v==0 for v in coefficients[1::2])
            assert max(abs(v).bit_length() for v in coefficients)>0
            uniqueness=simple_root_rational(p,a,b)
        else:
            rur=json.loads((ROOT/'S7_value_RUR_EXACT.json').read_text())
            assert coefficients==rur['parameter_polynomial']
            declared=json.loads((ROOT/'S7_minimal_polynomial_and_interval.json').read_text())
            assert coefficients==declared['polynomial_coefficients'] and interval==declared['isolating_interval']
            assert max(abs(v).bit_length() for v in coefficients)==2641
            uniqueness=simple_root_taylor(p,a,b)
        H=fmpz_poly(rur['parameter_polynomial'])
        assert fmpz_poly(rur['common_denominator'])==H.derivative()
        proofs[str(n)]={'degree':p.degree(),'interval':interval,'unique_root':uniqueness}
    h5=fmpz_poly([677,-414,-936,-1026,675])
    assert primitive(h5)==h5
    f=nmod_poly(h5.coeffs(),7)*5
    assert list(map(int,f.coeffs()))==[4,2,3,1,1]
    irreducible_mod(f,7)
    proofs['5']={'degree':4,'interval':[[218,100],[219,100]],
                 'unique_root':simple_root_rational(h5,fmpq(218,100),fmpq(219,100)),
                 'irreducible_modulus':7}
    print('POLYNOMIAL_TEXT_AND_ALL_PRINTED_INTERVALS_VERIFIED',flush=True)
    return proofs

def check_updated_interfaces():
    """A fast, explicitly incomplete regression mode; never a full proof run."""
    from chart6_geometry import identities
    from certify_s7_value_rur import real_branch
    from critical_bridge import check_printed_box
    from irreducibility import verify_square_substitution,verify_sieve
    assert identities()
    check_printed_box(6)
    raw6=json.loads((ROOT/'target6_rur_exact.json').read_text())
    verify_square_substitution(raw6['parameter_polynomial'],json.loads((ROOT/'S6_minimal_polynomial_irreducibility.json').read_text()))
    raw7=json.loads((ROOT/'S7_value_RUR_EXACT.json').read_text())
    p7=json.loads((ROOT/'S7_minimal_polynomial_and_interval.json').read_text())
    verify_sieve(raw7['parameter_polynomial'],p7['irreducibility'])
    assert real_branch(raw7)['accepted']
    print('UPDATED_INTERFACES_VERIFIED; full covers and integer identities were not rerun by this mode.',flush=True)

def run_program(name,args=()):
    print('RUN_FULL_COMPONENT',name,flush=True)
    subprocess.run([sys.executable,str(ROOT/name),*args],cwd=ROOT,check=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--part',choices=['all','global','algebraic','interfaces'],default='all')
    args=parser.parse_args(); os.chdir(ROOT); start=time.monotonic()
    if args.part in ('all','global'):
        run_program('Sidon5_Exact_Closure_Certificate_Revised.py')
        run_program('Sidon6_7_Exact_Closure_Certificate.py')
    if args.part in ('all','algebraic','interfaces'):
        check_polynomial_text()
        if args.part=='interfaces':
            check_updated_interfaces()
        else:
            run_program('verify_s6_algebraic.py')
            run_program('certify_s7_value_rur.py',('S7_value_RUR_EXACT.json','--verify-only'))
        from odd_degree_square import produce as square
        from galois_symmetric_certificate import produce as galois
        source=ROOT/'S7_minimal_polynomial_and_interval.json'
        square(source); galois(source)
    message=('ALL_THREE_EXACT_CLOSURES_VERIFIED' if args.part=='all'
             else 'SELECTED_LAYER_VERIFIED: '+args.part)
    print(message,'seconds',round(time.monotonic()-start,3),flush=True)

if __name__=='__main__':main()
