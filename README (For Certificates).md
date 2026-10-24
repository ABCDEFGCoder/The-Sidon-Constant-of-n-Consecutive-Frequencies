# Exact Sidon constants S5, S6 and S7

This release accompanies *Exact Sidon Constants in Degrees Five, Six and
Seven and Strict Monotonicity in Every Degree* (24 October 2026).
The coefficients in the extremal problem are unrestricted complex numbers.

| Constant | Minimal polynomial degree | Unique rational root interval |
|---|---:|---|
| S5 | 4 | (218/100, 219/100) |
| S6 | 372 | (2362563342647160829/10^18, 2362563342647160830/10^18) |
| S7 | 855 | (2551658061766431151/10^18, 2551658061766431152/10^18) |

The complete coefficient lists for S6 and S7 are in
`S6_S7_Minimal_Polynomials.txt`, ordered from constant to leading term.
For S5 the polynomial is `675*X**4-1026*X**3-936*X**2-414*X+677`.

## Run the complete verification

Use Python 3.10 or later without `-O` (assertions are proof predicates):

```sh
python3 -m pip install -r requirements.txt
python3 verify_all.py
```

The default rechecks all global phase covers, all exact integer identities,
the real branch bridges, all three unique root intervals, irreducibility,
and the degree-seven field and Galois corollaries. No saved success flag
can bypass a check. The verifier itself uses no network and no optimizer.
The longest part is the degree-seven integer RUR verification. Keep these
large identities sequential; a machine with at least 8 GiB RAM is advisable.
Execution time depends on the host (typically several to tens of minutes).

Optional layers are explicit:

```sh
python3 verify_all.py --part global
python3 verify_all.py --part algebraic
python3 verify_all.py --part interfaces
```

`interfaces` is an incomplete regression check. It rechecks polynomial text,
printed intervals, modular irreducibility, the updated chart boxes and real
bridges, and structural corollaries. It does **not** rerun phase covers or
the expensive integer RUR identities, and does not claim full closure by
itself. The manuscript gives the analytic proof of strict monotonicity;
that theorem needs no numerical certificate.

## Proof chain

1. The independent standard-library global scripts establish feasible
   candidates, locally sharp Toeplitz bounds, full phase covers, and reversal.
2. `critical_bridge.py` maps their roots into the six-variable boxes printed
   in the manuscript. Its preconditioners are exact rational Jacobian
   inverses at those printed rational centers.
3. `verify_s6_algebraic.py` checks the degree-186 squared-value RUR, seven
   integer polynomial identities, root isolation and its global bridge. A
   nonsquare norm doubles the value degree to 372.
4. `certify_s7_value_rur.py` checks seven characteristic-zero RUR identities,
   the unique real branch and degree-855 irreducibility.
5. `verify_all.py` compares every long-polynomial coefficient with its RUR
   source and checks the field/Galois corollaries.

The proof matrices A,B,L,W and both stationary boxes are printed in the
manuscript. The large global-cover payloads, six coordinate numerator arrays
for each RUR, and modular factor coefficients are machine-readable witnesses
whose entire meaning and acceptance rules are specified in the paper.
They are not unverified numerical output. Quotients proposed modulo a large
prime power are accepted only after exact integer multiplication.

## Reproducibility

`verification_logs/` records the full rechecks carried out for this revision.
The original global scripts and their complete data are retained unchanged;
the updated algebraic interfaces were also checked separately after curation.
`MANIFEST_SHA256.json` records file hashes for integrity. The manifest and
logs are not substitutes for running the verifier.

The trusted arithmetic is Python arbitrary-precision integers and rationals,
FLINT exact polynomials, and outward Arb enclosures through python-flint.
SymPy supplies exact symbolic expressions. This is a computer-assisted proof,
not a proof-assistant formalization.

## Body-only adversarial review

The manuscript was reviewed and tightened after the initial full verification.
The updated 38-page paper, `paper/Sidon_5_to_7_Body_Audit.txt` and
`paper/Body_Audit_Changes.diff` are synchronized in this archive.
No certificate was rerun for this body review. All proof programs, RUR and
other algebraic witnesses, long polynomial lists, and prior verification
logs remain unchanged; the earlier logs document the initial full rechecks.
