# OPAC-009 code (Python 3.11; pycddlib 3 for exact LPs, scipy for floating LPs)

Core: `schurcone.py` — partitions, dominance, Pieri/Jacobi–Trudi Schur products, White's interlacing/nested test,
enumeration of (nested) multisets, exact extremality LP (`is_extreme`), separating functionals.

## Scripts used for the main theorem (session 2)
* `verify_levelid.py N` — sanity check of the Level Identity Theorem: for every non-nested A with weight ≤ N,
  finds an exact decomposition and verifies (i) n_I(B) ≥ n_I(A) for all B in the support and (ii) the polynomial
  identity K_A P_A = Σ c_B K_B P_B, for every interval I (run: N ≤ 9, all OK).
* `oddmono.py N` — checks the combinatorial lemma of the proof: for every nested A (weight ≤ N) and every nested
  competitor B ≠ A with φ(B) ⊵ φ(A) satisfying all interval inequalities, the monomial of Step 3/4 has
  coefficient 0 (run: N ≤ 18, 0 failures).
* `levelid.py`, `lstar.py`, `epsjet.py` — polynomial cone test P_A ∉ cone{P_B} and the restricted competitor set
  L*(A); `poly_of` expands ∏(δ_b − δ_{a+1}).
* `oddcase.py`, `leveldual.py`, `bottommono*.py`, `monofind.py`, `monotest.py` — exploration of the odd case.

## Session-1 scripts
* `formulaF.py`, `formulaF2.py`, `dcount.py`, `dtest.py`, `onebox.py` — one-box LR formula (Theorem F) and d vs c.
* `vertices.py`, `witness.py`, `check_ext.py`, `sweep.py`, `table*.py` — certificate polytopes, exact checks.
* `avoid*.py`, `ravoid*.py`, `hybrid*.py`, `closure.py`, `coverage.py`, `verify_thm1.py` — earlier reduction ideas
  (restricted avoidance etc.), mostly negative results; `fastlp.py` floating-point LP helpers.
* `signcrit.py`, `signtest.py`, `signsearch.py`, `tropical*.py`, `coarse*.py`, `chaintest.py`, `levelsign.py`,
  `levelid2.py`, `lstar2-4.py` — sign/tropical criteria tried before the monomial argument (weaker).
