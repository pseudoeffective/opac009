"""Test: does cone{s_B : B nested, rho not in B} avoid the subspace s_rho * Lambda_{N-|rho|}?
Farkas: exists F in ker(s_rho^perp) with F(s_B) >= 1 for all such B."""
import sys, time
from fractions import Fraction
import cdd.gmp as cdd
from schurcone import *

def avoidance(N, rho, restrict_above=None):
    parts = list(partitions(N)); idx = {p:i for i,p in enumerate(parts)}; d = len(parts)
    gens = [B for B in nested_multisets(N) if rho not in B]
    if restrict_above is not None:
        gens = [B for B in gens if dominates(phi(B), restrict_above)]
    rows = []
    # equalities: <F, s_rho s_nu> = 0
    for nu in partitions(N - sum(rho)):
        v = vec(mul_schur2({nu:1}, rho), idx)
        rows.append([Fraction(0)] + [Fraction(x) for x in v])
    neq = len(rows)
    for B in gens:
        v = vec(schur_product(B), idx)
        rows.append([Fraction(-1)] + [Fraction(x) for x in v])   # -1 + <F,s_B> >= 0
    mat = cdd.matrix_from_array(rows, lin_set=range(neq), rep_type=cdd.RepType.INEQUALITY,
                                obj_type=cdd.LPObjType.MIN, obj_func=[Fraction(0)]*(d+1))
    lp = cdd.linprog_from_matrix(mat); cdd.linprog_solve(lp)
    return lp.status == cdd.LPStatusType.OPTIMAL, len(gens)

for N in range(4, int(sys.argv[1])+1):
    for p in range(1, N//2 + 1):
        t = time.time()
        ok, ng = avoidance(N, (p,p))
        print(f"N={N} rho=({p},{p}): avoidance {'HOLDS' if ok else 'FAILS'} ({ng} gens, {time.time()-t:.1f}s)", flush=True)
