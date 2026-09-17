"""Restricted avoidance: for nested A and rho in A, does there exist F in ker(s_rho^perp) (F in Lambda_N)
with F(s_B) >= 1 for all nested B with rho not in B and phi(B) ⊵ phi(A)?"""
import sys, time
from fractions import Fraction
import cdd.gmp as cdd
from schurcone import *

def restricted_avoidance(A, rho):
    lam = phi(A); N = sum(lam)
    parts = list(partitions(N)); idx = {p:i for i,p in enumerate(parts)}; d = len(parts)
    gens = [B for B in nested_multisets(N) if rho not in B and dominates(phi(B), lam)]
    rows = []
    for nu in partitions(N - sum(rho)):
        v = vec(mul_schur2({nu:1}, rho), idx)
        rows.append([Fraction(0)] + [Fraction(x) for x in v])
    neq = len(rows)
    for B in gens:
        v = vec(schur_product(B), idx)
        rows.append([Fraction(-1)] + [Fraction(x) for x in v])
    mat = cdd.matrix_from_array(rows, lin_set=range(neq), rep_type=cdd.RepType.INEQUALITY,
                                obj_type=cdd.LPObjType.MIN, obj_func=[Fraction(0)]*(d+1))
    lp = cdd.linprog_from_matrix(mat); cdd.linprog_solve(lp)
    return lp.status == cdd.LPStatusType.OPTIMAL

if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(4, Nmax+1):
        gens = nested_multisets(N)
        t = time.time(); nA = 0; fails = []
        good_some = 0
        for A in gens:
            if len(A) < 2: continue
            nA += 1
            oks = {}
            for rho in set(A):
                oks[rho] = restricted_avoidance(A, rho)
            if any(oks.values()): good_some += 1
            for rho, ok in oks.items():
                if not ok: fails.append((A, rho))
        print(f"N={N}: {nA} A's with >=2 parts; {good_some} have some rho with restricted avoidance; "
              f"{len(fails)} (A,rho) failures  [{time.time()-t:.0f}s]", flush=True)
        for A, rho in fails[:12]:
            print("    FAIL", A, rho)
