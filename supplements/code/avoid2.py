"""Avoidance with F supported on partitions mu with mu_2 < p (automatically in ker s_(p,p)^perp)."""
import sys, time
from fractions import Fraction
import cdd.gmp as cdd
from schurcone import *

def avoidance_small_support(N, p, verbose=False):
    rho = (p,p)
    sup = [mu for mu in partitions(N) if len(mu) < 2 or mu[1] < p]
    idx = {m:i for i,m in enumerate(sup)}; d = len(sup)
    gens = [B for B in nested_multisets(N) if rho not in B]
    rows = []
    for B in gens:
        f = schur_product(B)
        rows.append([Fraction(-1)] + [Fraction(f.get(m,0)) for m in sup])
    # objective: minimise sum of |F| -- use F = F+ - F-; simpler: minimise sum F (with F >= -M?) skip; just feasibility
    mat = cdd.matrix_from_array(rows, rep_type=cdd.RepType.INEQUALITY,
                                obj_type=cdd.LPObjType.MIN, obj_func=[Fraction(0)]*(d+1))
    lp = cdd.linprog_from_matrix(mat); cdd.linprog_solve(lp)
    ok = lp.status == cdd.LPStatusType.OPTIMAL
    sol = None
    if ok:
        sol = {sup[i]: lp.primal_solution[i] for i in range(d) if lp.primal_solution[i] != 0}
    return ok, sol, gens

if __name__ == "__main__":
    for N in range(4, int(sys.argv[1])+1):
        for p in range(1, N//2+1):
            ok, sol, gens = avoidance_small_support(N, p)
            print(f"N={N} p={p}: {'OK' if ok else 'NO'}  {sol if ok and N<=8 else ''}", flush=True)
