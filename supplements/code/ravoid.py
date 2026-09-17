import sys, time
import numpy as np
from scipy.optimize import linprog
from schurcone import *
from hybrid import in_R

def R_avoidance_fast(A, rho, verbose=False):
    A = canon(A); lam = phi(A); N = sum(lam)
    parts = list(partitions(N)); idx = {p: i for i, p in enumerate(parts)}; d = len(parts)
    gens = [B for B in nested_multisets(N) if rho not in B and dominates(phi(B), lam) and in_R(B, A)]
    A_eq = np.array([vec(mul_schur2({nu: 1}, rho), idx) for nu in partitions(N - sum(rho))], dtype=float)
    b_eq = np.zeros(len(A_eq))
    if not gens: return True, 0
    A_ub = -np.array([vec(schur_product(B), idx) for B in gens], dtype=float); b_ub = -np.ones(len(gens))
    res = linprog(c=np.zeros(d), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)] * d, method="highs")
    return res.status == 0, len(gens)

if __name__ == "__main__":
    A = canon(eval(sys.argv[1]))
    for rho in sorted(set(A), reverse=True):
        t = time.time(); ok, ng = R_avoidance_fast(A, rho)
        print(A, "rho =", rho, "R-avoidance:", ok, f"({ng} R-competitors, {time.time()-t:.1f}s)", flush=True)
