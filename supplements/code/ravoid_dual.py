import sys
import numpy as np
from scipy.optimize import linprog
from schurcone import *
from hybrid import in_R
A = canon(eval(sys.argv[1])); rho = tuple(eval(sys.argv[2]))
lam = phi(A); N = sum(lam)
parts = list(partitions(N)); idx = {p: i for i, p in enumerate(parts)}; d = len(parts)
gens = [B for B in nested_multisets(N) if rho not in B and dominates(phi(B), lam) and in_R(B, A)]
nus = list(partitions(N - sum(rho)))
cols = [vec(schur_product(B), idx) for B in gens] + [[-x for x in vec(mul_schur2({nu: 1}, rho), idx)] for nu in nus]
M = np.array(cols, dtype=float).T
m, k = len(gens), len(nus)
A_eq = np.vstack([M, np.array([[1.0]*m + [0.0]*k])]); b_eq = np.zeros(d+1); b_eq[-1] = 1
bounds = [(0, None)]*m + [(None, None)]*k
res = linprog(c=np.zeros(m+k), A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
print("status", res.status)
if res.status == 0:
    for j in range(m):
        if res.x[j] > 1e-9: print(f"  {res.x[j]:.4f} * {gens[j]}")
    print(" = s_rho * (")
    for j in range(k):
        if abs(res.x[m+j]) > 1e-9: print(f"  {res.x[m+j]:.4f} s{nus[j]}")
    print(" )")
