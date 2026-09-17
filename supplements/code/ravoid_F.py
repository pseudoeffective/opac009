import sys
import numpy as np
from scipy.optimize import linprog
from schurcone import *
from hybrid import in_R
A = canon(eval(sys.argv[1])); rho = tuple(eval(sys.argv[2]))
lam = phi(A); N = sum(lam)
parts = list(partitions(N)); idx = {p: i for i, p in enumerate(parts)}; d = len(parts)
gens = [B for B in nested_multisets(N) if rho not in B and dominates(phi(B), lam) and in_R(B, A)]
print(f"A={A} rho={rho} lam={lam}; {len(gens)} R-competitors not containing rho:")
for B in gens: print("   ", B, phi(B))
A_eq = np.array([vec(mul_schur2({nu: 1}, rho), idx) for nu in partitions(N - sum(rho))], dtype=float)
A_ub = -np.array([vec(schur_product(B), idx) for B in gens], dtype=float)
# minimise sum |F| : F = P - Q, P,Q >= 0
c = np.ones(2*d)
A_eq2 = np.hstack([A_eq, -A_eq]); A_ub2 = np.hstack([A_ub, -A_ub])
res = linprog(c=c, A_ub=A_ub2, b_ub=-np.ones(len(gens)), A_eq=A_eq2, b_eq=np.zeros(len(A_eq)), bounds=[(0,None)]*(2*d), method="highs")
print("status", res.status)
if res.status == 0:
    F = res.x[:d] - res.x[d:]
    print("F (min L1):")
    for i in range(d):
        if abs(F[i]) > 1e-7: print(f"   {F[i]:+.4f} s{parts[i]}")
    vals = [(B, float(np.dot(vec(schur_product(B), idx), F))) for B in gens]
    print("values on R-competitors:", [(B, round(v,3)) for B, v in vals])
