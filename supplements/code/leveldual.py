"""Find a Farkas functional y on monomials with y.P_A >= 1, y.P_M <= 0 (level competitors), minimising L1 norm."""
import sys
import numpy as np
from scipy.optimize import linprog
from schurcone import *
from epsjet import inside_pairs
from levelid import poly_of

def dual(A, b0, a0):
    A = canon(A); lam = phi(A); N = sum(lam)
    comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
    T = list(range(b0, a0 + 2))
    AI = tuple(sorted(inside_pairs(A, b0, a0), reverse=True)); n = len(AI)
    Ms = sorted(set(tuple(sorted(inside_pairs(B, b0, a0), reverse=True)) for B in comps if len(inside_pairs(B, b0, a0)) == n))
    PA = poly_of(AI, T); PMs = [poly_of(M, T) for M in Ms]
    monos = sorted(set(PA) | set(m for p in PMs for m in p)); idx = {m: i for i, m in enumerate(monos)}; d = len(monos)
    def v(p):
        x = np.zeros(d)
        for m, c in p.items(): x[idx[m]] = c
        return x
    # variables y = yp - yn >= 0 ; constraints: -v(PA).y <= -1 ; v(PM).y <= 0
    rows = [-v(PA)] + [v(p) for p in PMs]; rhs = [-1.0] + [0.0]*len(PMs)
    Aub = np.array([np.concatenate([r, -r]) for r in rows]); c = np.ones(2*d)
    res = linprog(c=c, A_ub=Aub, b_ub=np.array(rhs), bounds=[(0,None)]*(2*d), method="highs")
    print(f"A={A} I=[{b0},{a0}] A_I={AI}; level M's:")
    for M in Ms: print("    ", M)
    if res.status != 0:
        print("  no functional (P_A in cone)"); return
    y = res.x[:d] - res.x[d:]
    print("  functional y (monomial -> weight):")
    for i in range(d):
        if abs(y[i]) > 1e-9: print(f"     {y[i]:+.3f} * delta{monos[i]}")
    print("  values: A ->", round(float(v(PA)@y),3), "; M's ->", [round(float(v(p)@y),3) for p in PMs])

if __name__ == "__main__":
    A = eval(sys.argv[1]); b0 = int(sys.argv[2]); a0 = int(sys.argv[3])
    dual(A, b0, a0)
