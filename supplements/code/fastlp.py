"""Floating-point LP helpers (scipy HiGHS) for exploration; use exact cdd for final verification."""
import numpy as np
from scipy.optimize import linprog
from schurcone import *


def restricted_avoidance_fast(A, rho, N=None, gens=None):
    """Exists F on Lambda_N with <F, s_rho s_nu> = 0 for all nu, <F, s_B> >= 1 for nested B not containing rho,
    phi(B) ⊵ phi(A)?  Returns (feasible, F dict or None)."""
    A = canon(A); lam = phi(A)
    if N is None: N = sum(lam)
    parts = list(partitions(N)); idx = {p: i for i, p in enumerate(parts)}; d = len(parts)
    if gens is None:
        gens = [B for B in nested_multisets(N) if rho not in B and dominates(phi(B), lam)]
    A_eq = []
    for nu in partitions(N - sum(rho)):
        A_eq.append(vec(mul_schur2({nu: 1}, rho), idx))
    A_eq = np.array(A_eq, dtype=float); b_eq = np.zeros(len(A_eq))
    A_ub = -np.array([vec(schur_product(B), idx) for B in gens], dtype=float); b_ub = -np.ones(len(gens))
    res = linprog(c=np.zeros(d), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)] * d,
                  method="highs")
    if res.status == 0:
        return True, {parts[i]: res.x[i] for i in range(d) if abs(res.x[i]) > 1e-9}
    return False, None


def separating_fast(A, support, gens):
    """Find f supported on `support` with f(s_A)=1, f(s_B) <= 0 for B in gens (B != A). Float."""
    A = canon(A)
    parts = list(support); d = len(parts)
    def rv(f): return [f.get(p, 0) for p in parts]
    A_eq = np.array([rv(schur_product(A))], dtype=float); b_eq = np.array([1.0])
    rows = [rv(schur_product(B)) for B in gens if canon(B) != A]
    A_ub = np.array(rows, dtype=float); b_ub = np.zeros(len(rows))
    res = linprog(c=np.zeros(d), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)] * d, method="highs")
    if res.status == 0:
        return {parts[i]: res.x[i] for i in range(d) if abs(res.x[i]) > 1e-9}
    return None


def is_extreme_fast(A, gens, N):
    """s_A in cone of other gens? Float LP feasibility. Returns True if (numerically) extreme."""
    A = canon(A)
    parts = list(partitions(N)); idx = {p: i for i, p in enumerate(parts)}
    others = [B for B in gens if canon(B) != A]
    M = np.array([vec(schur_product(B), idx) for B in others], dtype=float).T
    t = np.array(vec(schur_product(A), idx), dtype=float)
    res = linprog(c=np.zeros(len(others)), A_eq=M, b_eq=t, bounds=[(0, None)] * len(others), method="highs")
    return res.status != 0
