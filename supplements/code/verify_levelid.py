"""Sanity check of the Level Identity Theorem on actual decompositions.

For non-nested A (which are not extreme), find an exact decomposition s_A = sum c_B s_B over ALL generators
B != A (nested or not), then verify for every interval I=[x,y], 1<=x<=y<=N:
  (i)  n_I(B) >= n_I(A) for every B in the support;
  (ii) K_A * P_{A_I} = sum_{n_I(B)=n_I(A)} c_B K_B P_{B_I}  (as polynomials in delta_k, k in T=[x,y+1]),
       K_B = prod_{pairs outside I} G0(a,b) * prod_{pairs inside I} g(a)g(b) * prod_{singletons} g(c),
       with the base d0 strictly decreasing except constant on T.
"""
import sys, math, time
from schurcone import *
from epsjet import inside_pairs
from levelid import poly_of


def base(N, x, y, step=0.07):
    d = {}
    for k in range(1, N + 2):
        d[k] = -step * (x if x <= k <= y + 1 else k)
    return d


def gfun(d, N):
    u = {0: 0.0}
    for k in range(1, N + 2):
        u[k] = u[k - 1] + d[k]
    return {k: math.exp(u[k]) for k in range(0, N + 2)}


def K_of(B, x, y, g):
    K = 1.0
    for r in B:
        if len(r) == 1:
            K *= g[r[0]]
        else:
            a, b = r
            if x <= b <= a <= y:
                K *= g[a] * g[b]
            else:
                K *= g[a] * g[b] - g[a + 1] * g[b - 1]
    return K


def check(A, wit, N, verbose=False):
    A = canon(A)
    problems = []
    for x in range(1, N + 1):
        for y in range(x, N + 1):
            nA = len(inside_pairs(A, x, y))
            # (i)
            for B in wit:
                if len(inside_pairs(B, x, y)) < nA:
                    problems.append(("(i)", (x, y), B))
            # (ii)
            T = list(range(x, y + 2))
            g = gfun(base(N, x, y), N)
            lhs = {m: c * K_of(A, x, y, g) for m, c in poly_of(inside_pairs(A, x, y), T).items()}
            rhs = {}
            for B, cB in wit.items():
                if len(inside_pairs(B, x, y)) != nA: continue
                KB = K_of(B, x, y, g)
                for m, c in poly_of(inside_pairs(B, x, y), T).items():
                    rhs[m] = rhs.get(m, 0.0) + float(cB) * KB * c
            monos = set(lhs) | set(rhs)
            scale = max(abs(v) for v in lhs.values())
            for m in monos:
                if abs(lhs.get(m, 0.0) - rhs.get(m, 0.0)) > 1e-9 * scale:
                    problems.append(("(ii)", (x, y), m, lhs.get(m, 0.0), rhs.get(m, 0.0)))
    return problems


if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(2, Nmax + 1):
        gens = all_multisets(N)
        nonnested = [A for A in gens if not is_nested(A)]
        t = time.time(); nprob = 0; ndec = 0
        for A in nonnested:
            ext, wit = is_extreme(A, generators=gens, N=N)
            assert not ext, A
            ndec += 1
            pr = check(A, wit, N)
            if pr:
                nprob += 1
                print("   PROBLEM", A, pr[:3])
        print(f"N={N}: {len(nonnested)} non-nested A decomposed, problems {nprob} [{time.time()-t:.0f}s]", flush=True)
