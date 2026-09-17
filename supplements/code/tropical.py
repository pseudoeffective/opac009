"""Tropical criterion on the level identity P_A = sum c_B P_B (full-range I).
For a weight w on T (injective), E_w(B) = sum_{arcs} max(w_b, w_{a+1}), sigma_w(B) = prod sign(w_b - w_{a+1}).
Criterion: (1) for each level E > E_w(A), all level-competitors B with E_w(B)=E have equal sigma;
           (2) all B with E_w(B) = E_w(A) have sigma = -sigma_w(A)."""
import sys, itertools, time
from collections import defaultdict
from schurcone import *
from epsjet import inside_pairs

def E_sigma(M, w):
    E = 0; s = 1
    for (a, b) in M:
        E += max(w[b], w[a + 1]); s *= 1 if w[b] > w[a + 1] else -1
    return E, s

def tropical_test(A, comps, exhaustive_limit=9):
    lam = phi(A); b0, a0 = lam[-1], lam[0]
    AI = inside_pairs(A, b0, a0); n = len(AI)
    Ms = set(tuple(sorted(inside_pairs(B, b0, a0), reverse=True)) for B in comps if len(inside_pairs(B, b0, a0)) == n)
    T = list(range(b0, a0 + 2))
    if len(T) > exhaustive_limit: return None
    for perm in itertools.permutations(range(len(T))):
        w = {T[i]: perm[i] for i in range(len(T))}
        EA, sA = E_sigma(AI, w)
        levels = defaultdict(set)
        for M in Ms:
            E, s = E_sigma(M, w); levels[E].add(s)
        ok = all(len(v) == 1 for E, v in levels.items() if E > EA) and all(s == -sA for s in levels.get(EA, set()))
        if ok: return w
    return False

if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(4, Nmax + 1):
        gens = nested_multisets(N); t = time.time(); fails = []; skipped = 0
        for A in gens:
            if len(A) < 2: continue
            lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
            r = tropical_test(A, comps)
            if r is None: skipped += 1
            elif r is False: fails.append(A)
        print(f"N={N}: {len(gens)} nested, {len(fails)} fail tropical criterion, {skipped} skipped (T too big) [{time.time()-t:.0f}s]", flush=True)
        for A in fails[:10]: print("      ", A)
