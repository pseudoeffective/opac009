"""For each nested A: is there an interval I and a single monomial m with coef_A(m) != 0 and
sign(coef_A(m)) * coef_M(m) <= 0 for all level competitors M?"""
import sys, time
from schurcone import *
from epsjet import inside_pairs
from levelid import poly_of

def find_mono(A, comps):
    lam = phi(A); vals = sorted(set(lam))
    for b0 in vals:
        for a0 in vals:
            if a0 < b0: continue
            AI = inside_pairs(A, b0, a0); n = len(AI)
            if n == 0: continue
            T = list(range(b0, a0 + 2))
            Ms = set(tuple(sorted(inside_pairs(B, b0, a0), reverse=True)) for B in comps if len(inside_pairs(B, b0, a0)) == n)
            if tuple(sorted(AI, reverse=True)) in Ms: continue
            PA = poly_of(AI, T); PMs = [poly_of(M, T) for M in Ms]
            for mono, cA in PA.items():
                s = 1 if cA > 0 else -1
                if all(s * p.get(mono, 0) <= 0 for p in PMs):
                    return (b0, a0, mono, cA)
    return None

if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(4, Nmax + 1):
        gens = nested_multisets(N); t = time.time(); fails = []; kinds = {}
        for A in gens:
            if len(A) < 2: continue
            lam = phi(A)
            comps = [B for B in gens if B != A and dominates(phi(B), lam)]
            r = find_mono(A, comps)
            if r is None: fails.append(A)
        print(f"N={N}: {len(gens)} nested, {len(fails)} without single-monomial certificate [{time.time()-t:.0f}s]", flush=True)
        for A in fails[:10]: print("      ", A)
