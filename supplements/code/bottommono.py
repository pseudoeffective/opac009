"""Test the bottom-monomial criterion: coefficient of prod_{(a,b) in A} delta_b in P_M <= 0 for all level competitors M,
I = [lam_n, lam_1]."""
import sys, time
from schurcone import *
from epsjet import inside_pairs
from levelid import poly_of

def bottom_test(A, comps):
    lam = phi(A); b0, a0 = lam[-1], lam[0]
    AI = inside_pairs(A, b0, a0); n = len(AI)
    T = list(range(b0, a0 + 2))
    mono = tuple(sorted(b for (a, b) in AI))
    PA = poly_of(AI, T)
    assert PA.get(mono, 0) == 1
    bad = []
    for B in comps:
        BI = inside_pairs(B, b0, a0)
        if len(BI) != n: continue
        c = poly_of(BI, T).get(mono, 0)
        if c > 0: bad.append((tuple(sorted(BI, reverse=True)), c))
    return bad

if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(4, Nmax + 1):
        gens = nested_multisets(N); t = time.time(); fails = []
        for A in gens:
            if len(A) < 2: continue
            lam = phi(A)
            comps = [B for B in gens if B != A and dominates(phi(B), lam)]
            bad = bottom_test(A, comps)
            if bad: fails.append((A, bad))
        print(f"N={N}: {len(gens)} nested, {len(fails)} fail bottom-monomial criterion [{time.time()-t:.0f}s]", flush=True)
        for A, bad in fails[:8]: print("      ", A, "->", bad[:3])
