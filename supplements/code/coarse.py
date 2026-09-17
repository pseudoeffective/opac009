"""Coarse-point criterion: ordered set partitions of T (classes with values c_1 > c_2 > ... ),
P evaluated: arc (b -> a+1) contributes 0 if same class, +1 if class(b) earlier(larger value), -1 otherwise.
Need P_A != 0 and every level competitor has P_B = 0 or sign(P_B) = -sign(P_A)."""
import sys, itertools, time
from schurcone import *
from epsjet import inside_pairs

def ordered_set_partitions(elems):
    elems = list(elems)
    n = len(elems)
    # assign each element a block index 0..k-1 with all blocks nonempty; order of blocks = value order
    for k in range(1, n + 1):
        for assign in itertools.product(range(k), repeat=n):
            if set(assign) == set(range(k)):
                yield dict(zip(elems, assign))

def sgn(M, cls):
    s = 1
    for (a, b) in M:
        cb, ch = cls[b], cls[a + 1]
        if cb == ch: return 0
        s *= 1 if cb < ch else -1   # smaller block index = larger value
    return s

def coarse_test(A, comps, limit=8):
    lam = phi(A); b0, a0 = lam[-1], lam[0]
    AI = inside_pairs(A, b0, a0); n = len(AI)
    Ms = set(tuple(sorted(inside_pairs(B, b0, a0), reverse=True)) for B in comps if len(inside_pairs(B, b0, a0)) == n)
    T = list(range(b0, a0 + 2))
    if len(T) > limit: return None
    for cls in ordered_set_partitions(T):
        sA = sgn(AI, cls)
        if sA == 0: continue
        if all(sgn(M, cls) in (0, -sA) for M in Ms):
            return cls
    return False

if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(4, Nmax + 1):
        gens = nested_multisets(N); t = time.time(); fails = []; skipped = 0
        for A in gens:
            if len(A) < 2: continue
            lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
            r = coarse_test(A, comps)
            if r is None: skipped += 1
            elif r is False: fails.append(A)
        print(f"N={N}: {len(gens)} nested, {len(fails)} fail coarse criterion, {skipped} skipped [{time.time()-t:.0f}s]", flush=True)
        for A in fails[:10]: print("      ", A)
