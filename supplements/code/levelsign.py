"""Point-evaluation version: find I and a linear order on T with sign(P_{A_I}) = + and sign(P_{B_I}) = - for all
level competitors (n_I(B) = n_I(A))."""
import sys, itertools, time
from schurcone import *
from epsjet import inside_pairs, sign_of

def sign_test(A, b0, a0, comps, exhaustive_limit=9):
    AI = inside_pairs(A, b0, a0); n = len(AI)
    Ms = set()
    for B in comps:
        BI = inside_pairs(B, b0, a0)
        if len(BI) == n: Ms.add(tuple(sorted(BI, reverse=True)))
    T = list(range(b0, a0 + 2))
    if len(T) > exhaustive_limit: return None
    for perm in itertools.permutations(T):
        e = {k: i for i, k in enumerate(perm)}
        if sign_of(AI, e) != 1: continue
        if all(sign_of(list(M), e) == -1 for M in Ms):
            return perm
    return False

if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(4, Nmax + 1):
        gens = nested_multisets(N); t = time.time(); unproved = []; skipped = 0
        for A in gens:
            if len(A) < 2: continue
            lam = phi(A)
            comps = [B for B in gens if B != A and dominates(phi(B), lam)]
            ok = None
            vals = sorted(set(lam))
            for b0 in vals:
                for a0 in vals:
                    if a0 < b0: continue
                    r = sign_test(A, b0, a0, comps)
                    if r: ok = (b0, a0, r); break
                if ok: break
            if not ok: unproved.append(A)
        print(f"N={N}: {len(gens)} nested, {len(unproved)} not proved by sign-point criterion [{time.time()-t:.0f}s]", flush=True)
        for A in unproved[:10]: print("      ", A)
