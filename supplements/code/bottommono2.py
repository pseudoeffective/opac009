import sys, time
from schurcone import *
from bottommono import bottom_test
Nmax = int(sys.argv[1])
for N in range(4, Nmax + 1):
    gens = nested_multisets(N); t = time.time(); fails = []; cnt = 0
    for A in gens:
        if len(A) < 2 or any(len(r) == 1 for r in A): continue
        cnt += 1
        lam = phi(A)
        comps = [B for B in gens if B != A and dominates(phi(B), lam)]
        bad = bottom_test(A, comps)
        if bad: fails.append((A, bad))
    print(f"N={N}: {cnt} nested singleton-free A with >=2 pairs, {len(fails)} fail bottom-monomial [{time.time()-t:.0f}s]", flush=True)
    for A, bad in fails[:8]: print("      ", A, "->", bad[:3])
