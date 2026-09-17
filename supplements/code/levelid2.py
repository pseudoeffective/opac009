import sys, time
from schurcone import *
from levelid import level_test
Nmax = int(sys.argv[1]); Nmin = int(sys.argv[2]) if len(sys.argv) > 2 else 4
for N in range(Nmin, Nmax + 1):
    gens = nested_multisets(N); t = time.time(); unproved_full = []; unproved_any = []
    for A in gens:
        if len(A) < 2: continue
        lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
        good, why, Ms = level_test(A, lam[-1], lam[0], comps)
        if not good:
            unproved_full.append(A)
            vals = sorted(set(lam)); ok = False
            for b0 in vals:
                for a0 in vals:
                    if a0 < b0: continue
                    g2, _, _ = level_test(A, b0, a0, comps)
                    if g2: ok = True; break
                if ok: break
            if not ok: unproved_any.append(A)
    print(f"N={N}: {len(gens)} nested; full-range level test fails for {len(unproved_full)}; no interval works for {len(unproved_any)} [{time.time()-t:.0f}s]", flush=True)
    for A in unproved_full[:6]: print("     full-range fail:", A)
    for A in unproved_any[:6]: print("     ALL fail:", A)
