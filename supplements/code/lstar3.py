import sys, time
from schurcone import *
from lstar import Lstar, bottom_ok, top_ok
from epsjet import inside_pairs
from levelid import poly_of
Nmax = int(sys.argv[1]); Nmin = int(sys.argv[2]) if len(sys.argv)>2 else 6
for N in range(Nmin, Nmax + 1):
    gens = nested_multisets(N); t = time.time(); fb = []; cnt = 0; strict = 0
    for A in gens:
        if len(A) < 2 or any(len(r)==1 for r in A): continue
        cnt += 1
        lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
        Ls = Lstar(A, comps)
        if not bottom_ok(A, Ls):
            b0, a0 = lam[-1], lam[0]; T = list(range(b0, a0+2)); AI = inside_pairs(A, b0, a0); mono = tuple(sorted(b for (a,b) in AI))
            bad = [(B, poly_of(inside_pairs(B,b0,a0), T).get(mono,0)) for B in Ls if poly_of(inside_pairs(B,b0,a0), T).get(mono,0) > 0]
            fb.append((A, bad))
    print(f"N={N}: {cnt} singleton-free A; bottom-monomial vs L* fails {len(fb)} [{time.time()-t:.0f}s]", flush=True)
    for A, bad in fb[:6]: print("      ", A, "->", bad[:3])
