import sys, time
from schurcone import *
from lstar import Lstar, bottom_ok, top_ok
Nmax = int(sys.argv[1]); Nmin = int(sys.argv[2]) if len(sys.argv)>2 else 6
for N in range(Nmin, Nmax + 1):
    gens = nested_multisets(N); t = time.time(); fboth = []; fb = []; ft = []; cnt = 0
    for A in gens:
        if len(A) < 2: continue
        cnt += 1
        lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
        Ls = Lstar(A, comps)
        b = bottom_ok(A, Ls); tp = top_ok(A, Ls)
        if not b: fb.append(A)
        if not tp: ft.append(A)
        if not (b or tp): fboth.append(A)
    print(f"N={N}: {cnt} A; bottom fails {len(fb)}, top fails {len(ft)}, both fail {len(fboth)} [{time.time()-t:.0f}s]", flush=True)
    for A in fboth[:8]: print("      both fail:", A)
    if N <= 10:
        print("   bottom fails:", fb)
        print("   top fails:", ft)
