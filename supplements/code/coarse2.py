import sys, time
from schurcone import *
from coarse import coarse_test
from tropical import tropical_test
Nmax = int(sys.argv[1])
for N in range(6, Nmax + 1):
    gens = nested_multisets(N); t = time.time(); fails_c = []; fails_t = []; cnt = 0; sk = 0
    for A in gens:
        if len(A) < 2 or any(len(r) == 1 for r in A): continue
        cnt += 1
        lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
        rc = coarse_test(A, comps, limit=8)
        rt = tropical_test(A, comps, exhaustive_limit=8)
        if rc is None or rt is None: sk += 1
        if rc is False: fails_c.append(A)
        if rt is False: fails_t.append(A)
    print(f"N={N}: {cnt} singleton-free A; coarse fails {len(fails_c)}, tropical fails {len(fails_t)}, skipped {sk} [{time.time()-t:.0f}s]", flush=True)
    for A in fails_c[:5]: print("     coarse fail:", A)
    for A in fails_t[:5]: print("     tropical fail:", A)
