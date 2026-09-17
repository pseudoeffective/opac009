import sys, time
from schurcone import *
from dcount import d_count, assignments
Nmax = int(sys.argv[1]); which = sys.argv[2]
mism = 0; tested = 0
for N in range(2, Nmax+1):
    gens = nested_multisets(N)
    for lam in partitions(N):
        if len(lam) < 2: continue
        Bs = [B for B in gens if phi(B) == lam]
        rho = lam_plus(lam) if which == "plus" else lam_plusplus(lam)
        mus = [m for m in interval(lam, rho) if m != lam]
        for B in Bs:
            f = schur_product(B)
            asg = assignments(B, lam)
            if len(asg) > 40: asg = asg[:40]
            for mu in mus:
                c = f.get(mu, 0)
                ds = [d_count(mu, lam, [list(b) for b in bl]) for bl in asg]
                tested += 1
                if any(d != c for d in ds):
                    mism += 1
                    print(f"MISMATCH lam={lam} B={B} mu={mu}: c={c} d={sorted(set(ds))}")
    print(f"N={N} done, tested={tested}, mismatches={mism}", flush=True)
