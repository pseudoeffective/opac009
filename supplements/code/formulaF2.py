"""Test product formula: for mu in [lam, lam+] with several segments, c_B^mu = prod over segments of F."""
import sys
from schurcone import *
from dcount import assignments
from formulaF import count_F

def segments(lam, mu):
    n = max(len(lam), len(mu))
    l = list(lam) + [0]*(n-len(lam)); m = list(mu) + [0]*(n-len(mu))
    Pl = Pm = 0; segs = []; cur = None
    for k in range(1, n+1):
        Pl += l[k-1]; Pm += m[k-1]
        d = Pm - Pl
        if d not in (0,1): return None
        if d == 1 and cur is None: cur = k
        if d == 0 and cur is not None: segs.append((cur, k)); cur = None
    return segs

Nmax = int(sys.argv[1]); bad = 0; tested = 0; multi = 0
for N in range(2, Nmax+1):
    for lam in partitions(N):
        if len(lam) < 2: continue
        Bs = nested_with_phi(lam)
        for mu in interval(lam, lam_plus(lam)):
            segs = segments(lam, mu)
            if segs is None or len(segs) == 0: continue
            for B in Bs:
                c = schur_product(B).get(mu, 0)
                bl = assignments(B, lam)[0]
                pairs = [tuple(x) for x in bl if len(x) == 2]
                F = 1
                for (i, j) in segs:
                    F *= count_F(lam, pairs, i, j)
                tested += 1
                if len(segs) > 1: multi += 1
                if F != c:
                    bad += 1
                    if bad < 10: print("MISMATCH", lam, B, mu, segs, c, F)
    print(f"N={N}: tested {tested} (multi-segment {multi}), mismatches {bad}", flush=True)
