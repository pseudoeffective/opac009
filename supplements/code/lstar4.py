import sys
from collections import Counter
from schurcone import *
from lstar import Lstar
from epsjet import inside_pairs
from levelid import poly_of
Nmax = int(sys.argv[1])
hist = Counter(); examples = {}
for N in range(6, Nmax + 1):
    gens = nested_multisets(N)
    for A in gens:
        if len(A) < 2 or any(len(r)==1 for r in A): continue
        lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
        b0, a0 = lam[-1], lam[0]; T = list(range(b0, a0+2)); AI = inside_pairs(A, b0, a0); mono = tuple(sorted(b for (a,b) in AI))
        for B in Lstar(A, comps):
            c = poly_of(inside_pairs(B,b0,a0), T).get(mono, 0)
            # count S-terms: subsets S with the multiset condition
            BI = inside_pairs(B, b0, a0); betaB = sorted(b for (a,b) in BI)
            same_bottom = (betaB == list(mono))
            hist[(c, same_bottom)] += 1
            if c < 0 and (c, same_bottom) not in examples: examples[(c, same_bottom)] = (A, B)
            if same_bottom and (c, same_bottom) not in examples: examples[(c, same_bottom)] = (A, B)
print("histogram of (coefficient, beta(B)==beta(A)):", dict(hist))
for k, v in examples.items(): print(k, v)
