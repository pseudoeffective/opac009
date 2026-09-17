import sys, time
from schurcone import *
A = eval(sys.argv[1])
lam = phi(A); N = sum(lam)
gens = nested_multisets(N)
above = [B for B in gens if dominates(phi(B), lam)]
t=time.time()
ext, wit = is_extreme(A, above, N)
print("extreme:", ext, f"{time.time()-t:.1f}s")
if not ext:
    print("witness: s_A = sum of")
    for B, c in wit.items():
        print("   ", c, "*", B)
    # verify exactly
    tot = {}
    for B, c in wit.items():
        for m, v in schur_product(B).items():
            tot[m] = tot.get(m, 0) + c * v
    tot = {m: v for m, v in tot.items() if v != 0}
    print("verified:", tot == {m: Fraction(v) for m, v in schur_product(A).items()})
