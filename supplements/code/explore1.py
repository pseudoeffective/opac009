from schurcone import *
from fractions import Fraction

def show(f):
    return " + ".join(f"{c}·s{m}" for m, c in sorted(f.items(), key=lambda kv: kv[0], reverse=True))

A = ((2,1),(2,1),(2,1))
lam = phi(A)
N = sum(lam)
gens = nested_multisets(N)
print("lam", lam, "lam+", lam_plus(lam), "lam++", lam_plusplus(lam))
print("SSP_lam:", nested_with_phi(lam))
for rho in [lam_plus(lam), lam_plusplus(lam)]:
    I = interval(lam, rho)
    print("interval", rho, ":", I)
    Bs = [B for B in gens if phi(B) in I]
    for B in Bs:
        f = schur_product(B)
        print("  ", B, "->", {m: f[m] for m in I if m in f})
    sep = separating_functional(A, support=I, generators=Bs, N=N)
    print("  separator:", sep)
