import sys
from schurcone import *
from epsjet import inside_pairs
from levelid import poly_of
A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam)
comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
b0, a0 = lam[-1], lam[0]
AI = inside_pairs(A, b0, a0); n = len(AI); T = list(range(b0, a0 + 2))
Ms = sorted(set(tuple(sorted(inside_pairs(B, b0, a0), reverse=True)) for B in comps if len(inside_pairs(B, b0, a0)) == n))
PA = poly_of(AI, T); PMs = {M: poly_of(M, T) for M in Ms}
print("A =", A, " arcs:", sorted((b, a+1) for (a,b) in AI))
print("level competitors:", Ms)
print("P_A =", PA)
good = []
for mono, cA in PA.items():
    s = 1 if cA > 0 else -1
    if all(s * p.get(mono, 0) <= 0 for p in PMs.values()):
        good.append((mono, cA, [p.get(mono,0) for p in PMs.values()]))
print("working monomials:")
for g in good: print("   ", g)

print("\ncoefficients of A's monomials across competitors:")
for mono, cA in PA.items():
    row = [(M, p.get(mono,0)) for M, p in PMs.items() if p.get(mono,0) != 0]
    print(f"  {mono} (A: {cA:+d}):", row)
