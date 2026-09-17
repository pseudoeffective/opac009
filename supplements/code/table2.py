import sys
from schurcone import *
from formulaF import first_last
def move(lam, a, b):
    first, last = first_last(lam); l = list(lam); l[first[a]-1] += 1; l[last[b]-1] -= 1
    return tuple(x for x in l if x > 0)
A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam)
I = interval(lam, lam_plus(lam)); I.sort(reverse=True)
gens = [B for B in nested_multisets(N) if phi(B) in I]
print("A =", A, "lam =", lam, "interval:", I)
print(f"{'':>28}" + "".join(f"{str(m):>16}" for m in I))
for B in gens:
    f = schur_product(B)
    print(f"{str(B):>28}" + "".join(f"{f.get(m,0):>16}" for m in I))
