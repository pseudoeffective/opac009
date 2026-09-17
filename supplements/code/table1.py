import sys
from schurcone import *
from formulaF import one_box_moves
A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam)
moves = one_box_moves(lam)
sup = [lam] + [mu for (a,b,i,j,mu) in moves]
names = {lam: "lam"}
for (a,b,i,j,mu) in moves: names[mu] = f"[{a}>{b}]"
I = interval(lam, lam_plus(lam))
gens = [B for B in nested_multisets(N) if phi(B) in sup]
print("A =", A, " lam =", lam)
print(f"{'':>36}" + "".join(f"{names[m]:>7}" for m in sup))
for B in gens:
    f = schur_product(B)
    print(f"{str(B):>36}" + "".join(f"{f.get(m,0):>7}" for m in sup))
