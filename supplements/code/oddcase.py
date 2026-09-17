import sys
from schurcone import *
from lstar import Lstar
from epsjet import inside_pairs
from levelid import poly_of
from leveldual import dual
A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam)
comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
b0, a0 = lam[-1], lam[0]; T = list(range(b0, a0+2))
Ls = Lstar(A, comps)
print("A =", A, "P_A =", poly_of(inside_pairs(A,b0,a0), T))
print("L*(A):")
for B in Ls:
    print("   ", B, "  P =", poly_of(inside_pairs(B,b0,a0), T))
