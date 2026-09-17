import sys
from schurcone import *
from hybrid import in_R

def distinct(A):
    l = phi(A); return len(set(l)) == len(l)

def rigid(A):
    lam = phi(A)
    return all(B == A for B in nested_with_phi(lam) if in_R(B, A))

def strip_pp(A):
    return canon([r for r in A if not (len(r) == 2 and r[0] == r[1])])

Nmax = int(sys.argv[1])
for N in range(2, Nmax+1):
    gens = nested_multisets(N)
    white = [A for A in gens if distinct(A)]
    G = [A for A in gens if distinct(strip_pp(A))]
    rig = [A for A in gens if rigid(A)]
    rigG = set(rig) | set(G)
    # closure under Thm1: A extreme if strip_pp(A) in (white ∪ rigid)  -- rigid is proven for any A
    cov = [A for A in gens if strip_pp(A) in set(white) or rigid(strip_pp(A)) or rigid(A)]
    print(f"N={N}: nested {len(gens)}, White {len(white)}, G(=White+Thm1) {len(G)}, rigid {len(rig)}, "
          f"G∪rigid {len(rigG)}, covered(closure) {len(cov)}, uncovered {len(gens)-len(cov)}")
    if N <= 10:
        print("    uncovered:", [A for A in gens if A not in set(cov)])
