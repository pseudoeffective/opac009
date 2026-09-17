import sys
from schurcone import *
from vertices import certificate_polyhedron, pretty
A = canon(eval(sys.argv[1])); assert is_nested(A)
which = sys.argv[2] if len(sys.argv) > 2 else "lam+"
lam = phi(A); N = sum(lam)
rho = lam_plus(lam) if which == "lam+" else (lam_plusplus(lam) if which == "lam++" else eval(which))
I = interval(lam, rho); I.sort(reverse=True)
gens = [B for B in nested_multisets(N) if phi(B) in I]
print(f"A = {A}, lam = {lam}, SSP_lam\\A = {[B for B in gens if phi(B)==lam and B!=A]}")
print(f"interval [lam,{rho}]: {len(I)} partitions, {len(gens)} nested B")
if "-t" in sys.argv:
    hdr = "".join(f"{''.join(map(str,m)):>9}" for m in I)
    print(f"{'':>40}{hdr}")
    for B in gens:
        f = schur_product(B)
        print(f"{str(B):>40}" + "".join(f"{f.get(m,0):>9}" for m in I))
verts, rays = certificate_polyhedron(A, I, gens)
print(len(verts), "vertices,", len(rays), "rays")
for v in verts: print("  V:", pretty(v))
