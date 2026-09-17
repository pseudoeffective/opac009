import sys
from schurcone import *
from vertices import certificate_polyhedron, pretty
from formulaF import one_box_moves
A = canon(eval(sys.argv[1])); assert is_nested(A)
lam = phi(A); N = sum(lam)
sup = [lam] + [mu for (a,b,i,j,mu) in one_box_moves(lam)]
sup = sorted(set(sup), reverse=True)
I = interval(lam, lam_plus(lam))
gens = [B for B in nested_multisets(N) if phi(B) in I]
print(f"A = {A}, lam = {lam}; one-box support: {len(sup)} shapes; {len(gens)} nested B in [lam,lam+]")
verts, rays = certificate_polyhedron(A, sup, gens)
print(len(verts), "vertices,", len(rays), "rays")
for v in verts: print("  V:", pretty(v))
