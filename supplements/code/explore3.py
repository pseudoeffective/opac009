import sys
from schurcone import *
from vertices import certificate_polyhedron, pretty

def study(A, verbose=False):
    A = canon(A)
    assert is_nested(A), A
    lam = phi(A); N = sum(lam)
    gens = nested_multisets(N)
    out = [f"A = {A}, lam = {lam}, SSP_lam = {[B for B in nested_with_phi(lam) if B != A]}"]
    for name, rho in [("lam+", lam_plus(lam)), ("lam++", lam_plusplus(lam))]:
        I = interval(lam, rho); I.sort(reverse=True)
        Bs = [B for B in gens if phi(B) in I]
        verts, rays = certificate_polyhedron(A, I, Bs)
        out.append(f"  on [lam,{name}] ({len(I)} parts, {len(Bs)} B): {len(verts)} vertices, {len(rays)} rays")
        for v in verts: out.append("     V: " + pretty(v))
        if verts: break
    print("\n".join(out), flush=True)

for A in [((4,1),(2,2)), ((3,3),(2,1)), ((5,1),(3,3)), ((5,2),(3,3),(1,)), ((4,2),(3,3),(1,)),
          ((3,2),(3,1)), ((4,3),(4,1),(5,)), ((3,1),(2,1)), ((4,2),(3,2)), ((5,1),(4,3),(4,2))]:
    study(A)
