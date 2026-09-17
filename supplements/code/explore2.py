import sys
from schurcone import *
from vertices import certificate_polyhedron, pretty

def study(A, rho=None, show_table=True):
    lam = phi(A); N = sum(lam)
    gens = nested_multisets(N)
    if rho is None: rho = lam_plusplus(lam)
    I = interval(lam, rho)
    I.sort(reverse=True)
    Bs = [B for B in gens if phi(B) in I]
    print(f"A = {A}, lam = {lam}, interval up to {rho}: {len(I)} partitions, {len(Bs)} nested B")
    if show_table:
        hdr = "".join(f"{''.join(map(str,m)):>9}" for m in I)
        print(f"{'':>34}{hdr}")
        for B in Bs:
            f = schur_product(B)
            row = "".join(f"{f.get(m,0):>9}" for m in I)
            print(f"{str(B):>34}{row}")
    verts, rays = certificate_polyhedron(A, I, Bs)
    print(len(verts), "vertices,", len(rays), "rays")
    for v in verts: print("  V:", pretty(v))
    return verts, rays

if __name__ == "__main__":
    A = eval(sys.argv[1])
    rho = eval(sys.argv[2]) if len(sys.argv) > 2 else None
    study(A, rho)
