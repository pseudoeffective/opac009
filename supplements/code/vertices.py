"""Enumerate all vertices of the polyhedron of separating functionals supported on a given interval,
normalised by <f, s_A> = 1.  (cdd exact.)"""
from schurcone import *
from fractions import Fraction
import cdd.gmp as cdd

def certificate_polyhedron(A, support, gens):
    A = canon(A)
    parts = list(support)
    d = len(parts)
    def rv(f): return [Fraction(f.get(p,0)) for p in parts]
    rows = [[Fraction(-1)] + rv(schur_product(A))]
    for B in gens:
        if canon(B) == A: continue
        rows.append([Fraction(0)] + [-x for x in rv(schur_product(B))])
    mat = cdd.matrix_from_array(rows, lin_set=[0], rep_type=cdd.RepType.INEQUALITY)
    poly = cdd.polyhedron_from_matrix(mat)
    gen = cdd.copy_generators(poly)
    verts, rays = [], []
    for i, row in enumerate(gen.array):
        v = {parts[j]: row[j+1] for j in range(d) if row[j+1] != 0}
        if row[0] == 1: verts.append(v)
        else: rays.append(v)
    return verts, rays

def pretty(v):
    return " ".join(f"{'+' if c>0 else '-'}{abs(c) if abs(c)!=1 else ''}s{''.join(map(str,m))}" for m, c in sorted(v.items(), key=lambda kv: kv[0], reverse=True))

if __name__ == "__main__":
    import sys
    A = ((2,1),(2,1),(2,1))
    lam = phi(A); N = sum(lam)
    gens = nested_multisets(N)
    I = interval(lam, lam_plusplus(lam))
    Bs = [B for B in gens if phi(B) in I]
    verts, rays = certificate_polyhedron(A, I, Bs)
    print(len(verts), "vertices,", len(rays), "rays")
    for v in verts: print("  V:", pretty(v))
    for r in rays: print("  R:", pretty(r))
