import sys, time
from schurcone import *
from vertices import certificate_polyhedron, pretty
from hybrid import in_R
from hybrid_verts import movename
A = canon(eval(sys.argv[1])); assert is_nested(A)
lam = phi(A); N = sum(lam)
sup = interval(lam, lam_plusplus(lam)); sup.sort(reverse=True)
t=time.time()
gens = [B for B in nested_multisets(N) if phi(B) in sup and in_R(B, A)]
print(f"A = {A}, lam = {lam}; [lam,lam++] has {len(sup)} shapes; R-competitors: {len(gens)}  ({time.time()-t:.0f}s)")
sep = separating_functional(A, support=sup, generators=gens, N=N)
print("one certificate:", None if sep is None else " ".join(f"{'+' if c>0 else '-'}{abs(c) if abs(c)!=1 else ''}s{movename(lam,m) if m!=lam else 'λ'}" for m,c in sorted(sep.items(), key=lambda kv: kv[0], reverse=True)))
if "-v" in sys.argv:
    verts, rays = certificate_polyhedron(A, sup, gens)
    print(len(verts), "vertices,", len(rays), "rays")
    for v in verts[:15]:
        print("  V:", " ".join(f"{'+' if c>0 else '-'}{abs(c) if abs(c)!=1 else ''}s{movename(lam,m) if m!=lam else 'λ'}" for m,c in sorted(v.items(), key=lambda kv: kv[0], reverse=True)))
