import sys
from schurcone import *
from vertices import certificate_polyhedron, pretty
from formulaF import one_box_moves, first_last
from hybrid import in_R

def movename(lam, mu):
    # express mu - lam as moves
    n = max(len(lam), len(mu)); l = list(lam)+[0]*(n-len(lam)); m = list(mu)+[0]*(n-len(mu))
    diff = [m[k]-l[k] for k in range(n)]
    plus = [k for k in range(n) if diff[k] == 1]; minus = [k for k in range(n) if diff[k] == -1]
    if all(d in (-1,0,1) for d in diff) and len(plus) == len(minus):
        return "".join(f"[{l[p]}>{l[q]}]" for p, q in zip(plus, minus))
    return str(mu)

A = canon(eval(sys.argv[1])); assert is_nested(A)
which = sys.argv[2] if len(sys.argv) > 2 else "onebox"
lam = phi(A); N = sum(lam)
if which == "onebox":
    sup = sorted(set([lam] + [mu for (a,b,i,j,mu) in one_box_moves(lam)]), reverse=True)
else:
    sup = interval(lam, lam_plus(lam))
gens = [B for B in nested_multisets(N) if phi(B) in sup and in_R(B, A)]
print(f"A = {A}, lam = {lam}; R-competitors: {len(gens)}; SSP_lam∩R\\A = {[B for B in gens if phi(B)==lam and B!=A]}")
verts, rays = certificate_polyhedron(A, sup, gens)
print(len(verts), "vertices,", len(rays), "rays")
def pretty2(v):
    return " ".join(f"{'+' if c>0 else '-'}{abs(c) if abs(c)!=1 else ''}s{movename(lam,m) if m!=lam else 'λ'}" for m, c in sorted(v.items(), key=lambda kv: kv[0], reverse=True))
for v in verts: print("  V:", pretty2(v))
