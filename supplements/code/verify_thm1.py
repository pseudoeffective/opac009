"""Numerical sanity check of the specialisation in Theorem 1 and of the two-tops/two-bottoms certificates."""
from fractions import Fraction
from schurcone import *
import math

def make_g(p, K=12):
    # delta(t) = -t for t<=p, -(t-1) for t>=p+1  (strictly decreasing except delta(p)=delta(p+1))
    u = {0: 0}
    for t in range(1, K+2):
        d = -t if t <= p else -(t-1)
        u[t] = u[t-1] + d
    return lambda k: math.exp(u[k]) if k >= 0 else 0.0

for p in range(1, 6):
    g = make_g(p)
    bad = []
    for a in range(1, 11):
        for b in range(1, a+1):
            G = g(a)*g(b) - g(a+1)*g(b-1)
            if (a, b) == (p, p):
                if abs(G) > 1e-12: bad.append(((a,b), G))
            elif G <= 0: bad.append(((a,b), G))
    print(f"p={p}: specialisation OK: {not bad}", bad[:3])

# the interval lemma for I=[b,a]
def make_gI(b, a, K=12):
    u = {0: 0}
    for t in range(1, K+2):
        d = -t if t <= b else (-b if t <= a+1 else -(t - (a+1-b)))
        u[t] = u[t-1] + d
    return lambda k: math.exp(u[k]) if k >= 0 else 0.0
for (b, a) in [(1,3), (2,5), (3,3), (1,1), (2,2)]:
    g = make_gI(b, a); bad = []
    for x in range(1, 11):
        for y in range(1, x+1):
            G = g(x)*g(y) - g(x+1)*g(y-1)
            inside = b <= y <= x <= a
            if inside and abs(G) > 1e-12: bad.append(((x,y),G))
            if not inside and G <= 0: bad.append(((x,y),G))
    print(f"I=[{b},{a}]: interval lemma OK: {not bad}", bad[:3])

# certificates for two-tops / two-bottoms
from formulaF import first_last
def move(lam, a, b):
    first, last = first_last(lam); l = list(lam); l[first[a]-1] += 1; l[last[b]-1] -= 1
    return tuple(x for x in l if x > 0)
def check(A, f):
    lam = phi(A); N = sum(lam); I = interval(lam, lam_plus(lam))
    gens = [B for B in nested_multisets(N) if phi(B) in I]
    vals = {B: sum(c*schur_product(B).get(m,0) for m,c in f.items()) for B in gens}
    ok = vals[A] > 0 and all(v <= 0 for B,v in vals.items() if B != A)
    return ok, vals[A], max(v for B,v in vals.items() if B != A)
for (v,x,y) in [(3,2,1),(5,3,1),(6,5,2),(4,3,2),(7,2,1)]:
    A = canon([(v,x),(v,y)]); lam = phi(A)
    f = {lam:1, move(lam,v,v):1, move(lam,x,y):1, move(lam,v,y):-1}
    print("two tops", A, check(A, f))
for (x,y,v) in [(3,2,1),(5,3,1),(6,5,2),(4,3,2),(7,6,1),(4,2,1)]:
    A = canon([(x,v),(y,v)]); lam = phi(A)
    f = {lam:1, move(lam,v,v):1, move(lam,x,y):1, move(lam,x,v):-1}
    print("two bottoms", A, check(A, f))
for (j,i) in [(2,1),(3,1),(3,2),(5,2),(6,1),(4,3)]:
    A = canon([(j,i),(j,i)]); lam = phi(A)
    f = {lam:1, move(lam,j,j):1, move(lam,i,i):1, move(lam,j,i):-1}
    print("(j,i)^2", A, check(A, f))
