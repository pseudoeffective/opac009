import sys, itertools
from collections import defaultdict
from schurcone import *
from epsjet import inside_pairs
from tropical import E_sigma
A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam)
comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
b0, a0 = lam[-1], lam[0]; AI = inside_pairs(A, b0, a0); n = len(AI)
Ms = sorted(set(tuple(sorted(inside_pairs(B, b0, a0), reverse=True)) for B in comps if len(inside_pairs(B, b0, a0)) == n))
T = list(range(b0, a0 + 2))
print("A =", A, "arcs", sorted((b,a+1) for (a,b) in AI), "T =", T)
print("level competitors:", Ms)
good = []
for perm in itertools.permutations(range(len(T))):
    w = {T[i]: perm[i] for i in range(len(T))}
    EA, sA = E_sigma(AI, w)
    levels = defaultdict(set)
    for M in Ms:
        E, s = E_sigma(M, w); levels[E].add(s)
    if all(len(v) == 1 for E, v in levels.items() if E > EA) and all(s == -sA for s in levels.get(EA, set())):
        # describe w as the order of T by increasing w
        order = sorted(T, key=lambda k: w[k])
        good.append((tuple(order), sA, sorted(levels.items())))
print(len(good), "working weight orders (T listed in increasing w):")
for g in good[:40]: print("   ", g[0], "sigma(A)=", g[1])
