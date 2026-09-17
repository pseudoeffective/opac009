import sys, time
from schurcone import *
from avoid3 import restricted_avoidance

def distinct(A):
    l = phi(A); return len(set(l)) == len(l)

Nmax = int(sys.argv[1])
good = {}   # A -> reason
for N in range(1, Nmax+1):
    gens = nested_multisets(N)
    hard = []
    for A in gens:
        if distinct(A) or len(A) == 1:
            good[A] = "White"; continue
        reason = None
        for rho in set(A):
            rest = tuple(x for x in A)
            rest = list(A); rest.remove(rho); rest = canon(rest)
            if rest not in good: continue
            if len(rho) == 2 and rho[0] == rho[1]:
                reason = f"Thm1 remove {rho}"; break
            if restricted_avoidance(A, rho):
                reason = f"avoid remove {rho}"; break
        if reason: good[A] = reason
        else: hard.append(A)
    print(f"N={N}: {len(gens)} nested, {len(hard)} not provable by closure:", flush=True)
    for A in hard: print("     ", A, " phi =", phi(A))
