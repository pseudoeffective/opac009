"""Hybrid: restrict competitors to R = {B nested : for every pair (a,b) in A, B has a pair (x,y) with b<=y<=x<=a},
then test for one-box / [lam,lam+] certificates."""
import sys, time
from schurcone import *
from vertices import certificate_polyhedron, pretty
from formulaF import one_box_moves

def in_R(B, A):
    for rho in A:
        if len(rho) < 2: continue
        a, b = rho
        if not any(len(s) == 2 and b <= s[1] <= s[0] <= a for s in B):
            return False
    return True

def hybrid(A, which="lam+"):
    lam = phi(A); N = sum(lam)
    if which == "onebox":
        sup = sorted(set([lam] + [mu for (a,b,i,j,mu) in one_box_moves(lam)]), reverse=True)
    elif which == "lam+":
        sup = interval(lam, lam_plus(lam))
    else:
        sup = interval(lam, lam_plusplus(lam))
    gens = [B for B in nested_multisets(N) if phi(B) in sup and in_R(B, A)]
    sepR = separating_functional(A, support=sup, generators=gens, N=N)
    return sepR, gens

if __name__ == "__main__":
    Nmax = int(sys.argv[1]); which = sys.argv[2] if len(sys.argv) > 2 else "lam+"
    for N in range(4, Nmax+1):
        gens = nested_multisets(N)
        t = time.time(); bad = []
        for A in gens:
            if len(A) < 2: continue
            sep, R = hybrid(A, which)
            if sep is None: bad.append(A)
        print(f"N={N}: {len(bad)} nested A without hybrid certificate on {which} [{time.time()-t:.0f}s]", flush=True)
        for A in bad: print("     ", A)
