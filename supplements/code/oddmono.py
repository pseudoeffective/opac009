"""Combinatorial Lemma check (even and odd cases).

For nested A with lam=phi(A), l=lam_n, u=lam_1, I=[l,u], m = #pairs of A:
  level competitors  L(A) = nested B != A, phi(B) ⊵ lam, n_I(B) = m
  restricted         L*(A) = {B in L(A): n_{I'}(B) >= n_{I'}(A) for all I' ⊆ I}
Even case (no singleton): monomial mu = prod_{(a,b) in A} delta_b, coefficient 1 in P_A;
    claim: coefficient 0 in P_B for every B in L*(A).
Odd case (singleton c): A = A_dn ⊎ A_eq ⊎ A_up (tops < c / pairs (c,c) / bottoms > c),
    claim 1: every B in L*(A) has M_eq = A_eq and pairs split the same way;
    claim 2: mu = prod_{A_dn} delta_{a+1} * prod_{A_up} delta_b has coefficient 0 in P_{M_dn} P_{M_up}
             for every B in L*(A)  (coefficient (-1)^{m_dn} in P_{A_dn} P_{A_up}).
"""
import sys, time
from schurcone import *
from epsjet import inside_pairs
from levelid import poly_of
from lstar import Lstar


def split(M, c):
    dn = [(a, b) for (a, b) in M if a < c]
    eq = [(a, b) for (a, b) in M if a == b == c]
    up = [(a, b) for (a, b) in M if b > c]
    return dn, eq, up


def check_A(A, gens):
    A = canon(A); lam = phi(A); N = sum(lam); l, u = lam[-1], lam[0]; T = list(range(l, u + 2))
    comps = [B for B in gens if B != A and dominates(phi(B), lam)]
    Ls = Lstar(A, comps)
    pairsA = [r for r in A if len(r) == 2]; sing = [r[0] for r in A if len(r) == 1]
    if not sing:
        mu = tuple(sorted(b for (a, b) in pairsA))
        assert poly_of(pairsA, T).get(mu, 0) == 1
        bad = [B for B in Ls if poly_of([r for r in B if len(r) == 2], T).get(mu, 0) != 0]
        return bad, len(Ls)
    c = sing[0]
    Adn, Aeq, Aup = split(pairsA, c)
    assert len(Adn) + len(Aeq) + len(Aup) == len(pairsA)
    mu = tuple(sorted([a + 1 for (a, b) in Adn] + [b for (a, b) in Aup]))
    assert poly_of(Adn + Aup, T).get(mu, 0) == (-1) ** len(Adn)
    bad = []
    for B in Ls:
        M = [r for r in B if len(r) == 2]
        Mdn, Meq, Mup = split(M, c)
        if not (len(Mdn) == len(Adn) and len(Meq) == len(Aeq) and len(Mup) == len(Aup)):
            bad.append(("split", B)); continue
        if poly_of(Mdn + Mup, T).get(mu, 0) != 0:
            bad.append(("coef", B))
    return bad, len(Ls)


if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(2, Nmax + 1):
        gens = nested_multisets(N); t = time.time()
        tot = 0; allbad = []; maxL = 0
        for A in gens:
            bad, nL = check_A(A, gens); tot += 1; maxL = max(maxL, nL)
            if bad: allbad.append((A, bad))
        print(f"N={N}: {tot} nested A, max |L*|={maxL}, failures {len(allbad)} [{time.time()-t:.0f}s]", flush=True)
        for A, bad in allbad[:5]: print("     ", A, bad[:3])
