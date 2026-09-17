"""Sign criterion via initial forms.

A total order on indices 1..N+1 is given by a list `rank` (rank[i] = position of i; smaller = earlier, "d_i smaller").
Pair (a,b) (a>=b>=1) is positive iff rank[b] < rank[a+1]   (initial term h_a h_b)
                     negative iff rank[b] > rank[a+1]   (initial term -h_{a+1} h_{b-1}).
nu*(B) = phi(B) with negative pairs pushed; n_-(B) = number of negative pairs.

Criterion: (1') for every nu != nu*(A): all competitors B with nu*(B)=nu have equal parity of n_-;
           (2) every competitor with nu*(B)=nu*(A) has parity != parity(A).
Competitors: nested B != A, weight N, phi(B) ⊵ phi(A).
"""
import sys
from collections import defaultdict
from schurcone import *


def nustar(B, rank):
    parts = []
    neg = 0
    for r in B:
        if len(r) == 1:
            parts.append(r[0])
        else:
            a, b = r
            if rank[b] > rank[a + 1]:
                neg += 1
                parts.append(a + 1)
                if b - 1 > 0:
                    parts.append(b - 1)
            else:
                parts.extend([a, b])
    return tuple(sorted(parts, reverse=True)), neg


def check(A, rank, competitors, verbose=False):
    nuA, nA = nustar(A, rank)
    classes = defaultdict(list)
    for B in competitors:
        nu, n = nustar(B, rank)
        classes[nu].append((B, n % 2))
    ok = True
    problems = []
    for nu, lst in classes.items():
        pars = {p for _, p in lst}
        if nu == nuA:
            if any(p == nA % 2 for _, p in lst):
                ok = False
                problems.append(("(2)", nu, [B for B, p in lst if p == nA % 2]))
        elif len(pars) > 1:
            ok = False
            problems.append(("(1')", nu, lst))
    if verbose:
        print(f"  nu*(A) = {nuA}, parity(A) = {nA % 2}; {len(classes)} classes")
        for pr in problems[:6]:
            print("   PROBLEM", pr)
    return ok, problems


def block_reversal_rank(M, p, q):
    """natural order 1..M except the block [p,q] is reversed."""
    order = list(range(1, p)) + list(range(q, p - 1, -1)) + list(range(q + 1, M + 1))
    rank = {v: i for i, v in enumerate(order)}
    return rank


def rank_from_order(order):
    return {v: i for i, v in enumerate(order)}


if __name__ == "__main__":
    A = canon(eval(sys.argv[1]))
    lam = phi(A); N = sum(lam); M = N + 1
    comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
    print(f"A = {A}, lam = {lam}, {len(comps)} competitors")
    found = []
    for p in range(1, M + 1):
        for q in range(p, M + 1):
            rank = block_reversal_rank(M, p, q)
            ok, pr = check(A, rank, comps)
            if ok:
                found.append((p, q))
    print("block reversals [p,q] that work:", found)
    if len(sys.argv) > 2:
        p, q = eval(sys.argv[2])
        check(A, block_reversal_rank(M, p, q), comps, verbose=True)


def check_levels(A, rank, competitors, dvals=None, verbose=False):
    """Level-aware criterion: with d_k = dvals[k] (default rank[k]+1), w_p = sum_{k<=p} d_k, w(nu) = sum w_p.
    (1) for nu with w(nu) < w(nu*(A)): classes parity-homogeneous; (2) class of nu*(A) has parity != parity(A)."""
    M = max(rank.keys())
    if dvals is None:
        dvals = {k: rank[k] + 1 for k in rank}
    w = {0: 0}
    for k in range(1, M + 1):
        w[k] = w[k - 1] + dvals[k]
    def wt(nu): return sum(w[p] for p in nu)
    nuA, nA = nustar(A, rank)
    wA = wt(nuA)
    classes = defaultdict(list)
    for B in competitors:
        nu, n = nustar(B, rank)
        classes[nu].append((B, n % 2))
    problems = []
    for nu, lst in classes.items():
        if nu == nuA:
            bad = [B for B, p in lst if p == nA % 2]
            if bad: problems.append(("(2)", nu, bad))
        elif wt(nu) < wA:
            if len({p for _, p in lst}) > 1:
                problems.append(("(1)", nu, lst))
    if verbose:
        print(f"  nu*(A) = {nuA}, parity(A) = {nA % 2}, w(nu*A) = {wA}")
        for pr in problems[:8]: print("   PROBLEM", pr)
    return (not problems), problems
