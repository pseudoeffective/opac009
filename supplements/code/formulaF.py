"""Formula (F): path-count formula for c_B^{lam[i->j]}.
Positions 1..n of lam; B given as list of position-pairs (u<v) (+ optional singleton, ignored).
"""
from itertools import combinations
from schurcone import *
from dcount import assignments


def first_last(lam):
    n = len(lam)
    first = {}; last = {}
    for p, v in enumerate(lam, 1):
        first.setdefault(v, p); last[v] = p
    return first, last


def count_F(lam, pairs, i, j):
    """pairs: list of (u,v) positions u<v. i<j positions with i=first(lam_i), j=last(lam_j)."""
    n = len(lam)
    inner = list(range(i + 1, j))  # candidates for J besides i
    total = 0
    for k in range(len(inner) + 1):
        for extra in combinations(inner, k):
            J = [i] + list(extra)
            Jset = set(J)
            def succ(u):
                for r in J:
                    if r > u:
                        return r
                return j
            ok = True
            # validity
            for r in J:
                if r == i:
                    continue
                if succ(r) >= r + 2 and not (lam[r - 1] > lam[r]):
                    ok = False; break
            if not ok:
                continue
            # lattice
            for (u, v) in pairs:
                if u < i or v > j:
                    continue
                if (u == i or u not in Jset) and v == succ(u):
                    ok = False; break
            if ok:
                total += 1
    return total


def one_box_moves(lam):
    first, last = first_last(lam)
    res = []
    vals = sorted(set(lam), reverse=True)
    for a in vals:
        for b in vals:
            if b > a:
                continue
            i, j = first[a], last[b]
            if i >= j:
                continue
            mu = list(lam); mu[i - 1] += 1; mu[j - 1] -= 1
            mu = tuple(x for x in mu if x > 0)
            res.append((a, b, i, j, mu))
    return res


if __name__ == "__main__":
    import sys
    Nmax = int(sys.argv[1])
    bad = 0; tested = 0
    for N in range(2, Nmax + 1):
        for lam in partitions(N):
            if len(lam) < 2:
                continue
            for B in nested_with_phi(lam):
                f = schur_product(B)
                asg = assignments(B, lam)
                for (a, b, i, j, mu) in one_box_moves(lam):
                    c = f.get(mu, 0)
                    for bl in asg[:30]:
                        pairs = [tuple(x) for x in bl if len(x) == 2]
                        F = count_F(lam, pairs, i, j)
                        tested += 1
                        if F != c:
                            bad += 1
                            if bad <= 20:
                                print(f"MISMATCH lam={lam} B={B} pairs={pairs} (a,b)=({a},{b}) mu={mu}: c={c} F={F}")
        print(f"N={N}: tested {tested}, mismatches {bad}", flush=True)
