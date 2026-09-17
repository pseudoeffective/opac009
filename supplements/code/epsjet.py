"""epsilon-jet of the interval specialisation.

Base: psi_I kills pairs inside I=[b0,a0] (d constant on T=[b0,a0+1]); perturb d_k = c + eps*e_k on T.
Then psi_eps(s_B) = eps^{n_I(B)} L_B + O(eps^{n_I(B)+1}), sign(L_B) = prod_{(a,b) in B_I} sign(e_b - e_{a+1}).
Criterion (given a linear order on T encoded by e):
  (a) sign(L_A) = +;
  (b) every competitor B with n_I(B) = n_I(A) has sign(L_B) = -;
  (c) for every m < n_I(A), all competitors with n_I(B) = m have the same sign.
"""
import sys, itertools, random
from collections import defaultdict
from schurcone import *


def inside_pairs(B, b0, a0):
    return [r for r in B if len(r) == 2 and b0 <= r[1] <= r[0] <= a0]


def sign_of(pairs, e):
    s = 1
    for (a, b) in pairs:
        s *= 1 if e[b] > e[a + 1] else -1
    return s


def criterion(A, comps, b0, a0, e):
    AI = inside_pairs(A, b0, a0); nA = len(AI)
    if sign_of(AI, e) != 1:
        return False, "a"
    levels = defaultdict(set)
    for B in comps:
        BI = inside_pairs(B, b0, a0)
        levels[len(BI)].add(sign_of(BI, e))
    if 1 in levels.get(nA, set()):
        return False, "b"
    for m in range(nA):
        if len(levels.get(m, set())) > 1:
            return False, f"c{m}"
    return True, "ok"


def search(A, b0, a0, comps, tries=20000, exhaustive=False):
    T = list(range(b0, a0 + 2))
    if exhaustive:
        it = itertools.permutations(T)
    else:
        def gen():
            for _ in range(tries):
                p = T[:]; random.shuffle(p); yield tuple(p)
        it = gen()
    reasons = defaultdict(int)
    for perm in it:
        e = {k: i for i, k in enumerate(perm)}   # e-value = position in perm
        ok, why = criterion(A, comps, b0, a0, e)
        reasons[why] += 1
        if ok:
            return perm, reasons
    return None, reasons


if __name__ == "__main__":
    A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam)
    comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
    b0 = int(sys.argv[2]) if len(sys.argv) > 2 else lam[-1]
    a0 = int(sys.argv[3]) if len(sys.argv) > 3 else lam[0]
    print(f"A={A} lam={lam} I=[{b0},{a0}] T size {a0-b0+2}; {len(comps)} competitors; n_I(A)={len(inside_pairs(A,b0,a0))}")
    exh = (a0 - b0 + 2) <= 8
    perm, reasons = search(A, b0, a0, comps, exhaustive=exh)
    print("exhaustive" if exh else "random", "search; found order (increasing e):", perm)
    print("failure reasons:", dict(reasons))
