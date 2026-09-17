"""White's d^mu_A(alpha): number of SSYT of shape mu, content lam (sorted), such that for each block
alpha_u of positions, the reading word restricted to letters in alpha_u is a lattice word.
Compare with the true LR coefficient c^mu_A."""
from itertools import product, permutations
from schurcone import *


def ssyt(shape, content):
    """generate SSYT of given shape and content (content = tuple of multiplicities of letters 1..n)
    as list of rows (lists)."""
    n = len(content)
    shape = list(shape)
    rows = [[None] * shape[r] for r in range(len(shape))]
    # fill letter by letter: letter k occupies a horizontal strip
    def rec(k, cur):  # cur = current shape filled with letters < k (list of row lengths)
        if k > n:
            if cur == shape:
                yield [row[:] for row in rows]
            return
        m = content[k - 1]
        # choose horizontal strip of size m added to cur inside shape
        L = len(cur)
        def strip(r, remaining, new):
            if r == len(shape):
                if remaining == 0:
                    yield new
                return
            lo = cur[r] if r < L else 0
            hi = shape[r]
            if r > 0:
                # horizontal strip: new[r] <= cur[r-1] (cells in row r must be under cells filled by < k in row r-1)
                prev = cur[r - 1] if r - 1 < L else 0
                hi = min(hi, prev)
            for a in range(0, min(remaining, hi - lo) + 1):
                yield from strip(r + 1, remaining - a, new + [lo + a])
        for new in strip(0, m, []):
            new_t = [x for x in new]
            # write letter k
            for r in range(len(shape)):
                lo = cur[r] if r < L else 0
                for c in range(lo, new_t[r]):
                    rows[r][c] = k
            # trim
            yield from rec(k + 1, new_t)
    yield from rec(1, [])


def reading_word(T):
    w = []
    for row in T:
        w.extend(reversed(row))
    return w


def is_lattice(word, letters):
    """word restricted to `letters` (sorted list) is lattice: at each prefix, #letters[t] >= #letters[t+1]"""
    cnt = {l: 0 for l in letters}
    for x in word:
        if x in cnt:
            cnt[x] += 1
            idx = letters.index(x)
            if idx > 0 and cnt[letters[idx - 1]] < cnt[x]:
                return False
    return True


def d_count(mu, lam, blocks):
    """blocks: list of sorted lists of positions (1-based letters)"""
    total = 0
    for T in ssyt(mu, lam):
        w = reading_word(T)
        if all(is_lattice(w, b) for b in blocks):
            total += 1
    return total


def assignments(B, lam):
    """all assignments of positions of lam to the pairs of B (each pair (a,b) gets a position of value a
    and a position of value b), as list of blocks; deduplicated."""
    positions = {}
    for p, v in enumerate(lam, 1):
        positions.setdefault(v, []).append(p)
    # for each value v, positions[v]; the pairs of B use value v some number of times; we need to distribute
    # positions of value v among the "slots" (pair, which part) with that value.
    slots = []  # (pair index, part index, value)
    for idx, rho in enumerate(B):
        for t, v in enumerate(rho):
            slots.append((idx, t, v))
    res = set()
    per_value = {}
    for s in slots:
        per_value.setdefault(s[2], []).append(s)
    values = list(per_value)
    def rec(vi, assign):
        if vi == len(values):
            blocks = [[] for _ in B]
            for (idx, t, v), p in assign.items():
                blocks[idx].append(p)
            res.add(tuple(tuple(sorted(b)) for b in blocks))
            return
        v = values[vi]
        sl = per_value[v]
        for perm in permutations(positions[v]):
            a2 = dict(assign)
            for s, p in zip(sl, perm):
                a2[s] = p
            rec(vi + 1, a2)
    rec(0, {})
    return sorted(res)


if __name__ == "__main__":
    import sys
    lam = eval(sys.argv[1])
    N = sum(lam)
    Bs = nested_with_phi(lam)
    mus = [m for m in interval(lam, lam_plus(lam)) if m != lam]
    for B in Bs:
        f = schur_product(B)
        for mu in mus:
            c = f.get(mu, 0)
            ds = [(bl, d_count(mu, lam, [list(b) for b in bl])) for bl in assignments(B, lam)]
            flag = "" if all(d == c for _, d in ds) else "   <-- MISMATCH"
            print(f"B={B} mu={mu}: c={c}, d={[d for _, d in ds]}{flag}")
            if flag:
                for bl, d in ds:
                    print("      ", bl, d)
