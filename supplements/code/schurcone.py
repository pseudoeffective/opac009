"""
Exact toolkit for the (N,2)-Schur cone / cone of log-concavity (OPAC-009).

Symmetric functions of degree N are dicts {partition(tuple): int} in the Schur basis.
"""
from functools import lru_cache
from fractions import Fraction
from itertools import combinations
import cdd.gmp as cdd


# ---------- partitions ----------
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def dominates(lam, mu):
    """lam ⊵ mu"""
    a = b = 0
    for i in range(max(len(lam), len(mu))):
        a += lam[i] if i < len(lam) else 0
        b += mu[i] if i < len(mu) else 0
        if a < b:
            return False
    return True


# ---------- Pieri ----------
@lru_cache(maxsize=None)
def pieri(lam, k):
    """s_lam * h_k as tuple of partitions (each with multiplicity 1) -- horizontal strips."""
    if k == 0:
        return (lam,)
    res = []
    n = len(lam)
    lam_ext = lam + (0,)
    # choose a_i = number of boxes added to row i (i=0..n), with a_i <= lam_{i-1}-lam_i (i>=1), sum = k
    def rec(i, remaining, cur):
        if i == n + 1:
            if remaining == 0:
                mu = tuple(x for x in cur if x > 0)
                res.append(mu)
            return
        cap = remaining if i == 0 else min(remaining, lam[i - 1] - lam_ext[i])
        for a in range(cap, -1, -1):
            rec(i + 1, remaining - a, cur + [lam_ext[i] + a])
    rec(0, k, [])
    return tuple(res)


def mul_h(f, k):
    out = {}
    for lam, c in f.items():
        for mu in pieri(lam, k):
            out[mu] = out.get(mu, 0) + c
    return {m: c for m, c in out.items() if c != 0}


def mul_schur2(f, rho):
    """multiply f by s_rho, rho a partition with <= 2 parts, via Jacobi-Trudi."""
    if len(rho) == 1:
        return mul_h(f, rho[0])
    a, b = rho
    t1 = mul_h(mul_h(f, a), b)
    t2 = mul_h(mul_h(f, a + 1), b - 1)
    out = dict(t1)
    for m, c in t2.items():
        out[m] = out.get(m, 0) - c
    return {m: c for m, c in out.items() if c != 0}


def canon(A):
    return tuple(sorted((tuple(r) for r in A), reverse=True))


def schur_product(A):
    return _schur_product(canon(A))


@lru_cache(maxsize=None)
def _schur_product(A):
    """A: sorted tuple of partitions (tuples). Returns dict of Schur expansion of s_A."""
    f = {(): 1}
    for rho in A:
        f = mul_schur2(f, rho)
    return f


# ---------- nested multisets ----------
def interlaced(lam, mu):
    """White's conditions (i)-(iii) on an unordered pair."""
    def cond(l, m):
        if len(l) == 2 and len(m) == 2:
            return l[0] > m[0] >= l[1] > m[1]
        if len(l) == 2 and len(m) == 1:
            return l[0] > l[1] and l[0] >= m[0] >= l[1]
        if len(l) == 1 and len(m) == 1:
            return True
        return False
    return cond(lam, mu) or cond(mu, lam)


def is_nested(A):
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            if interlaced(A[i], A[j]):
                return False
    return True


def phi(A):
    return tuple(sorted((p for rho in A for p in rho), reverse=True))


def all_multisets(N):
    """all multisets of partitions with <=2 parts and total weight N, as sorted tuples (descending)."""
    parts = [p for n in range(1, N + 1) for p in partitions(n, None) if len(p) <= 2]
    parts.sort(reverse=True)
    res = []
    def rec(idx, remaining, cur):
        if remaining == 0:
            res.append(tuple(cur))
            return
        for j in range(idx, len(parts)):
            p = parts[j]
            if sum(p) <= remaining:
                rec(j, remaining - sum(p), cur + [p])
    rec(0, N, [])
    return res


def nested_multisets(N):
    return [A for A in all_multisets(N) if is_nested(A)]


def nested_with_phi(lam):
    """all nested A with phi(A) = lam"""
    N = sum(lam)
    return [A for A in nested_multisets(N) if phi(A) == lam]


# ---------- vectors ----------
def vec(f, index):
    v = [0] * len(index)
    for m, c in f.items():
        v[index[m]] = c
    return v


# ---------- LP / extremality (exact, via cdd) ----------
def is_extreme(A, generators=None, N=None):
    """Exact test: is s_A NOT a nonneg combination of the other nested s_B?
    Returns (extreme:bool, witness) where witness is a nonneg combination (dict B->coef) if not extreme,
    or None."""
    A = canon(A)
    if N is None:
        N = sum(phi(A))
    if generators is None:
        generators = nested_multisets(N)
    others = [B for B in generators if canon(B) != A]
    parts = list(partitions(N))
    index = {p: i for i, p in enumerate(parts)}
    target = vec(schur_product(A), index)
    cols = [vec(schur_product(B), index) for B in others]
    # LP: find x >= 0 with sum_B x_B cols_B = target.  cdd H-representation: rows are b + A x >= 0 (ineq) / = 0 (eq)
    rows = []
    d = len(parts)
    m = len(others)
    # equalities: -target_i + sum_B cols_B[i] x_B = 0
    for i in range(d):
        rows.append([Fraction(-target[i])] + [Fraction(cols[j][i]) for j in range(m)])
    neq = d
    for j in range(m):
        r = [Fraction(0)] * (m + 1)
        r[j + 1] = Fraction(1)
        rows.append(r)
    mat = cdd.matrix_from_array(rows, lin_set=range(neq), rep_type=cdd.RepType.INEQUALITY,
                                obj_type=cdd.LPObjType.MIN, obj_func=[Fraction(0)] * (m + 1))
    lp = cdd.linprog_from_matrix(mat)
    cdd.linprog_solve(lp)
    if lp.status == cdd.LPStatusType.OPTIMAL:
        sol = lp.primal_solution
        wit = {others[j]: sol[j] for j in range(m) if sol[j] != 0}
        return False, wit
    elif lp.status == cdd.LPStatusType.INCONSISTENT:
        return True, None
    else:
        raise RuntimeError(f"LP status {lp.status}")


def separating_functional(A, support=None, generators=None, N=None, objective=None):
    """Find phi (dict partition->Fraction), supported on `support` (list of partitions; default all),
    with phi(s_A) = 1 and phi(s_B) <= 0 for all other nested B.  Returns dict or None.
    objective: optional dict partition->Fraction to minimise (default: minimise sum of |phi| via
    split variables is not implemented; default objective 0)."""
    A = canon(A)
    if N is None:
        N = sum(phi(A))
    if generators is None:
        generators = nested_multisets(N)
    parts = list(partitions(N)) if support is None else list(support)
    index = {p: i for i, p in enumerate(parts)}
    d = len(parts)
    def restricted_vec(f):
        return [Fraction(f.get(p, 0)) for p in parts]
    tA = restricted_vec(schur_product(A))
    rows = []
    # equality phi . tA = 1  ->  -1 + tA.phi = 0
    rows.append([Fraction(-1)] + tA)
    for B in generators:
        if canon(B) == A:
            continue
        tB = restricted_vec(schur_product(B))
        # -tB.phi >= 0
        rows.append([Fraction(0)] + [-x for x in tB])
    if objective is None:
        obj = [Fraction(0)] * (d + 1)
    else:
        obj = [Fraction(0)] + [Fraction(objective.get(p, 0)) for p in parts]
    mat = cdd.matrix_from_array(rows, lin_set=[0], rep_type=cdd.RepType.INEQUALITY,
                                obj_type=cdd.LPObjType.MIN, obj_func=obj)
    lp = cdd.linprog_from_matrix(mat)
    cdd.linprog_solve(lp)
    if lp.status == cdd.LPStatusType.OPTIMAL:
        sol = lp.primal_solution
        return {parts[i]: sol[i] for i in range(d) if sol[i] != 0}
    return None


def interval(lam, rho):
    """partitions mu with lam ⊴ mu ⊴ rho"""
    N = sum(lam)
    return [mu for mu in partitions(N) if dominates(mu, lam) and dominates(rho, mu)]


def lam_plus(lam):
    l = list(lam)
    l[0] += 1
    l[-1] -= 1
    return tuple(x for x in l if x > 0)


def lam_plusplus(lam):
    l = list(lam)
    m = len(l)
    k = m // 2
    for i in range(k):
        l[i] += 1
    for i in range(m - k, m):
        l[i] -= 1
    return tuple(x for x in l if x > 0)


if __name__ == "__main__":
    # sanity: White's list of the 13 extreme rays for N=6
    N = 6
    gens = nested_multisets(N)
    ext = [A for A in gens if is_extreme(A, gens, N)[0]]
    print(len(gens), "nested;", len(ext), "extreme")
    for A in ext:
        print(A)
    # check s_(3,1) s_(2) = s_(3,2)s_(1) + s_(1,1)s_(4)
    lhs = schur_product(((3, 1), (2,)))
    r1 = schur_product(((3, 2), (1,)))
    r2 = schur_product(((4,), (1, 1)))
    tot = dict(r1)
    for k, v in r2.items():
        tot[k] = tot.get(k, 0) + v
    print("syzygy ok:", lhs == tot)
