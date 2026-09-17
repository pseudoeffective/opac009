"""Single-pair intervals: for (a,b) in A with n_{[b,a]}(A)=1, level-1 competitors give arcs b'->a'+1;
A extreme if no directed path from b to a+1 using competitor arcs (excluding the arc b->a+1 itself only if
no competitor has exactly the pair (a,b) inside... if some competitor has (a,b) as its unique inside pair, fail)."""
import sys, time
from collections import defaultdict
from schurcone import *
from epsjet import inside_pairs

def chain_test(A, comps):
    results = []
    for (a, b) in set(r for r in A if len(r) == 2):
        if len(inside_pairs(A, b, a)) != 1: continue
        arcs = set()
        for B in comps:
            BI = inside_pairs(B, b, a)
            if len(BI) == 1:
                (a2, b2) = BI[0]; arcs.add((b2, a2 + 1))
        # path from b to a+1?
        adj = defaultdict(list)
        for (s, t) in arcs: adj[s].append(t)
        seen = {b}; stack = [b]; reach = False
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y == a + 1: reach = True
                if y not in seen: seen.add(y); stack.append(y)
        results.append(((a, b), not reach, sorted(arcs)))
    return results

if __name__ == "__main__":
    if sys.argv[1].startswith("N="):
        Nmax = int(sys.argv[1][2:])
        for N in range(4, Nmax + 1):
            gens = nested_multisets(N); t = time.time(); fails = []
            for A in gens:
                if len(A) < 2: continue
                lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
                res = chain_test(A, comps)
                if not any(ok for _, ok, _ in res): fails.append(A)
            print(f"N={N}: {len(gens)} nested, {len(fails)} not proved by single-pair chain test [{time.time()-t:.0f}s]", flush=True)
            for A in fails[:12]: print("      ", A)
    else:
        A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam)
        comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
        for pair, ok, arcs in chain_test(A, comps):
            print(A, "pair", pair, "no chain:" , ok, "arcs:", arcs)
