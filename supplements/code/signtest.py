import sys
from schurcone import *
from signcrit import *
A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam); M = N + 1
comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
print(f"A = {A}, lam = {lam}, {len(comps)} competitors")
found = []
for p in range(1, M + 1):
    for q in range(p, M + 1):
        rank = block_reversal_rank(M, p, q)
        ok, pr = check_levels(A, rank, comps)
        if ok: found.append((p, q))
print("block reversals [p,q] that work (level-aware):", found)
if len(sys.argv) > 2:
    p, q = eval(sys.argv[2]); check_levels(A, block_reversal_rank(M, p, q), comps, verbose=True)
