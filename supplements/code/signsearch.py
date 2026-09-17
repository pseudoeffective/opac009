import sys, random, itertools
from schurcone import *
from signcrit import *
A = canon(eval(sys.argv[1])); lam = phi(A); N = sum(lam); M = N + 1
comps = [B for B in nested_multisets(N) if B != A and dominates(phi(B), lam)]
print(f"A = {A}, lam = {lam}, {len(comps)} competitors, M = {M}")
mode = sys.argv[2] if len(sys.argv) > 2 else "random"
found = 0
if mode == "exhaustive":
    it = itertools.permutations(range(1, M + 1))
else:
    def it_():
        while True:
            o = list(range(1, M + 1)); random.shuffle(o); yield o
    it = it_()
tried = 0
for order in it:
    tried += 1
    rank = rank_from_order(order)
    ok, pr = check_levels(A, rank, comps)
    if ok:
        found += 1
        print("FOUND order:", order); 
        check_levels(A, rank, comps, verbose=True)
        if found >= 3: break
    if mode != "exhaustive" and tried >= int(sys.argv[3]) if len(sys.argv) > 3 else tried >= 200000: break
print("tried", tried, "found", found)
