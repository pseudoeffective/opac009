import sys, time
from schurcone import *
A = eval(sys.argv[1])
print("nested:", is_nested(A))
lam = phi(A); N = sum(lam)
t=time.time()
gens = nested_multisets(N)
above = [B for B in gens if dominates(phi(B), lam)]
print(len(gens), "nested gens;", len(above), "with phi(B) ⊵ lam", f"{time.time()-t:.1f}s")
sup = [m for m in partitions(N) if dominates(m, lam)]
t=time.time()
sep = separating_functional(A, support=sup, generators=above, N=N)
print("separator from above:", sep, f"{time.time()-t:.1f}s")
