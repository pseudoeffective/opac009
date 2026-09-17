import time
from schurcone import *
for N in range(7, 13):
    t = time.time()
    gens = nested_multisets(N)
    bad = [A for A in gens if not is_extreme(A, gens, N)[0]]
    print(N, len(gens), "nested; non-extreme:", bad, f"{time.time()-t:.1f}s", flush=True)
