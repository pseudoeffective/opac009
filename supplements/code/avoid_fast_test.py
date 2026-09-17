import sys, time
from schurcone import *
from fastlp import restricted_avoidance_fast
A = canon(eval(sys.argv[1]))
for rho in sorted(set(A), reverse=True):
    t = time.time()
    ok, F = restricted_avoidance_fast(A, rho)
    print(A, "rho =", rho, "restricted avoidance:", ok, f"({time.time()-t:.1f}s)", flush=True)
