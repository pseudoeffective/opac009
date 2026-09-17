import sys, time
from schurcone import *
from avoid3 import restricted_avoidance
A = canon(eval(sys.argv[1])); 
for rho in set(A):
    t = time.time()
    print(A, "rho =", rho, "restricted avoidance:", restricted_avoidance(A, rho), f"({time.time()-t:.0f}s)", flush=True)
