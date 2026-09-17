"""Restricted level competitors L*(A): B in L(A) with n_{I'}(B) >= n_{I'}(A) for all intervals I'.
Test simple criteria against L*(A)."""
import sys, time, itertools
from schurcone import *
from epsjet import inside_pairs
from levelid import poly_of
from coarse import ordered_set_partitions, sgn

def Lstar(A, comps):
    lam = phi(A); b0, a0 = lam[-1], lam[0]
    AI = inside_pairs(A, b0, a0); n = len(AI)
    vals = sorted(set(lam))
    L = [B for B in comps if len(inside_pairs(B, b0, a0)) == n]
    out = []
    for B in L:
        ok = True
        for x in vals:
            for y in vals:
                if y < x: continue
                if len(inside_pairs(B, x, y)) < len(inside_pairs(A, x, y)): ok = False; break
            if not ok: break
        if ok: out.append(B)
    return out

def bottom_ok(A, Ls):
    lam = phi(A); b0, a0 = lam[-1], lam[0]; T = list(range(b0, a0+2))
    AI = inside_pairs(A, b0, a0); mono = tuple(sorted(b for (a,b) in AI))
    return all(poly_of(inside_pairs(B, b0, a0), T).get(mono, 0) <= 0 for B in Ls)

def top_ok(A, Ls):
    lam = phi(A); b0, a0 = lam[-1], lam[0]; T = list(range(b0, a0+2))
    AI = inside_pairs(A, b0, a0); mono = tuple(sorted(a+1 for (a,b) in AI)); s = (-1)**len(AI)
    return all(s * poly_of(inside_pairs(B, b0, a0), T).get(mono, 0) <= 0 for B in Ls)

if __name__ == "__main__":
    Nmax = int(sys.argv[1])
    for N in range(4, Nmax + 1):
        gens = nested_multisets(N); t = time.time(); fb = []; ft = []; fboth = []; cnt = 0
        for A in gens:
            if len(A) < 2: continue
            cnt += 1
            lam = phi(A); comps = [B for B in gens if B != A and dominates(phi(B), lam)]
            Ls = Lstar(A, comps)
            b = bottom_ok(A, Ls); tp = top_ok(A, Ls)
            if not b: fb.append(A)
            if not tp: ft.append(A)
            if not (b or tp): fboth.append(A)
        print(f"N={N}: {cnt} A; bottom fails {len(fb)}, top fails {len(ft)}, both fail {len(fboth)} [{time.time()-t:.0f}s]", flush=True)
        for A in fboth[:8]: print("      both fail:", A)
