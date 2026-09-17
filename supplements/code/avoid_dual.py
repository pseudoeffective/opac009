"""Find a nonneg combination of nested s_B (B not containing rho, phi(B) ⊵ lam) lying in s_rho*Lambda."""
import sys
from fractions import Fraction
import cdd.gmp as cdd
from schurcone import *
A = canon(eval(sys.argv[1])); rho = tuple(eval(sys.argv[2]))
lam = phi(A); N = sum(lam)
parts = list(partitions(N)); idx = {p:i for i,p in enumerate(parts)}; d = len(parts)
gens = [B for B in nested_multisets(N) if rho not in B and dominates(phi(B), lam)]
nus = list(partitions(N - sum(rho)))
# variables: c_B >= 0 (len gens), y_nu free (len nus); constraint: sum c_B s_B - sum y_nu s_rho s_nu = 0; sum c_B = 1
m = len(gens); k = len(nus)
cols = [vec(schur_product(B), idx) for B in gens] + [[-x for x in vec(mul_schur2({nu:1}, rho), idx)] for nu in nus]
rows = []
for i in range(d):
    rows.append([Fraction(0)] + [Fraction(cols[j][i]) for j in range(m+k)])
rows.append([Fraction(-1)] + [Fraction(1)]*m + [Fraction(0)]*k)
neq = len(rows)
for j in range(m):
    r = [Fraction(0)]*(m+k+1); r[j+1] = Fraction(1); rows.append(r)
mat = cdd.matrix_from_array(rows, lin_set=range(neq), rep_type=cdd.RepType.INEQUALITY,
                            obj_type=cdd.LPObjType.MIN, obj_func=[Fraction(0)]*(m+k+1))
lp = cdd.linprog_from_matrix(mat); cdd.linprog_solve(lp)
print("status", lp.status)
if lp.status == cdd.LPStatusType.OPTIMAL:
    sol = lp.primal_solution
    print("combination:")
    for j in range(m):
        if sol[j] != 0: print("   ", sol[j], gens[j])
    print(" = s_rho * (")
    for j in range(k):
        if sol[m+j] != 0: print("   ", sol[m+j], "s", nus[j])
    print(" )")
