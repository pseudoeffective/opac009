# OPAC-009 — working notes (2026-09-09, updated after session 2)

**Headline (session 2): White's conjecture is proved in full — s_A is extreme iff A is nested.**
Proof = Level Identity Theorem (specialisation ε-jet) + monomial-coefficient lemmas + dominance; see memo §§2–5
(refereed twice by independent subagents; all issues raised were editorial and have been fixed).
Everything in sections 1–4 below is now superseded by / a corollary of the main theorem, but retained.

## 0. The proof of the main theorem (summary; full details in memo/opac009-memo.tex §§2–5)

Setup: A nested, λ=φ(A), n=ℓ(λ), m = #pairs of A (n ∈ {2m, 2m+1}), w=λ_n, u=λ_1, I=[w,u], T=[w,u+1].
Suppose s_A = Σ_{B∈S} c_B s_B with B nested ≠ A, φ(B) ⊵ λ (Lemma red(a)).

**Level Identity Theorem.** For ANY decomposition s_A = Σ c_B s_B (c_B>0, B arbitrary multisets) and any
interval I=[x,y], T=[x,y+1]: base d⁰ constant on T, strictly decreasing at every step leaving T; perturb
d = d⁰ + t·e (e ∈ ℝ^T); then ψ_t(s_B) = t^{n_I(B)} P_{B_I}(e) R_B(t), P_M(δ)=∏_{(a,b)∈M}(δ_b − δ_{a+1}),
R_B entire, R_B(0)=K_B>0 independent of e.  Comparing coefficients of t^k:
  (i) n_I(B) ≥ n_I(A) for all B in the support (lowest-order term at strictly decreasing e is >0);
  (ii) K_A P_{A_I} = Σ_{n_I(B)=n_I(A)} c_B K_B P_{B_I} as polynomials in δ_k, k∈T.
Verified on real decompositions of all non-nested A, N ≤ 9, all intervals (verify_levelid.py).

**Step 1 (shape).** (i) for I: 2m ≤ 2n_I(B) ≤ ℓ(φ(B)) ≤ n ≤ 2m+1 ⟹ B = M_B ⊎ {(c'_B)} (c'_B = N − ‖M_B‖ ≥ 0,
singleton only if c'_B ≥ 1), M_B = m pairs inside I.  (i) for all [x,y] ⊆ I: n_{[x,y]}(M_B) ≥ n_{[x,y]}(A);
hence β(M_B)_{≤x} ≤ β(A)_{≤x}, τ(M_B)_{≤y} ≥ τ(A)_{≤y}.
**Step 2.** (ii): P_A = Σ c''_B P_{M_B}.

**Lemmas.** (Sum) X_{≤x} ≤ Y_{≤x} ∀x, |X|=|Y| ⟹ ΣX ≥ ΣY, equality iff X=Y.
(Bottom) β(M)_{≤x} ≤ β(A')_{≤x} ⟹ coeff of ∏_{A'} δ_b in P_M = [β(M)=β(A')].
(Top) τ(M)_{≤y} ≥ τ(A')_{≤y} ⟹ coeff of ∏_{A'} δ_{a+1} in P_M = (−1)^{|M|}[τ(M)=τ(A')].
(Unique) a nested multiset of pairs is determined by (τ, β): max bottom b* with mult r pairs with the r smallest
tops ≥ b* (interlacing (i) forces every (a,b) with b<b*≤a to have a ≥ a*); induct.

**Step 3 (n=2m).** Bottom monomial: survivors have β(M_B)=β(A); weight ⟹ Στ equal ⟹ τ(M_B)=τ(A) ⟹ M_B=A. ⊥.

**Step 4 (n=2m+1, singleton c).** A = A↓ ⊎ A= ⊎ A↑ (tops < c / (c,c) / bottoms > c) by White's condition (ii).
(i) on [w,c−1],[c,c],[c+1,u] ⟹ M_B = M↓ ⊎ A= ⊎ M↑ (same sizes), with τ(M↓)_{≤y} ≥ τ(A↓)_{≤y},
β(M↓)_{≤x} ≤ β(A↓)_{≤x}, β(M↑)_{≤x} ≤ β(A↑)_{≤x}, τ(M↑)_{≤y} ≥ τ(A↑)_{≤y}.  Cancel P_{A=}.
Monomial θ = ∏_{A↓} δ_{a+1} · ∏_{A↑} δ_b (variables ≤ c and ≥ c+1 separate) has coefficient (−1)^{m↓} in
P_{A↓}P_{A↑} and (−1)^{m↓}[τ(M↓)=τ(A↓), β(M↑)=β(A↑)] in P_{M↓}P_{M↑}.  Survivors: p = Σβ(M↓)−Σβ(A↓) ≥ 0,
q = Στ(A↑)−Στ(M↑) ≥ 0, c'_B = c − p + q.  Dominance on the top block (2m↑ largest parts): c'_B ≤ c and q=0
⟹ M↑=A↑; on the bottom block (2m↓ smallest parts of φ(B) padded): c'_B = c, p = 0 ⟹ M↓ = A↓.  So B = A. ⊥.
Dominance is essential (A={(5,4),(3),(2,1)}, B={(4,4),(4),(2,1)} survives everything else).
Checked: oddmono.py — for all nested A, N ≤ 18, every competitor satisfying the interval inequalities has
coefficient 0 at the relevant monomial.

**Not obtained:** an explicit separating linear functional (White's Conjecture 11 on [λ,λ⁺⁺] stays open);
possible follow-up: supports of the Taylor-coefficient functionals L_{k,e}.

Code: `code/schurcone.py` (exact Schur expansions via Pieri/Jacobi–Trudi, nested enumeration, exact LP via
pycddlib), `vertices.py` (certificate polytopes), `dcount.py`/`formulaF.py` (LR path formula tests),
`avoid*.py` (avoidance LPs), `hybrid*.py`, `coverage.py`.

Sanity: reproduces White's 13 extreme rays for N=6; all nested A extreme for N ≤ 12 (exact LP).

## 1. Specialisation argument — PROVEN

Ring homomorphisms ψ_g: Λ → ℝ, h_k ↦ g(k) (h's algebraically independent), g(0)=1.
ψ_g(s_{(a,b)}) = g(a)g(b) − g(a+1)g(b−1) =: G(a,b).  With g = exp(u), u(k)=Σ_{t≤k} δ(t):
G(a,b) > 0 ⟺ δ(b) > δ(a+1);  G(a,b)=0 ⟺ δ(b)=δ(a+1).

**Lemma S (interval lemma).** For an interval I=[b,a] choose δ strictly decreasing except constant on
{b,…,a+1}. Then ψ(s_{(x,y)}) = 0 if b ≤ y ≤ x ≤ a, and > 0 otherwise; ψ(h_c) > 0.
Consequence: if s_A = Σ c_B s_B (c_B ≥ 0) and A has a pair inside I, every B with c_B > 0 has a pair inside I.

**Theorem 1.** s_{A''} extreme in C^2_{N''} ⟹ s_{A''∪{(p,p)}} extreme in C^2_{N''+2p}.
Proof: I=[p,p]; only pairs inside are (p,p); so all B with c_B>0 contain (p,p); cancel s_{(p,p)} (Λ a domain);
contradiction with extremality of A''.

**Corollary.** Family G = {nested A : after deleting all pairs (v,v), the remaining parts are distinct} is extreme
(White + Thm 1). Nestedness: (v,v) never interlaces with anything. G covers 98/127 nested A at N=12, 180/251 at N=14.
The complement ("hard core") = nested A with some value that is the top of ≥2 pairs or the bottom of ≥2 pairs.

Restricted avoidance (∃F ∈ ker s_ρ^⊥ with F>0 on nested B∌ρ, φ(B)⊵λ) holds numerically for many (A,ρ) with ρ=(a,b),
a>b, but fails whenever ρ∪B'' non-nested with φ ⊵ λ exists (e.g. singleton (|A''|) when |A''| ≤ a). No positive
log-specialisation can kill exactly ρ≠(p,p) (chain obstruction δ(b)>δ(x)>δ(a+1)).

Theorem 2 (interval rigidity: A unique in SSP_λ with pairs inside all its intervals ⟹ extreme) is true but adds
nothing beyond G numerically (N ≤ 14).

## 2. Formula (F) for one-box-move LR coefficients — PROVEN (proof in memo)

λ = (λ_1 ≥ … ≥ λ_n), positions 1..n; first(a), last(b). μ = λ[i→j] := λ + e_i − e_j with i=first(λ_i), j=last(λ_j), i<j.
B ∈ SP^2_λ with an arbitrary assignment of positions (pairs (u,v), u<v; singleton ignored).
A *path* is J ⊆ [i,j−1] with i ∈ J; succ_J(u) = min{r ∈ J∪{j}: r>u}. J is *valid* if every r ∈ J\{i} with
succ_J(r) ≥ r+2 satisfies λ_r > λ_{r+1}. A pair (u,v) of B with i ≤ u < v ≤ j is *violated* by J iff
(u ∉ J or u = i) and succ_J(u) = v.
   c_B^{λ[i→j]} = #{valid J with no violated pair}.
Verified exactly for all nested B, all λ ⊢ N ≤ 10, all assignments. Multi-segment version (μ ∈ [λ,λ⁺]):
product over segments (verified N ≤ 10).
Proof: s_B = Σ_S (−1)^{|S|} h_{ν_S} (Jacobi–Trudi); K_{μ,ν_S} ≠ 0 only for S with distinct top values, distinct
bottom values, pairwise disjoint λ-intervals [first(a),last(b)) ⊆ [i,j); K_{μ,ν} for P(μ)−P(ν)∈{0,1}
counted by "intruder paths" (one intruder cell per level in each segment; stay allowed only at strict corners);
sign-reversing bijection (S,J') ↔ (J, S ⊆ V(J)): J = J' ∪ {i} ∪ ⋃_{(u,v)∈S} [v, last(λ_v)) ∩ [i,j−1].
Special case: White's Lemma 15 / Theorem 17 (distinct λ).

Note: White's d-count (sorted content, alphabets) equals c for all μ ∈ [λ,λ⁺] (numerically, N ≤ 11) but NOT for
two-box shapes like λ⁺⁺ (this is exactly the Gaetz et al. error).

## 3. Applications of (F)

Four-term certificates on [λ,λ⁺]:
* A={(j,i),(j,i)}: f = s_λ + s_{λ[j→j]} + s_{λ[i→i]} − s_{λ[j→i]}  (f(A)=1, f(B_0)=0, f(C)=0 on SSP_{λ[j→j]},
  SSP_{λ[i→i]}, f ≤ 0 elsewhere).  ⟹ extreme.  (Gaetz Thm 5.1, now proven.)
* A={(v,x),(v,y)} (x>y): f = s_λ + s_{λ[v→v]} + s_{λ[x→y]} − s_{λ[v→y]}.
* A={(x,v),(y,v)} (x>y): f = s_λ + s_{λ[v→v]} + s_{λ[x→y]} − s_{λ[x→v]}.
(Second and third verified numerically; proofs analogous — to be written if included.)

## 4. Earlier 'what remains' (now superseded by the main theorem)

Hybrid criterion (Lemma S restriction to R + one-box or [λ,λ⁺] certificate) covers all nested A with N ≤ 12 except
those with a value shared by ≥3 pairs (as tops or bottoms). But for N=20, {(5,2),(5,1),(4,3)} needs two-segment
shapes ([5→4][2→1], [5→5][3→1]); certificate:
  −s[5>1] + s[5>4][2>1] + s[5>4] + s[5>5][3>1] + s[5>5] + s[4>1] + s[3>1]  (hybrid, [λ,λ⁺]).
Uniform certificates for "two tops (v,x),(v,y) + distinct rest" not yet found.
White's Conjecture 11 ([λ,λ⁺⁺] suffices) remains the natural frame for the hard core.
