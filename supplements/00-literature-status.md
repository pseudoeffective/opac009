# OPAC-009 — literature status (2026-09-09)

## White, arXiv:0903.2831v2 (Nov 2014) — read in full (arXiv HTML)

Solid, self-contained. Results:

* Thm 4: not nested ⇒ not extreme (four Jacobi–Trudi syzygies (1)–(4)).
* Farkas reformulation (Thm 7): A extreme iff ∃ f with ⟨f,s_A⟩>0, ⟨f,s_B⟩≤0 for all *nested* B≠A.
  (It suffices to separate from nested B since the non-nested s_B are positive combos of nested ones.)
* Lemma 8 / Lemma 10: it suffices to separate "from above" (only B with φ(B) ⊵ λ) and even only
  on a dominance interval [λ, ρ]; the rest is fixed by a greedy dual-order-ideal correction.
* Thm 12: λ=φ(A) distinct, A nested ⇒ separator exists on [λ, λ⁺], λ⁺=(λ1+1,λ2,…,λm−1).
  Proof: chain X_0 ⊇ … ⊇ X_m={A} indexed by an inside-out ordering of the pairs of A;
  f_i = (c_A^{λ[ρ]}+1) s_λ − s_{λ[ρ]} partially separates, via Lemma 15 (c_A^{λ[ρ]}+1 = c_B^{λ[ρ]}
  when ρ∈A, ρ∉B, A,B agree within ρ).
* Lemma 15 comes from Thm 17: for λ distinct and shape λ[ν], c^{λ[ν]}_A = d^{λ[ν]}_A where d counts
  SSYT with *sorted* content λ and each alphabet-subword lattice — proved by tableau switching
  (Benkart–Sottile–Stroomer). White stresses d ≠ c in general (example d=15, c=13).
* Remark: for λ=2³1³, A={(2,1)³}, NO separator exists on [λ,λ⁺]. Conjecture 11: a separator always
  exists on [λ, λ⁺⁺] (λ⁺⁺ = add 1 to first ⌊m/2⌋ parts, subtract 1 from last ⌊m/2⌋ parts). Verified N≤20.
* §7: nested A with |φ(A)|=2m parts ↔ plane partitions of shape (m,m) with the given parts (bijection).

## OPAC blog post (White, Sept 2019) — realopacblog.wordpress.com/2019/09/22/...

Same statements as the paper (Thm: not nested ⇒ not extreme; Thm: nested + distinct ⇒ extreme;
OPAC-009 = nested ⇒ extreme; OPAC-010 = k ≥ 3). No mention of Gaetz et al., no comments, no
further families. Confirms the problem is open exactly as in the handoff.

## Gaetz–Meyer–Tam–Wimberley–Yao–Zhu, arXiv:1409.4859 — **WITHDRAWN (v2, Sept 2016)**

Comment: "withdrawn by the authors due to a misinterpretation of the generalized Littlewood–Richardson
rule in several proofs". Their Thm 2.9 states the multi-LR rule with content = *sorted* φ(A) and lattice
condition on the alphabet subwords — that is White's d^λ_A, not c^λ_A. So everything resting on it is
unproven:

* Thm 4.1 (A extreme ⇒ A∪{(p,p)} extreme), Cor 4.2 (completely separated nested A extreme),
* Thm 5.1 ({(j,i),(j,i)} extreme),
* Conjectures 6.1 (juxtaposition: λ_n > μ_1, at most one of n,m odd), 6.2 (flanking by ρ=(ρ1>ρ2) with
  ρ1 ≥ λ1, λn ≥ ρ2) and Lemma 6.3 (6.1+6.2 (+4.1?) ⇒ White's conjecture): the conjectures are still
  meaningful; the reduction should be re-derived if used.

Consequence (confirmed with David 2026-09-09): the only proven cases of OPAC-009 are White's
distinct-parts theorem. Any rigorously proven infinite family with repeated parts is new, including
the families Gaetz et al. claimed. White v2 is taken as given.
