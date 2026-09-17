# OPAC-009 session transcript (2026-09-09)

Conversation between David Anderson and Claude (Cowork, model claude-fable-5-1).

Note on completeness: David's prompts are reproduced verbatim. The conversation was compacted once
(context limit) between the two sessions, so Claude's earlier messages are reproduced from notes and
from the delivered documents rather than verbatim; they are marked [condensed]. The final two Claude
messages are verbatim. Tool activity (reading, coding, compiling, refereeing) is summarised in
square brackets.

---

**David (session 1, opening):**

Look in the Dropbox folder attached, and start working on the problem described there. Remember, your
goal is a mathematical argument giving a solution or counterexample; focus on constructing proofs.
Write and run code *only* as needed. Use the included `mathmemo.sty` to write LaTeX notes when
headline results are achieved; keep intermediate notes in MD or whatever form you prefer. After
looking at the materials, check in with me if you see any ambiguity or have questions. Let's see if
you can knock down this open problem!

[Claude read `opac-009-handoff.md/.tex` and `mathmemo.sty`; read White, arXiv:0903.2831v2, in full
via the built-in browser (arXiv is blocked from the shell); checked the status of Gaetz–Meyer–Tam–
Wimberley–Yao–Zhu, arXiv:1409.4859, and found it withdrawn because its multi-factor
Littlewood–Richardson rule used sorted content (White's d-count) instead of the true coefficient c.]

**Claude (check-in) [condensed]:**

Two points of ambiguity before starting. (1) The handoff treats the families of Gaetz et al.
(A ∪ {(p,p)} extreme when A is, {(j,i),(j,i)} extreme, and their Conjectures 6.1/6.2 with the
reduction Lemma 6.3) as known; the paper is withdrawn and every one of those statements rests on the
flawed LR rule, so they are unproven. Should I treat them as open targets (so a rigorous proof of any
of them, or of any infinite family with repeated parts, counts as new), or as background facts?
(2) White's v2 paper (Theorem 4: non-nested ⇒ not extreme; Theorem 12: distinct parts ⇒ extreme;
Prop. 6 on dominance) looks solid — take it as given, or re-verify it first?

**David:**

Treat as open targets (Recommended).
Take White v2 as given (Recommended).

**David (mid-session):**

See also the more recent (2019) post here, by White:
https://realopacblog.wordpress.com/2019/09/22/the-schur-cone-and-the-cone-of-log-concavity/

[Claude read the post: same statements as the paper, no further families; problem open exactly as in
the handoff. Session 1 work: built an exact toolkit (`code/schurcone.py`: Schur products, White's
interlacing test, enumeration of nested multisets, exact LPs via cddlib); reproduced White's data for
N ≤ 12; explored several reduction ideas; found the specialisation argument; proved and refereed the
results below; wrote `memo/opac009-memo.pdf` (9 pp.), `notes/00-literature-status.md`,
`notes/01-findings.md`.]

**Claude (session-1 report) [condensed]:**

Three proven results, written up with `mathmemo.sty` and refereed by an independent pass.

Theorem ((p,p)-insertion). If s_{A''} is extreme in C²_{N''} then s_{A''∪{(p,p)}} is extreme in
C²_{N''+2p}. Proof: the ring homomorphism h_k ↦ exp(u(k)) with u(k)−u(k−1) strictly decreasing except
constant on {p, p+1} kills s_{(p,p)} and keeps every other generator positive, so every competitor in
a decomposition contains (p,p); cancel s_{(p,p)}. Corollary: every nested A that has distinct parts
after deleting its (v,v) pairs is extreme — a new infinite family with repeated parts, strictly
containing White's family (98 of the 127 nested multisets of weight 12). The complement (a value that
is the top of two or more pairs, or the bottom of two or more) is the remaining hard core.

Theorem (one-box LR formula). For μ = λ[i→j] (one box moved up), c^μ_B counts "valid paths"
J ⊆ [i, j−1] with i ∈ J and no violated pair of B; valid for repeated parts, specialising to White's
Lemma 15/Theorem 17 for distinct λ. This is the correct replacement for the flawed rule in Gaetz et
al.; c = d on one-box shapes but not on two-box shapes (their fatal error). Proof by Jacobi–Trudi
inclusion–exclusion, an "intruder path" description of the relevant Kostka numbers, and a
sign-reversing bijection.

Proposition. Explicit four-term separating functionals f = s_λ + s_{λ[v→v]} + s_{λ[x→y]} − s_{λ⁺}
prove extremality of {(j,i),(j,i)} (the withdrawn Theorem 5.1), {(v,x),(v,y)} and {(x,v),(y,v)}.

What remains: the hard core; computations show separators on [λ,λ⁺] with one-box shapes exist for all
N ≤ 12 except when a value is shared by three or more pairs, where two-box shapes ([λ,λ⁺⁺], White's
Conjecture 11) are needed.

**David (session 2):**

Good progress. Go ahead and attack the next stage of the problem, see if you can finish it off!

[Session 2 work. Several reduction ideas failed (restricted "avoidance" LPs, sign/initial-form
criteria, tropical/coarse criteria, a "splitting theorem" whose proof only controlled the lowest
level). The breakthrough was the Level Identity: perturbing the flat block of the specialisation and
comparing Taylor coefficients turns any decomposition into interval-count inequalities plus a
polynomial identity among products of arc forms ∏(δ_b − δ_{a+1}). The singleton-free case then
followed from a bottom-monomial argument; the context was compacted at this point. After compaction:
the one-singleton case was settled with a mixed monomial plus two dominance comparisons; the
combinatorial claim was checked for every nested A with N ≤ 18 (`oddmono.py`) and the level identity
on real decompositions of every non-nested A with N ≤ 9 (`verify_levelid.py`); the memo was
rewritten around the main theorem (old results moved to appendices), refereed by two independent
subagent passes (no mathematical errors; editorial fixes applied), and committed to the Dropbox
folder and the Project.]

**Claude (session-2 report, verbatim):**

The memo is finished and delivered. Here is the headline.

**Result.** White's conjecture (OPAC-009) is proved in full: for a multiset $A$ of partitions with at most two parts, $s_A$ is extreme in $\mathcal C^2_N$ if and only if $A$ is nested. The write-up is `memo/opac009-memo.pdf` in your Dropbox folder (tex alongside; also in the Project as `claude/opac009-memo.tex`), §§2–5, about six pages; the earlier session's results are now appendices. Two independent referee passes (subagents, adversarial, pen-and-paper) found no mathematical error, only editorial points, all of which I fixed.

**The proof in one paragraph.** Specialise $h_k\mapsto e^{u(k)}$ with $d(k)=u(k)-u(k-1)$ strictly decreasing except constant on a block $T=[x,y+1]$; this kills exactly the $s_{(a,b)}$ with $x\le b\le a\le y$ and keeps every other generator positive (that was the interval lemma from session 1). The new step is to perturb the flat block by $t\,e$, $e\in\mathbb R^T$, and expand in $t$: each $\psi_t(s_B)$ equals $t^{\,n_I(B)}\prod_{(a,b)\in B_I}(e_b-e_{a+1})\cdot R_B(t)$ with $R_B(0)>0$ independent of $e$. Applied to a decomposition $s_A=\sum c_Bs_B$ this gives (i) every $B$ in the support has at least as many pairs inside $[x,y]$ as $A$ does, for every interval simultaneously, and (ii) a polynomial identity $P_A=\sum c''_B P_{B}$ among products of the linear forms $\delta_b-\delta_{a+1}$ (taking $[x,y]=[\lambda_n,\lambda_1]$). The inequalities (i) confine the competitors so tightly (bottoms stochastically larger, tops stochastically smaller than $A$'s, interval by interval) that a single monomial does the rest: for singleton-free $A$, the monomial $\prod\delta_b$ has coefficient $1$ in $P_A$ and coefficient $0$ in every competitor's $P_B$ (a competitor producing it must have the same bottoms, then by weight the same tops, and a nested matching is determined by its tops and bottoms). With a singleton $(c)$, the pairs of $A$ split into those below $c$, the $(c,c)$'s, and those above $c$ (nestedness), competitors split the same way, and the mixed monomial $\prod_{\text{below}}\delta_{a+1}\prod_{\text{above}}\delta_b$ works, with the dominance condition $\phi(B)\trianglerighteq\lambda$ used twice to force the survivor to be $A$ itself. Dominance is essential there ($A=\{(5,4),(3),(2,1)\}$, $B=\{(4,4),(4),(2,1)\}$ survives everything else).

**Checks.** Beyond the refereeing: the level identity (both conclusions, all intervals) was verified on actual exact decompositions of every non-nested $A$ up to $N=9$; the monomial-coefficient claim was verified against every competitor of every nested $A$ up to $N=18$ (`code/oddmono.py`, `code/verify_levelid.py`, README updated).

**Caveats and what it does not do.** The reduction to nested competitors with $\phi(B)\trianglerighteq\phi(A)$ leans on White's Theorem 4 and Prop. 6 as given, per your instruction. The proof yields no explicit separating linear functional, so White's Conjecture 11 (separators on $[\lambda,\lambda^{++}]$) remains open as stated; the Taylor-coefficient functionals $L_{k,e}$ might be worth comparing with those intervals. I have not looked at whether the method reaches OPAC-010 ($k\ge3$); the specialisation of a $3\times3$ Jacobi–Trudi determinant no longer factors into a product of arc forms, so that would need a new idea.

Given the stakes, I'd suggest reading §§3–5 with a sceptical eye yourself before anything is circulated; the places I'd look hardest are Theorem 3.3 (the analytic step) and the two dominance comparisons in Step 4.

**David:**

Looks good -- well done! For my records, please prepare a minimally-formatted MD transcript of our
conversation. It must include all my prompts, and as much of your side as you're willing to provide
(without the note becoming gargantuan).

**Claude:**

[This file; saved as `notes/transcript-2026-09-09.md` in the Dropbox folder.]

---

Deliverables in `experiments/OPAC009_Fable`: `memo/opac009-memo.pdf` and `.tex` (main theorem §§2–5;
appendices: one-box LR formula, four-term certificates); `notes/00-literature-status.md`,
`notes/01-findings.md`, this transcript; `code/` with `README.md` (key scripts: `schurcone.py`,
`verify_levelid.py`, `oddmono.py`, `levelid.py`, `lstar.py`).
