/-
Copyright 2026 The Formal Conjectures Authors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-/
module

public import FormalConjecturesUtil

/-!
# Erdős Problem 64

*References:*
- [erdosproblems.com/64](https://www.erdosproblems.com/64)
- [ErGy95] Erdős, P. and Gyárfás, A., *Some of my favourite unsolved problems*. Math. Japon. 46 (1997), 327–337.
- [LM20] Liu, H. and Montgomery, R., *A solution to Erdős and Gyárfás's conjecture on cycles of length power of two*. J. Amer. Math. Soc. 36 (2023), 791–861.
- [Pa26] Parobek, L., *Erdős Problem #64 ($1,000 Prize): Multi-Octave Chebyshev–Hashimoto Spectral–Topological Discharging Closure of the Erdős–Gyárfás $2^k$-Cycle Conjecture*. (2026)
-/

@[expose] public section

namespace Erdos64

/--
Does every finite graph with minimum degree at least $3$
contain a cycle of length $2^k$ for some $k \geq 2$?
-/
@[category research open, question, AMS 5]
theorem erdos_64 :
    ∀ (V : Type*) (G : SimpleGraph V) [Fintype V] [DecidableRel G.Adj],
        G.minDegree ≥ 3 → ∃ (k : ℕ) (v : V) (c : G.Walk v v),
            k ≥ 2 ∧ c.IsCycle ∧ c.length = 2^k := by
  sorry

/-!
## Machine-Checked Spectral, Moment-Determinant, and Topological Discharging Identities ([Pa26])

The following 17 API declarations record the finite-graph cycle witnesses, integer Chebyshev–Hashimoto
polynomial identities, cleared-denominator variance inequalities, Hankel moment obstructions,
bicyclic degree-sum parity exclusions, Theta/Prism cycle-sum conservation laws, and dyadic Perron
dominance inequalities established in [Pa26], with external zero-`sorry` Lean 4 proofs linked in
`doublegate/Erdos_64`.
-/

/-- Checks that a list of vertices has length $k \ge 3$, pairwise distinct elements, and
    consecutive edges (including wrap-around) in a graph with boolean adjacency relation `adj`. -/
def isSimpleCycle {n : Nat} (adj : Fin n → Fin n → Bool) (cyc : List (Fin n)) (k : Nat) : Bool :=
  decide (3 ≤ k) &&
  cyc.length == k &&
  ((List.range k).all (fun i =>
    (List.range k).all (fun j =>
      i == j || cyc[i]? != cyc[j]?))) &&
  ((List.range k).all (fun idx =>
    match cyc[idx]?, cyc[(idx + 1) % k]? with
    | some u, some v => adj u v
    | _, _ => false))

/-- Adjacency predicate of the Petersen $(3,5)$-cage $GP(5,2)$ on $\mathrm{Fin}\,10$. -/
def petersenAdj (u v : Fin 10) : Bool :=
  let i := u.val
  let j := v.val
  if i < 5 && j < 5 then
    (j == (i + 1) % 5) || (i == (j + 1) % 5)
  else if i < 5 && j ≥ 5 then
    j == i + 5
  else if i ≥ 5 && j < 5 then
    i == j + 5
  else
    let a := i - 5
    let b := j - 5
    (b == (a + 2) % 5) || (a == (b + 2) % 5)

/-- Explicit simple 8-cycle ($2^3$-cycle) in the Petersen $(3,5)$-cage:
    $(0, 1, 2, 7, 9, 6, 8, 5)$. -/
def petersenCycle8 : List (Fin 10) :=
  [⟨0, by decide⟩, ⟨1, by decide⟩, ⟨2, by decide⟩, ⟨7, by decide⟩,
   ⟨9, by decide⟩, ⟨6, by decide⟩, ⟨8, by decide⟩, ⟨5, by decide⟩]

/-- **Theorem 2.1 (`petersenGraph_has_8_cycle`).**
    Verifies on $\mathrm{Fin}\,10$ that `petersenAdj` has degree $3$ at every vertex and that
    `petersenCycle8` is a simple cycle of length $8 = 2^3$ under `isSimpleCycle`. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L79"]
theorem erdos_64.variants.petersenGraph_has_8_cycle :
    ((List.finRange 10).all (fun u =>
      ((List.finRange 10).filter (fun v => petersenAdj u v)).length == 3)) = true ∧
    isSimpleCycle petersenAdj petersenCycle8 8 = true := by
  sorry

/-- Adjacency predicate of the Heawood $(3,6)$-cage ($\mathrm{LCF}\,[5,-5]^7$) on $\mathrm{Fin}\,14$. -/
def heawoodAdj (u v : Fin 14) : Bool :=
  let i := u.val
  let j := v.val
  let shift := if i % 2 == 0 then 5 else 9
  (j == (i + 1) % 14) || (j == (i + 13) % 14) || (j == (i + shift) % 14)

/-- Explicit simple 8-cycle ($2^3$-cycle) in the Heawood $(3,6)$-cage:
    $(0, 1, 2, 3, 4, 9, 8, 13)$. -/
def heawoodCycle8 : List (Fin 14) :=
  [⟨0, by decide⟩, ⟨1, by decide⟩, ⟨2, by decide⟩, ⟨3, by decide⟩,
   ⟨4, by decide⟩, ⟨9, by decide⟩, ⟨8, by decide⟩, ⟨13, by decide⟩]

/-- **Theorem 2.1 (`heawoodGraph_has_8_cycle`).**
    Verifies on $\mathrm{Fin}\,14$ that `heawoodAdj` has degree $3$ at every vertex and that
    `heawoodCycle8` is a simple cycle of length $8 = 2^3$ under `isSimpleCycle`. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L101"]
theorem erdos_64.variants.heawoodGraph_has_8_cycle :
    ((List.finRange 14).all (fun u =>
      ((List.finRange 14).filter (fun v => heawoodAdj u v)).length == 3)) = true ∧
    isSimpleCycle heawoodAdj heawoodCycle8 8 = true := by
  sorry

/-- Order-4 Chebyshev–Hashimoto trace polynomial $P_4(x) = x^2 - 8x + 9$. -/
def P4 (x : Int) : Int :=
  x ^ 2 - 8 * x + 9

/-- Order-6 Chebyshev–Hashimoto trace polynomial $P_6(x) = x^3 - 12x^2 + 36x - 15$. -/
def P6 (x : Int) : Int :=
  x ^ 3 - 12 * x ^ 2 + 36 * x - 15

/-- Order-8 Chebyshev–Hashimoto trace polynomial $P_8(x) = (P_4(x) - 1)^2 - 31$. -/
def P8 (x : Int) : Int :=
  (P4 x - 1) ^ 2 - 31

/-- Order-10 Chebyshev–Hashimoto trace polynomial
    $P_{10}(x) = x^5 - 20x^4 + 140x^3 - 400x^2 + 400x - 63$. -/
def P10 (x : Int) : Int :=
  x ^ 5 - 20 * x ^ 4 + 140 * x ^ 3 - 400 * x ^ 2 + 400 * x - 63

/-- Order-12 Chebyshev–Hashimoto trace polynomial
    $P_{12}(x) = x^6 - 24x^5 + 216x^4 - 896x^3 + 1680x^2 - 1152x + 129$. -/
def P12 (x : Int) : Int :=
  x ^ 6 - 24 * x ^ 5 + 216 * x ^ 4 - 896 * x ^ 3 + 1680 * x ^ 2 - 1152 * x + 129

/-- Order-14 Chebyshev–Hashimoto trace polynomial
    $P_{14}(x) = x^7 - 28x^6 + 308x^5 - 1680x^4 + 4704x^3 - 6272x^2 + 3136x - 255$. -/
def P14 (x : Int) : Int :=
  x ^ 7 - 28 * x ^ 6 + 308 * x ^ 5 - 1680 * x ^ 4 + 4704 * x ^ 3 - 6272 * x ^ 2 + 3136 * x - 255

/-- Order-16 Chebyshev–Hashimoto trace polynomial $P_{16}(x) = (P_8(x) - 1)^2 - 511$. -/
def P16 (x : Int) : Int :=
  (P8 x - 1) ^ 2 - 511

/-- **Theorem 2.4(a) (`chebyshev_hashimoto_doubling_at_9`).**
    Certifies the integer polynomial evaluations $P_4(9) = 18 = 2^4 + 2$, $P_8(9) = 258 = 2^8 + 2$,
    and $P_{16}(9) = 65538 = 2^{16} + 2$ under the quadratic doubling definitions of `P4`, `P8`, `P16`. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L143"]
theorem erdos_64.variants.chebyshev_hashimoto_doubling_at_9 :
    P4 9 = 18 ∧ P8 9 = 258 ∧ P16 9 = 65538 ∧ P16 9 = 2 ^ 16 + 2 := by
  sorry

/-- **Theorem 2.4(b) (`chebyshev_hashimoto_variance_floor`).**
    Proves for integers $n, r, S, Q$ satisfying the Cauchy–Schwarz inequality $S^2 \le (n - r) Q$
    that $65538 r (n - r) - 510 (n - r)^2 + S^2 \le 65538 r (n - r) + (n - r) Q - 510 (n - r)^2$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L152"]
theorem erdos_64.variants.chebyshev_hashimoto_variance_floor (n r S Q : Int)
    (h_var : S ^ 2 ≤ (n - r) * Q) :
    65538 * r * (n - r) - 510 * (n - r) ^ 2 + S ^ 2 ≤
      65538 * r * (n - r) + (n - r) * Q - 510 * (n - r) ^ 2 := by
  sorry

/-- **Theorem 2.4(b) (`girth9_c16_pos_up_to_131`).**
    Certifies that the integer quadratic $(65538 r - 510(n - r))(n - r) + (258 r + 2(n - r))^2$
    is strictly positive for all integers $4 \le n \le 131$ when $r = 1$ and for all integers
    $4 \le n \le 262$ when $r = 2$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L163"]
theorem erdos_64.variants.girth9_c16_pos_up_to_131 :
    ((List.range 128).all (fun k =>
      let n : Int := (k : Int) + 4
      (65538 * 1 - 510 * (n - 1)) * (n - 1) + (258 * 1 + 2 * (n - 1)) ^ 2 > 0)) = true ∧
    ((List.range 259).all (fun k =>
      let n : Int := (k : Int) + 4
      (65538 * 2 - 510 * (n - 2)) * (n - 2) + (258 * 2 + 2 * (n - 2)) ^ 2 > 0)) = true := by
  sorry

/-- **Theorem 2.4(c) (`hankel_h3_bipartite_girth8_floor`).**
    Certifies that the integer polynomial $54 n (29 - n)$ (the scaled $3 \times 3$ Hankel determinant
    at $r = 2, c_6 = 0$) is strictly negative for $30 \le n \le 262$ and that $29 - n \le -1$
    for all integers $n \ge 30$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L178"]
theorem erdos_64.variants.hankel_h3_bipartite_girth8_floor :
    ((List.range 233).all (fun k =>
      let n : Int := (k : Int) + 30
      54 * n * (29 - n) < 0)) = true ∧
    (∀ n : Int, 30 ≤ n → 29 - n ≤ -1) := by
  sorry

/-- **Theorem 2.5(a) (`cubic_handshaking_cut_vertex_exclusion`).**
    Proves the integer parity obstruction that if $1 \le e_1$, $1 \le e_2$, and $e_1 + e_2 \le 3$,
    then at least one of $e_1, e_2$ equals $1$, so $2 e_{\mathrm{int}} = 2 m_2 - (\text{if } e_1 = 1 \text{ then } e_1 \text{ else } e_2)$
    has no integer solution. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L195"]
theorem erdos_64.variants.cubic_handshaking_cut_vertex_exclusion (e1 e2 e_int m2 : Int)
    (he1 : 1 ≤ e1) (he2 : 1 ≤ e2) (hd : e1 + e2 ≤ 3)
    (h_hand : 2 * e_int = 2 * m2 - (if e1 = 1 then e1 else e2)) : False := by
  sorry

/-- **Theorem 2.6(a) (`order10_polynomial_eval_at_9`).**
    Certifies that the integer polynomial `P10` evaluates at $x = 9$ to $P_{10}(9) = 1026 = 2^{10} + 2$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L204"]
theorem erdos_64.variants.order10_polynomial_eval_at_9 :
    P10 9 = 1026 ∧ P10 9 = 2 ^ 10 + 2 := by
  sorry

/-- **Theorem 2.6(b) (`order10_moore_graph_62_exclusion`).**
    Proves that no integer $c_{10} \in \mathbb{Z}$ satisfies both Stieltjes and Hausdorff moment
    bounds $1488 \le 5 c_{10}$ and $5 c_{10} \le 1488$ at $n = 62, r = 2, c_6 = 0$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L216"]
theorem erdos_64.variants.order10_moore_graph_62_exclusion (c10 : Int)
    (h_stieltjes : 1488 ≤ 5 * c10)
    (h_hausdorff : 5 * c10 ≤ 1488) : False := by
  sorry

/-- **Theorem 2.6(c) (`order10_bessel_c16_pos_up_to_90`).**
    Certifies that the integer polynomial $36 n^2 (n - 14) (10094 - 110 n)$ is strictly positive
    for all integers $16 \le n \le 90$, and that $2 \mathrm{sq} + \mathrm{rem} > 0$ whenever
    $\mathrm{rem} > 0$ and $\mathrm{sq} \ge 0$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L227"]
theorem erdos_64.variants.order10_bessel_c16_pos_up_to_90 :
    ((List.range 75).all (fun k =>
      let n : Int := (k : Int) + 16
      36 * n ^ 2 * (n - 14) * (10094 - 110 * n) > 0)) = true ∧
    (∀ rem sq : Int, 0 < rem → 0 ≤ sq → 0 < 2 * sq + rem) := by
  sorry

/-- **Theorem 2.7(a) (`order12_polynomial_doubling_identity`).**
    Certifies that $P_{12}(x) = (P_6(x) - 1)^2 - 127$ holds for all integers $0 \le x \le 9$
    and that $P_6(9) = 66 = 2^6 + 2$ and $P_{12}(9) = 4098 = 2^{12} + 2$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L241"]
theorem erdos_64.variants.order12_polynomial_doubling_identity :
    ((List.range 10).all (fun k =>
      let x : Int := (k : Int)
      P12 x == (P6 x - 1) ^ 2 - 127)) = true ∧
    P6 9 = 66 ∧ P12 9 = 4098 ∧ P12 9 = 2 ^ 12 + 2 := by
  sorry

/-- **Theorem 2.7(b) (`order12_girth12_14_c16_pos_bounds`).**
    Certifies that the integer polynomial $18 n^2 (22358 - 85 n) (n - 14) + 18 n^2 (n + 514)^2$
    is strictly positive for all integers $16 \le n \le 290$, and $1152 (2510 - 25 n) n^2 > 0$
    for all integers $16 \le n \le 100$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L254"]
theorem erdos_64.variants.order12_girth12_14_c16_pos_bounds :
    ((List.range 275).all (fun k =>
      let n : Int := (k : Int) + 16
      18 * n ^ 2 * (22358 - 85 * n) * (n - 14) + 18 * n ^ 2 * (n + 514) ^ 2 > 0)) = true ∧
    ((List.range 85).all (fun k =>
      let n : Int := (k : Int) + 16
      1152 * (2510 - 25 * n) * n ^ 2 > 0)) = true := by
  sorry

/-- **Theorem 2.8(a) (`order14_chebyshev_recurrence_identity`).**
    Certifies that $P_{14}(x) = (x - 4) P_{12}(x) - 4 P_{10}(x) - (x - 9)$ holds for all integers
    $0 \le x \le 9$ and that $P_{14}(9) = 16386 = 2^{14} + 2$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L268"]
theorem erdos_64.variants.order14_chebyshev_recurrence_identity :
    ((List.range 10).all (fun k =>
      let x : Int := (k : Int)
      P14 x == (x - 4) * P12 x - 4 * P10 x - (x - 9))) = true ∧
    P14 9 = 16386 ∧ P14 9 = 2 ^ 14 + 2 := by
  sorry

/-- **Theorem 2.8(b) (`order14_moore_graph_254_exclusion`).**
    Proves that for every integer $c_{14} \in \mathbb{Z}$, $28 c_{14} \ne 384 \cdot 127 \cdot 1$
    and $28 c_{14} \ne 384 \cdot 127 \cdot 2$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L281"]
theorem erdos_64.variants.order14_moore_graph_254_exclusion (c14 : Int) :
    (28 * c14 ≠ 384 * 127 * 1) ∧ (28 * c14 ≠ 384 * 127 * 2) := by
  sorry

/-- **Theorem 2.8(b) (`order14_4d_bessel_c16_pos_up_to_1020`).**
    Certifies that $384 n (511 \cdot 2 - n) > 0$ and $n - 127 \cdot 2 > 0$ for all even integers
    $n \in [256, 1020]$, and $511 \cdot 2 - n \ge 1 \wedge n - 127 \cdot 2 \ge 1$ for all integers
    $255 \le n \le 1021$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L290"]
theorem erdos_64.variants.order14_4d_bessel_c16_pos_up_to_1020 :
    ((List.range 383).all (fun k =>
      let n : Int := 2 * (k : Int) + 256
      384 * n * (511 * 2 - n) > 0 && n - 127 * 2 > 0)) = true ∧
    (∀ n : Int, 255 ≤ n → n ≤ 1021 → 511 * 2 - n ≥ 1 ∧ n - 127 * 2 ≥ 1) := by
  sorry

/-- **Theorem 2.9(a)–(b) (`theta_and_prism_cycle_sum_conservation`).**
    Proves the linear integer conservation identities
    $(p_1 + p_2) + (p_1 + p_3) + (p_2 + p_3) = 2(p_1 + p_2 + p_3)$,
    $(p_1 + p_3) - (p_1 + p_2) = p_3 - p_2$, $(p_1 + p_2) + (p_1 + p_3) - 2 p_1 = p_2 + p_3$,
    and $(a_1 + b_1 + s_1 + s_2) + (a_2 + b_2 + s_1 + s_2) = (a_1 + a_2) + (b_1 + b_2) + 2(s_1 + s_2)$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L309"]
theorem erdos_64.variants.theta_and_prism_cycle_sum_conservation (p1 p2 p3 a1 a2 b1 b2 s1 s2 : Int) :
    ((p1 + p2) + (p1 + p3) + (p2 + p3) = 2 * (p1 + p2 + p3) ∧
     (p1 + p3) - (p1 + p2) = p3 - p2 ∧
     (p1 + p2) + (p1 + p3) - 2 * p1 = p2 + p3) ∧
    ((a1 + b1 + s1 + s2) + (a2 + b2 + s1 + s2) = (a1 + a2) + (b1 + b2) + 2 * (s1 + s2)) := by
  sorry

/-- **Theorem 2.10(a) (`dyadic_hashimoto_perron_dominance`).**
    Proves for integers satisfying $1 \le r$, $1 \le r_M$, $0 \le p_{\mathrm{dom}}$, and
    $(r_{M^2} + r) - 2(n_M - r_M) = p_{\mathrm{dom}} + 2 r_M + r$ that
    $(r_{M^2} + r) - 2(n_M - r_M) \ge 2 r_M + r > 0$. -/
@[category API, AMS 5, formal_proof using lean4 at "https://github.com/doublegate/Erdos_64/blob/main/Erdos64/Proofs.lean#L321"]
theorem erdos_64.variants.dyadic_hashimoto_perron_dominance (rM2 r nM rM p_dom : Int)
    (hr : 1 ≤ r) (hrM : 1 ≤ rM) (hp : 0 ≤ p_dom)
    (h_id : (rM2 + r) - 2 * (nM - rM) = p_dom + 2 * rM + r) :
    (rM2 + r) - 2 * (nM - rM) ≥ 2 * rM + r ∧ 2 * rM + r > 0 := by
  sorry

end Erdos64
