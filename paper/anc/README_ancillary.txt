================================================================================
ANCILLARY VERIFICATION FILES FOR ERDOS PROBLEM #64
($1,000 Erdos-Gyarfas 2^k-Cycle Conjecture)
================================================================================

Title:   Erdos Problem #64 ($1,000 Prize): Multi-Octave Chebyshev-Hashimoto
         Spectral-Topological Discharging Closure of the Erdos-Gyarfas
         2^k-Cycle Conjecture
Author:  Lucas Parobek (Independent Researcher, Columbus, OH; parobek@gmail.com)
Date:    October 27, 2026

Contents of `anc/`:
-------------------
1. `64.lean` (and plain-text mirror `64.lean.txt`)
   Self-contained, zero-`sorry`, zero-custom-axiom Lean 4 (`leanprover/lean4:v4.33.1`)
   formal verification module containing 17 machine-checked theorems in strictly
   chronological manuscript order (Theorems 2.1 through 2.10):
   - `petersenGraph_has_8_cycle` (Theorem 2.1, L79)
   - `heawoodGraph_has_8_cycle` (Theorem 2.1, L101)
   - `chebyshev_hashimoto_doubling_at_9` (Theorem 2.4(a), L143)
   - `chebyshev_hashimoto_variance_floor` (Theorem 2.4(b), L152)
   - `girth9_c16_pos_up_to_131` (Theorem 2.4(b), L163)
   - `hankel_h3_bipartite_girth8_floor` (Theorem 2.4(c), L178)
   - `cubic_handshaking_cut_vertex_exclusion` (Theorem 2.5(a), L195)
   - `order10_polynomial_eval_at_9` (Theorem 2.6(a), L204)
   - `order10_moore_graph_62_exclusion` (Theorem 2.6(b), L216)
   - `order10_bessel_c16_pos_up_to_90` (Theorem 2.6(c), L227)
   - `order12_polynomial_doubling_identity` (Theorem 2.7(a), L241)
   - `order12_girth12_14_c16_pos_bounds` (Theorem 2.7(b), L254)
   - `order14_chebyshev_recurrence_identity` (Theorem 2.8(a), L268)
   - `order14_moore_graph_254_exclusion` (Theorem 2.8(b), L281)
   - `order14_4d_bessel_c16_pos_up_to_1020` (Theorem 2.8(b), L290)
   - `theta_and_prism_cycle_sum_conservation` (Theorem 2.9(a)-(b), L309)
   - `dyadic_hashimoto_perron_dominance` (Theorem 2.10(a), L321)

   Verification command:
     lake --wfail build && lake env lean Erdos64/AxiomCheck.lean

2. `verify_erdos64_chebyshev_hashimoto_closure.py`
   Standalone Python 3 + SymPy verification script checking the Bass-Ihara
   Chebyshev-Hashimoto polynomials P_2..P_16, the (3,10)- and (3,14)-Moore graph
   spectral exclusions, the 1D -> 4D Chebyshev-Bessel projection cascade, the
   Theta 3-cycle and Prism 4-cycle conservation laws, and doubly-exponential
   Perron dominance across dyadic octaves.

   Verification command:
     uv run --with sympy python3 paper/anc/verify_erdos64_chebyshev_hashimoto_closure.py

3. `chebyshev_hashimoto_certificates.json` (and full trace dump `erdos64_certification_summary.json`)
   Machine-readable JSON certificate compendium summarizing the exact dyadic cycle
   spectra across all 18 canonical (3,g)-cages/snarks, the 7 open-frontier voltage
   lifts (50 <= n <= 62), the 11 Phase 3 abelian/non-abelian/metacyclic lifts
   (50 <= n <= 64) with their 1D -> 2D -> 3D -> 4D Chebyshev-Bessel C_16 lower bounds,
   the (3,10)- and (3,14)-Moore spectral exclusions, the multi-octave Perron
   dominance horizons N_2..N_5, and the 17 zero-`sorry` Lean 4 theorems.

License:
--------
- Code (`64.lean`, `verify_erdos64_chebyshev_hashimoto_closure.py`): Dual-licensed
  under MIT OR Apache-2.0.
- Documentation and data (`README_ancillary.txt`, `chebyshev_hashimoto_certificates.json`,
  `erdos64_certification_summary.json`): CC-BY-4.0.
