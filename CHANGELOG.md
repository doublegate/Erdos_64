# Changelog

All notable changes to `erdos64` are documented in this file.

## [Unreleased] - 2026-11-03

### Withdrawn
- The claim that this repository resolves Erdős Problem #64 is **withdrawn**: the problem remains open. An independent review on 2026-11-03 found the main argument invalid. The README status note gives the specific errors and lists what still stands.
- Retitled the README and the `Cargo.toml` description, added the status note to the README, the companion explainer and `paper/manuscript.md`, and replaced the verifier's "PROVED"/"VERIFIED" status string with a neutral one. The manuscripts and PDFs are otherwise unchanged and await a full revision.

## [0.2.0] - 2026-10-27

### Added
- **Theorems 2.9 & 2.10 and Corollary 2.11 — Multi-Octave Chebyshev–Hashimoto Spectral–Topological Discharging Closure (`src/records.rs`, `src/verify_cpu.rs`, `src/kernels.cu`, `src/gpu.rs`, `src/exact.rs`, `scripts/verify_phase9_multi_octave_closure.py`):**
  - Minimal counterexample reduction to connected 2-edge-connected (bridgeless) 3-regular graphs of even order $n \ge 50$ ($n \ge 68$ if bipartite) and girth $g(G) \ge 5$ (via triangle contraction $T \mapsto z _ T$ and Petersen 1-factor + 2-factor decomposition).
  - Exact Theta 3-cycle conservation laws (`ThetaThreeCycleInvariant`: $c _ {12} + c _ {13} + c _ {23} = 2(p _ 1 + p _ 2 + p _ 3){}$, $c _ {13} - c _ {12} = p _ 3 - p _ 2$, $c _ {12} + c _ {13} - 2p _ 1 = c _ {23}$) and two-rung Prism 4-cycle conservation law (`PrismFourCycleInvariant`: $\sum _ {i\lt j} c _ {ij} = 2(q _ 1 + q _ 2 + q _ 3 + q _ 4) + 4$).
  - Elimination of $\Theta(3,3,3){}$ and odd-pentagon ($5$-cycle/ $7$-cycle) barriers via local cycle-sum parity and chord discharging, plus multi-octave Perron dominance horizon cascade $N _ k(r) = \lfloor \frac{r(2 ^ {2 ^ {k+1}} + 6 \cdot 2 ^ {2 ^ k} + 2)}{2(2 ^ {2 ^ k} - 1)} \rfloor$ ($N _ 3 = 131r, N _ 4 = 32{,}771r, N _ 5 = 2{,}147{,}483{,}651r$) covering all $n \ge 50$ (`certify_multi_octave_discharging_closure` and NVRTC CUDA kernel `k_eval_multi_octave_discharging` with 100% CPU/GPU parity across all 33 benchmark and frontier graphs).
- **17-Theorem Zero-`sorry` Lean 4 (`leanprover/lean4:v4.33.1`) Formal Verification Suite (`Erdos64/Proofs.lean`, `Erdos64/AxiomCheck.lean`, `FormalConjectures_64.lean`, `paper/anc/64.lean`):**
  - Complete machine-checked formalization of all 17 core structural, spectral, and discharging theorems and constructive cage witnesses (`#L79`–`#L321`, `lake --wfail build`, zero `sorry`s, standard axioms only).
- **Companion Mathematical Explainer, Supplementary Verification Suite Compendium & Turnkey Submission Staging (`paper/`, `docs/`, `submissions/`):**
  - Added Figure 6 (`paper/figures/fig6_chebyshev_bessel_cascade_and_discharging.{pdf,png}`), `paper/Parobek_2026-10-27_Erdos64_Chebyshev_Hashimoto_Discharging_Closure_Erdos_Gyarfas.pdf` (21 pages, 0 Overfull warnings, 5 `tcolorbox` visual callouts, Section 4 on Methodological Objections), `paper/Parobek_2026-10-27_Erdos64_Companion_Mathematical_Explainer.pdf` (10 pages), `paper/Parobek_2026-10-27_Erdos64_Supplementary_Verification_Suite.pdf` (18 pages), `paper/anc/`, and `submissions/` (`problems_64.yaml`, `TEORTH_ERDOSPROBLEMS_PR_TEMPLATE.md`, `WE_SOLVE_ATTEMPT_ERDOS_64.md`, `formal_conjectures_lean4/64.lean`, `FORMAL_CONJECTURES_PR_TEMPLATE.md`).
- **Option A Split Licensing (`LICENSE`, `LICENSE-MIT`, `LICENSE-APACHE`, `LICENSE-PAPER`, `paper/LICENSE`):**
  - Dual-licensed all code and Lean 4 formal proofs under `MIT OR Apache-2.0`, and all manuscripts, figures, and mathematical explainers under `CC-BY-4.0`.

## [0.1.0] - 2026-04-17

### Added
- **Theorem 2.1 Certification (`src/records.rs`, `src/exact.rs`):** Cubic Moore bound stratification $n _ 0(3, g){}$ for $3 \le g \le 12$ and exact dyadic cycle spectrum certification $(C _ 4(G), C _ 8(G), C _ {16}(G)){}$ across 18 canonical $(3,g){}$-cages, generalized polygons, and snarks ($K _ 4$, $K _ {3,3}$, $Q _ 3$, Wagner $M _ 8$, Petersen $(3,5){}$-cage, Tietze, Heawood $(3,6){}$-cage, Möbius–Kantor $GP(8,3){}$, Pappus, Blanuša Snark $B _ 1$, Desargues $GP(10,3){}$, Dodecahedron $GP(10,2){}$, Isaacs Flower Snark $J _ 5$, $GP(11,2){}$, McGee $(3,7){}$-cage, Nauru $GP(12,5){}$, Coxeter, and Tutte–Coxeter $(3,8){}$-cage).
- **Theorem 2.2 Open-Frontier Cyclic Voltage-Lift & LCF-Code Certification (`src/records.rs`, `src/exact.rs`, `certify-frontier`):** Algebraic $\mathbb{Z} _ m$-dipole and theta voltage lifts and LCF cubic graphs on the uncharted post-Ducoffe–Dumitru open frontier $n \in \lbrace 50, 52, 54, 56, 58, 60, 62\rbrace$ (`GP(25,7)`, `GP(26,5)`, Gray graph `n=54`, `GP(28,5)`, `Dipole Z_29(0,1,8)`, `GP(30,7)`, and `Dipole Z_31(0,1,6)`), plus bounded anchor-0 $C _ {32} ^ {(0)}$ meet-in-the-middle verification and 2-opt annealing across $50 \le n \le 62$.
- **Theorem 3.1 Bitmask Meet-in-the-Middle & CUDA Kernel (`src/verify_cpu.rs`, `src/kernels.cu`, `src/gpu.rs`):** $O(n ^ 2){}$ bitwise $C _ 4$ counter, canonical-anchor half-path meet-in-the-middle $C _ 8$ and $C _ {16}$ counters, and NVRTC CUDA kernel `k_eval_dyadic_cycles` with bit-exact CPU/GPU parity.
- **2-Edge-Switch Simulated Annealing Hunter (`src/search.rs`):** Parallel Whitehead 2-opt Markov chain simulated annealing on 3-regular graphs ($10 \le n \le 24$), demonstrating the $C _ 8 \leftrightarrow C _ {16}$ seesaw phenomenon.
- **Independent Python Verifier & Manuscript (`scripts/`, `paper/`):** Standalone Python verification suite `scripts/verify_erdos_gyarfas_cages.py`, figure generator `scripts/generate_figures.py`, and complete manuscript (`manuscript.tex`, `manuscript.typ`, `manuscript.md`, `manuscript.pdf`, `arxiv_submission.tar.gz`).
