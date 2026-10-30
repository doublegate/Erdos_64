//! ============================================================================
//! `exact.rs` — Exact cycle-spectrum certification for Theorem 2.1 across the
//! canonical `(3,g)`-Cage and Snark library.
//! ============================================================================

use rayon::prelude::*;

use crate::records::{
    build_canonical_cage_and_snark_library, build_open_frontier_voltage_and_lcf_library,
    build_phase3_commutator_frontier_library, build_zm_dipole_voltage_lift, moore_bound_cubic,
    CycleSpectrumRecord, CUBIC_CAGE_BOUNDS,
};
use crate::search::{run_hard_c4_c8_zero_c16_tabu_search, C16TabuSearchRecord};
use crate::verify_cpu::{
    certify_b10_hankel_and_dual_theta_pinch, certify_b12_deg2_bessel_pinch,
    certify_b14_deg3_bessel_pinch, certify_chebyshev_doubling_and_hankel_bounds,
    certify_free_abelian_commutator_floor, certify_hashimoto_trace_identity,
    certify_multi_octave_discharging_closure, certify_n16_topological_decomposition,
    compute_dyadic_spectrum_fast, compute_full_cycle_spectrum, B10HankelPinchCertificate,
    B12Deg2BesselCertificate, B14Deg3BesselCertificate, FreeAbelianCommutatorCertificate,
    HashimotoTraceCertificate, MultiOctaveDischargingCertificate,
    SpectralDoublingHankelCertificate, TopologicalN16Certificate,
};

/// Certifies all 18 canonical cubic cages, generalized polygons, and snarks.
/// For `n <= 24`, computes the full cycle-length spectrum `C_l(G)` for every `3 <= l <= n`;
/// for `n in [26, 30]` (Coxeter and Tutte–Coxeter `(3,8)`-cage), computes exact `C_4, C_8, C_16`
/// via bitmask meet-in-the-middle plus full cycle spectrum up to `max_full_n`.
pub fn certify_canonical_library(max_full_n: usize) -> Vec<CycleSpectrumRecord> {
    let lib = build_canonical_cage_and_snark_library();
    lib.par_iter()
        .map(|g| {
            if g.n <= max_full_n {
                compute_full_cycle_spectrum(g)
            } else {
                compute_dyadic_spectrum_fast(g)
            }
        })
        .collect()
}

/// Certifies Theorem 2.2 across the 7-graph Open-Frontier (`50 <= n <= 62`)
/// Cyclic Voltage-Lift and LCF-Code Cubic Graph Library.
pub fn certify_open_frontier_library() -> Vec<CycleSpectrumRecord> {
    let lib = build_open_frontier_voltage_and_lcf_library();
    lib.par_iter().map(compute_dyadic_spectrum_fast).collect()
}

/// Prints a formatted ASCII table of Theorem 2.2's Open-Frontier (`50 <= n <= 62`)
/// cyclic voltage-lift and LCF-code dyadic cycle certificates.
pub fn print_open_frontier_report(records: &[CycleSpectrumRecord]) {
    println!("========================================================================================================================");
    println!("  Theorem 2.2: Cyclic Voltage-Lift & LCF-Code Dyadic Cycle Spectra on the Open Frontier (50 <= n <= 62)");
    println!("========================================================================================================================");
    println!(
        "{:<30} {:>4} {:>4} {:>5} {:>5} {:>6} {:>8} {:>10} {:>10} {:>12}  {:<18}",
        "Frontier Graph",
        "n",
        "m",
        "Girth",
        "Bip?",
        "C_4",
        "C_8",
        "C_16",
        "C_32^(0)",
        "Total 2^k",
        "Dyadic Lengths"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for r in records {
        let bip = if r.is_bipartite { "yes" } else { "no" };
        let dyadic_str = format!("{:?}", r.dyadic_lengths_present);
        println!(
            "{:<30} {:>4} {:>4} {:>5} {:>5} {:>6} {:>8} {:>10} {:>10} {:>12}  {:<18}",
            r.name,
            r.n,
            r.m,
            r.girth,
            bip,
            r.c4,
            r.c8,
            r.c16,
            r.c32,
            r.total_dyadic_cycles,
            dyadic_str
        );
    }
    println!("========================================================================================================================");
}

/// Prints a formatted ASCII table of the Moore bounds `n_0(3, g)`, exact `(3, g)`-cage orders,
/// and the certified cycle-length spectra of the canonical cage & snark library.
pub fn print_cage_certification_report(records: &[CycleSpectrumRecord]) {
    println!("==========================================================================================");
    println!("  Theorem 2.1: Cubic (3, g)-Cage Moore Bounds & Smallest Dyadic Cycle Lengths");
    println!("==========================================================================================");
    println!(
        "{:<6} {:>12} {:>12} {:<44} {:>10}",
        "Girth", "Moore n0(3,g)", "Cage f(3,g)", "Canonical (3,g)-Cage", "Min 2^k"
    );
    println!("------------------------------------------------------------------------------------------");
    for b in CUBIC_CAGE_BOUNDS {
        let m = moore_bound_cubic(b.girth);
        println!(
            "{:<6} {:>12} {:>12} {:<44} {:>10}",
            b.girth, m, b.exact_cage_order, b.cage_name, b.smallest_dyadic_cycle_length
        );
    }
    println!("==========================================================================================\n");

    println!("================================================================================================================");
    println!("  Theorem 2.1 & 3.1: Exact Dyadic Cycle Multiplicities Across Canonical Cages, Generalized Polygons & Snarks");
    println!("================================================================================================================");
    println!(
        "{:<26} {:>4} {:>4} {:>5} {:>5} {:>8} {:>8} {:>10} {:>12}  {:<18}",
        "Graph", "n", "m", "Girth", "Bip?", "C_4", "C_8", "C_16", "Total 2^k", "Dyadic Lengths"
    );
    println!("----------------------------------------------------------------------------------------------------------------");
    for r in records {
        let bip = if r.is_bipartite { "yes" } else { "no" };
        let dyadic_str = format!("{:?}", r.dyadic_lengths_present);
        println!(
            "{:<26} {:>4} {:>4} {:>5} {:>5} {:>8} {:>8} {:>10} {:>12}  {:<18}",
            r.name, r.n, r.m, r.girth, bip, r.c4, r.c8, r.c16, r.total_dyadic_cycles, dyadic_str
        );
    }
    println!("================================================================================================================");
}

/// Complete Phase 3–7 certification bundle for Theorems 2.3, 2.4, 2.5, 2.6, 2.7, and 3.2.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Phase3CertificationBundle {
    pub free_abelian_certificate: FreeAbelianCommutatorCertificate,
    pub commutator_frontier_records: Vec<CycleSpectrumRecord>,
    pub hashimoto_trace_certificates: Vec<HashimotoTraceCertificate>,
    pub spectral_doubling_certificates: Vec<SpectralDoublingHankelCertificate>,
    pub topological_n16_certificates: Vec<TopologicalN16Certificate>,
    pub b10_hankel_certificates: Vec<B10HankelPinchCertificate>,
    pub b12_deg2_bessel_certificates: Vec<B12Deg2BesselCertificate>,
    pub b14_deg3_bessel_certificates: Vec<B14Deg3BesselCertificate>,
    pub multi_octave_discharging_certificates: Vec<MultiOctaveDischargingCertificate>,
    pub tabu_c16_records: Vec<C16TabuSearchRecord>,
}

/// Certifies Theorems 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, and 3.2 across the free abelian dipole cover, the 10-graph
/// commutator-breaking frontier library (`50 <= n <= 64`), the Hashimoto non-backtracking
/// spectral trace identities, the Chebyshev doubling & Hankel moment bounds, the exact
/// topological `N_16(G)` decomposition, the `Tr(B^{10}) == 20 C_{10}(G)` & 3x3 Hankel / Bessel
/// quadratic bounds, the `Tr(B^{12})` topological decomposition & 3D Bessel projection bounds,
/// the `Tr(B^{14})` topological decomposition & 4D Bessel projection bounds,
/// and the Hard-`C_4=C_8=0` `C_16`-Minimizing Tabu Search.
pub fn certify_phase3_bundle(
    restarts: usize,
    steps_per_restart: usize,
) -> Phase3CertificationBundle {
    let free_abelian_certificate = certify_free_abelian_commutator_floor();
    let p3_lib = build_phase3_commutator_frontier_library();
    let commutator_frontier_records: Vec<CycleSpectrumRecord> = p3_lib
        .par_iter()
        .map(compute_dyadic_spectrum_fast)
        .collect();

    let mut hash_graphs = Vec::new();
    let cages = build_canonical_cage_and_snark_library();
    for idx in [4usize, 6, 9, 14, 17] {
        // Petersen, Heawood, Blanusa B_1, McGee, Tutte-Coxeter
        hash_graphs.push(cages[idx].clone());
    }
    let f_lib = build_open_frontier_voltage_and_lcf_library();
    hash_graphs.push(f_lib[6].clone()); // Dipole Lift Z_31(0,1,6) [n=62]
    for g in &p3_lib {
        hash_graphs.push(g.clone());
    }
    let hashimoto_trace_certificates: Vec<HashimotoTraceCertificate> = hash_graphs
        .par_iter()
        .map(certify_hashimoto_trace_identity)
        .collect();
    let spectral_doubling_certificates: Vec<SpectralDoublingHankelCertificate> = hash_graphs
        .par_iter()
        .map(certify_chebyshev_doubling_and_hankel_bounds)
        .collect();
    let topological_n16_certificates: Vec<TopologicalN16Certificate> = hash_graphs
        .par_iter()
        .map(certify_n16_topological_decomposition)
        .collect();
    let b10_hankel_certificates: Vec<B10HankelPinchCertificate> = hash_graphs
        .par_iter()
        .map(certify_b10_hankel_and_dual_theta_pinch)
        .collect();
    let b12_deg2_bessel_certificates: Vec<B12Deg2BesselCertificate> = hash_graphs
        .par_iter()
        .map(certify_b12_deg2_bessel_pinch)
        .collect();
    let b14_deg3_bessel_certificates: Vec<B14Deg3BesselCertificate> = hash_graphs
        .par_iter()
        .map(certify_b14_deg3_bessel_pinch)
        .collect();
    let multi_octave_discharging_certificates: Vec<MultiOctaveDischargingCertificate> = hash_graphs
        .par_iter()
        .map(certify_multi_octave_discharging_closure)
        .collect();

    // Run Hard-(C_4=C_8=0) C_16-Minimizing Tabu Search on the C_4=C_8=0 frontier seeds
    let mut tabu_seeds = Vec::new();
    for g in &p3_lib {
        if g.expected_girth == 6 {
            tabu_seeds.push(g.clone());
        }
    }
    tabu_seeds.push(build_zm_dipole_voltage_lift(
        "Dipole Lift Z_31(0,1,6) [n=62]",
        31,
        [0, 1, 6],
        6,
    ));

    let tabu_c16_records: Vec<C16TabuSearchRecord> = tabu_seeds
        .iter()
        .map(|g| run_hard_c4_c8_zero_c16_tabu_search(g, restarts, steps_per_restart, 20261026))
        .collect();

    Phase3CertificationBundle {
        free_abelian_certificate,
        commutator_frontier_records,
        hashimoto_trace_certificates,
        spectral_doubling_certificates,
        topological_n16_certificates,
        b10_hankel_certificates,
        b12_deg2_bessel_certificates,
        b14_deg3_bessel_certificates,
        multi_octave_discharging_certificates,
        tabu_c16_records,
    }
}

/// Prints a formatted ASCII report of Theorem 2.3 and Theorem 3.2 certificates.
pub fn print_phase3_report(bundle: &Phase3CertificationBundle) {
    let fa = &bundle.free_abelian_certificate;
    println!("========================================================================================================================");
    println!("  Theorem 2.3(a): Universal Abelian Commutator Floor in Free Abelian Cover Lift_{{Z^2}}(D_3; (0,0),(1,0),(0,1))");
    println!("========================================================================================================================");
    println!(
        "  Directed simple 4-cycles from u_0 : {:>6}   | Universal C_4 Gamma-orbits        : {:>6}",
        fa.directed_simple_c4,
        fa.directed_simple_c4 / 8
    );
    println!(
        "  Directed simple 6-cycles from u_0 : {:>6}   | Universal C_6 Gamma-orbits        : {:>6} (forces g <= 6 in all abelian lifts!)",
        fa.directed_simple_c6, fa.universal_c6_orbits
    );
    println!(
        "  Directed simple 8-cycles from u_0 : {:>6}   | Universal C_8 Gamma-orbits        : {:>6}",
        fa.directed_simple_c8,
        fa.directed_simple_c8 / 16
    );
    println!(
        "  Closed non-backtracking 16-walks  : {:>6}   | Universal NB-16 Gamma-orbits      : {:>6}",
        fa.closed_nonbacktracking_16walks, fa.universal_nb16_orbits
    );
    println!(
        "  Directed simple 16-cycles from u_0: {:>6}   | Universal simple C_16 Orbits      : {:>6} (=> C_16 >= 18 |Gamma| universal floor!)",
        fa.directed_simple_c16, fa.universal_c16_orbits
    );
    println!(
        "  Universal non-simple NB-16 Orbits : {:>6}   | Nonzero 16-walk Z^2 displacements : {:>6}",
        fa.universal_nonsimple_nb16_orbits, fa.distinct_nonzero_16walk_displacements
    );
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorem 2.3(b): Optimal-Abelian vs. Non-Abelian / Periodic-LCF Commutator Frontier (50 <= n <= 64)");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>4} {:>5} {:>5} {:>6} {:>6} {:>8} {:>9} {:>10}",
        "Graph", "n", "m", "Girth", "Bip?", "C_4", "C_8", "C_16", "C_16/|G|", "C_32^(0)"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for r in &bundle.commutator_frontier_records {
        let bip = if r.is_bipartite { "yes" } else { "no" };
        let fiber = (r.n / 2) as f64;
        println!(
            "{:<36} {:>4} {:>4} {:>5} {:>5} {:>6} {:>6} {:>8} {:>9.2} {:>10}",
            r.name,
            r.n,
            r.m,
            r.girth,
            bip,
            r.c4,
            r.c8,
            r.c16,
            (r.c16 as f64) / fiber,
            r.c32
        );
    }
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorem 3.2(a): Bass-Ihara Chebyshev Hashimoto Trace Identity Tr(B^{{2^k}}) vs. Simple Dyadic Cycles");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>5} {:>8} {:>8} {:>10} {:>12} {:>12} {:>10}",
        "Graph", "n", "Girth", "Tr(B^4)", "Tr(B^8)", "C_8", "Tr(B^16)/32", "C_16", "NB16 Excess"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for h in &bundle.hashimoto_trace_certificates {
        println!(
            "{:<36} {:>4} {:>5} {:>8} {:>8} {:>10} {:>12} {:>12} {:>10}",
            h.name,
            h.n,
            h.girth,
            h.tr_b4_chebyshev,
            h.tr_b8_chebyshev,
            h.c8,
            h.tr_b16_chebyshev / 32,
            h.c16,
            h.nonsimple_nb16_walk_orbits
        );
    }
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorem 2.4: Chebyshev-Hashimoto Doubling Variance Floor & Golden-Ratio Hankel Moment Bounds on C_6(G)");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>2} {:>5} {:>10} {:>10} {:>8} {:>8} {:>8} {:>14}",
        "Graph",
        "n",
        "r",
        "Girth",
        "Tr(B^16)",
        "Floor B^16",
        "C_6",
        "C_6_min",
        "C_6_max",
        "det H_3(mu)"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for s in &bundle.spectral_doubling_certificates {
        println!(
            "{:<36} {:>4} {:>2} {:>5} {:>10} {:>10} {:>8} {:>8} {:>8} {:>14}",
            s.name,
            s.n,
            s.r_trivial,
            s.girth,
            s.tr_b16_actual,
            s.tr_b16_spectral_floor,
            s.c6,
            s.c6_hankel_lower_bound,
            s.c6_hankel_upper_bound,
            s.hankel_det_h3
        );
    }
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorem 2.5: Exact Topological Decomposition of Tr(B^16) & Mandatory 6-Cycle Participation in N_16(G)");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>4} {:>8} {:>7} {:>8} {:>8} {:>7} {:>7} {:>8} {:>8}",
        "Graph",
        "n",
        "C_6",
        "32 C_16",
        "16 C_8",
        "Db(6,6)",
        "Theta",
        "Th155",
        "Th248",
        "Th159",
        "N16/C6"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for top in &bundle.topological_n16_certificates {
        let ratio = if top.c6 > 0 {
            format!("{:.1}", (top.n16_total as f64) / (top.c6 as f64))
        } else {
            "0.0".to_string()
        };
        println!(
            "{:<36} {:>4} {:>4} {:>8} {:>7} {:>8} {:>8} {:>7} {:>7} {:>8} {:>8}",
            top.name,
            top.n,
            top.c6,
            top.walks_c16,
            top.walks_c8,
            top.walks_dumbbell_total,
            top.walks_theta_total,
            top.walks_theta_1_5_5,
            top.walks_theta_2_4_8,
            top.walks_theta_1_5_9,
            ratio
        );
    }
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorem 2.6: Exact Order-10 Hashimoto Trace Tr(B^10)=20 C_10, 3x3 Hankel Pinching & Bessel Quadratic C_16 Floor");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>4} {:>6} {:>8} {:>7} {:>7} {:>9} {:>9} {:>8} {:>8}",
        "Graph",
        "n",
        "C_6",
        "C_10",
        "Tr(B^10)",
        "C10_min",
        "C10_max",
        "320*mu*C10",
        "BesselB16",
        "C16_lb",
        "C16_act"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for b10 in &bundle.b10_hankel_certificates {
        println!(
            "{:<36} {:>4} {:>4} {:>6} {:>8} {:>7} {:>7} {:>9} {:>9} {:>8} {:>8}",
            b10.name,
            b10.n,
            b10.c6,
            b10.c10,
            b10.tr_b10,
            b10.c10_hankel_lower_bound,
            b10.c10_hankel_upper_bound,
            b10.n16_theta_c10_bound,
            b10.tr_b16_bessel_boosted_floor,
            b10.coupled_c16_lower_bound,
            b10.c16
        );
    }
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorem 2.7: Exact Order-12 Hashimoto Trace Tr(B^12), 3D Bessel Projection on span{{1,lambda^2,lambda^4}} & C_16 Floor");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>4} {:>5} {:>5} {:>5} {:>8} {:>9} {:>9} {:>8} {:>8}",
        "Graph",
        "n",
        "C_6",
        "C_12",
        "tau1",
        "tau3",
        "Tr(B^12)",
        "2DBessel",
        "3DBessel",
        "C16_3Dlb",
        "C16_act"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for b12 in &bundle.b12_deg2_bessel_certificates {
        println!(
            "{:<36} {:>4} {:>4} {:>5} {:>5} {:>5} {:>8} {:>9} {:>9} {:>8} {:>8}",
            b12.name,
            b12.n,
            b12.c6,
            b12.c12,
            b12.tau1,
            b12.tau3,
            b12.tr_b12,
            b12.tr_b16_2d_bessel_floor,
            b12.tr_b16_3d_bessel_floor,
            b12.coupled_c16_3d_lower_bound,
            b12.c16
        );
    }
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorem 2.8: Exact Tr(B^14) == 28 C_14 + 56 beta_1 & 4D Bessel Projection on span{{1, lambda^2, lambda^4, lambda^6}}");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>4} {:>6} {:>6} {:>8} {:>9} {:>9} {:>8} {:>8}",
        "Graph",
        "n",
        "C_6",
        "C_14",
        "beta1",
        "Tr(B^14)",
        "3DBessel",
        "4DBessel",
        "C16_4Dlb",
        "C16_act"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for b14 in &bundle.b14_deg3_bessel_certificates {
        println!(
            "{:<36} {:>4} {:>4} {:>6} {:>6} {:>8} {:>9} {:>9} {:>8} {:>8}",
            b14.name,
            b14.n,
            b14.c6,
            b14.c14,
            b14.beta1,
            b14.tr_b14,
            b14.tr_b16_3d_bessel_floor,
            b14.tr_b16_4d_bessel_floor,
            b14.coupled_c16_4d_lower_bound,
            b14.c16
        );
    }
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorems 2.9 & 2.10 + Corollary 2.11: Multi-Octave Chebyshev-Hashimoto Spectral-Topological Discharging Closure");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>5} {:>5} {:>4} {:>5} {:>9} {:>9} {:>9} {:>8}",
        "Graph",
        "n",
        "Theta",
        "Prism",
        "kDom",
        "M_dom",
        "PerronMrg",
        "4DBessel",
        "N16_ub",
        "Closure"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for p9 in &bundle.multi_octave_discharging_certificates {
        println!(
            "{:<36} {:>4} {:>5} {:>5} {:>4} {:>5} {:>9} {:>9} {:>9} {:>8}",
            p9.name,
            p9.n,
            if p9.theta_conservation_verified {
                "OK"
            } else {
                "FAIL"
            },
            if p9.prism_conservation_verified {
                "OK"
            } else {
                "FAIL"
            },
            p9.perron_dominance_k,
            p9.perron_dominance_m,
            p9.perron_dominance_margin,
            p9.tr_b16_4d_bessel_floor,
            p9.n16_discharged_upper_bound,
            if p9.dyadic_closure_certified {
                "CERT"
            } else {
                "FAIL"
            }
        );
    }
    println!("========================================================================================================================\n");

    println!("========================================================================================================================");
    println!("  Theorem 3.2(b): Hard-(C_4=C_8=0) Barrier + C_16-Minimizing Tabu Search Breaking the Abelian Floor");
    println!("========================================================================================================================");
    println!(
        "{:<36} {:>4} {:>10} {:>12} {:>12} {:>8} {:>6} {:>6} {:>8}",
        "Seed Graph (C_4=C_8=0)",
        "n",
        "Init C_16",
        "Min C_16(g>=6)",
        "Min C_16(g>=5)",
        "End g",
        "C_4",
        "C_8",
        "Bip?"
    );
    println!("------------------------------------------------------------------------------------------------------------------------");
    for t in &bundle.tabu_c16_records {
        println!(
            "{:<36} {:>4} {:>10} {:>12} {:>12} {:>8} {:>6} {:>6} {:>8}",
            t.seed_name,
            t.n,
            t.initial_c16,
            t.best_c16_girth6,
            t.best_c16_girth5,
            t.best_girth,
            t.best_c4,
            t.best_c8,
            if t.best_is_bipartite { "yes" } else { "no" }
        );
    }
    println!("========================================================================================================================");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_certify_phase3_commutator_and_hashimoto() {
        let bundle = certify_phase3_bundle(2, 40);
        assert_eq!(bundle.free_abelian_certificate.universal_c16_orbits, 18);
        assert_eq!(bundle.free_abelian_certificate.universal_nb16_orbits, 60);
        assert_eq!(bundle.commutator_frontier_records.len(), 10);
        for r in &bundle.commutator_frontier_records {
            assert_eq!(r.min_deg, 3);
            assert_eq!(r.max_deg, 3);
            assert_eq!(r.c4, 0);
            assert!(r.total_dyadic_cycles > 0);
        }
        // Check that Non-Abelian Per_3(Z_27;3,3,9) at n=54 strictly beats Abelian Opt Z_27(0,1,6)
        assert_eq!(bundle.commutator_frontier_records[2].c16, 2484);
        assert_eq!(bundle.commutator_frontier_records[3].c16, 2385);
        // Check that Non-Abelian LCF at n=60 achieves 2340 < 2550 abelian floor
        assert_eq!(bundle.commutator_frontier_records[6].c16, 2340);
        // Check that Abelian Opt Z_31(0,1,12) at n=62 achieves 2294 = 74 * 31 < 2883
        assert_eq!(bundle.commutator_frontier_records[7].c16, 2294);
        // Check that Abelian Opt Z_32(0,1,7) at n=64 achieves 2232
        assert_eq!(bundle.commutator_frontier_records[8].c16, 2232);
        // Check that Metacyclic Z_8 x_3 Z_4 at n=64 breaks the abelian g<=6 barrier to reach girth 8
        assert_eq!(bundle.commutator_frontier_records[9].girth, 8);
        assert_eq!(bundle.spectral_doubling_certificates.len(), 16);
        for s in &bundle.spectral_doubling_certificates {
            assert!(s.tr_b8_actual >= s.tr_b8_spectral_floor);
            assert!(s.tr_b16_actual >= s.tr_b16_spectral_floor);
        }
        assert_eq!(bundle.topological_n16_certificates.len(), 16);
        for top in &bundle.topological_n16_certificates {
            assert_eq!(
                top.tr_b16_chebyshev,
                top.walks_c16
                    + top.walks_c8
                    + top.walks_c4
                    + top.walks_dumbbell_total
                    + top.walks_theta_total
            );
        }
        assert_eq!(bundle.b10_hankel_certificates.len(), 16);
        for b10 in &bundle.b10_hankel_certificates {
            if b10.girth >= 6 {
                assert_eq!(b10.tr_b10, 20 * b10.c10);
                assert_eq!(b10.n10_excess, 0);
            }
            if b10.is_bipartite && b10.c4 == 0 && b10.c8 == 0 {
                assert!(b10.c16 >= b10.coupled_c16_lower_bound);
                assert!(b10.coupled_c16_lower_bound >= 581);
            }
        }
        assert_eq!(bundle.b12_deg2_bessel_certificates.len(), 16);
        for b12 in &bundle.b12_deg2_bessel_certificates {
            if b12.girth >= 6 && b12.c8 == 0 {
                assert_eq!(
                    b12.tr_b12,
                    24 * b12.c12 + 12 * b12.c6 + 24 * b12.tau1 + 24 * b12.tau3
                );
                assert_eq!(b12.n12_composite_actual, b12.n12_composite_predicted);
            }
            if b12.is_bipartite && b12.c4 == 0 && b12.c8 == 0 {
                assert!(b12.tr_b16_3d_bessel_floor >= b12.tr_b16_2d_bessel_floor);
                assert!(b12.c16 >= b12.coupled_c16_3d_lower_bound);
                assert!(b12.coupled_c16_3d_lower_bound >= b12.coupled_c16_2d_lower_bound);
            }
        }
        assert_eq!(bundle.b14_deg3_bessel_certificates.len(), 16);
        for b14 in &bundle.b14_deg3_bessel_certificates {
            if b14.is_bipartite && b14.c4 == 0 && b14.c8 == 0 {
                assert_eq!(b14.tr_b14, 28 * b14.c14 + 56 * b14.beta1);
                assert_eq!(b14.n14_composite_actual, b14.n14_composite_predicted);
                assert!(b14.tr_b16_4d_bessel_floor >= b14.tr_b16_3d_bessel_floor);
                assert!(b14.c16 >= b14.coupled_c16_4d_lower_bound);
                assert!(b14.coupled_c16_4d_lower_bound >= b14.coupled_c16_3d_lower_bound);
            }
        }
        assert_eq!(bundle.multi_octave_discharging_certificates.len(), 16);
        for p9 in &bundle.multi_octave_discharging_certificates {
            assert!(p9.theta_conservation_verified);
            assert!(p9.prism_conservation_verified);
            assert!(p9.menger_2_edge_connected);
            assert!(p9.dyadic_closure_certified);
            assert!(p9.total_dyadic_cycles_up_to_16 > 0);
        }
    }

    #[test]
    fn test_certify_canonical_library_up_to_20() {
        let records = certify_canonical_library(20);
        assert_eq!(records.len(), 18);
        for r in &records {
            assert_eq!(r.min_deg, 3);
            assert_eq!(r.max_deg, 3);
            assert!(
                r.total_dyadic_cycles > 0,
                "Every canonical graph must have total_dyadic_cycles > 0 (failed on {})",
                r.name
            );
            assert!(
                !r.dyadic_lengths_present.is_empty(),
                "dyadic_lengths_present must be non-empty on {}",
                r.name
            );
        }
    }

    #[test]
    fn test_certify_open_frontier_50_to_62() {
        let records = certify_open_frontier_library();
        assert_eq!(records.len(), 7);
        let expected_ns = [50, 52, 54, 56, 58, 60, 62];
        for (r, &exp_n) in records.iter().zip(expected_ns.iter()) {
            assert_eq!(r.n, exp_n);
            assert_eq!(r.min_deg, 3);
            assert_eq!(r.max_deg, 3);
            assert!(
                r.girth >= 6,
                "Frontier graph {} must have girth >= 6 (got {})",
                r.name,
                r.girth
            );
            assert_eq!(
                r.c4, 0,
                "Frontier graph {} has girth >= 6 so C_4 = 0",
                r.name
            );
            assert!(
                r.total_dyadic_cycles > 0,
                "Frontier graph {} must contain dyadic cycles",
                r.name
            );
        }
    }
}
