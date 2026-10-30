//! ============================================================================
//! `erdos64` — Hybrid CPU+CUDA Certification & Counterexample Hunter for
//! Erdős Problem #64 ($1,000 Falsifiable Prize — Erdős–Gyárfás 2^k-Cycle Conjecture).
//! ============================================================================

pub mod exact;
pub mod gpu;
pub mod records;
pub mod search;
pub mod verify_cpu;

use anyhow::Result;
use clap::{Parser, Subcommand};

use crate::exact::{
    certify_canonical_library, certify_open_frontier_library, certify_phase3_bundle,
    print_cage_certification_report, print_open_frontier_report, print_phase3_report,
};
use crate::gpu::Gpu;
use crate::records::{
    build_canonical_cage_and_snark_library, build_open_frontier_voltage_and_lcf_library,
    build_phase3_commutator_frontier_library,
};
use crate::search::run_dyadic_annealing_search;
use crate::verify_cpu::{
    eval_b10_hankel_pinch_cpu, eval_b12_deg2_bessel_pinch_cpu, eval_b14_deg3_bessel_pinch_cpu,
    eval_dyadic_cycles_cpu, eval_hashimoto_b4_b8_cpu, eval_multi_octave_discharging_cpu,
    eval_n16_bicyclic_walks_cpu,
};

#[derive(Parser, Debug)]
#[command(
    name = "erdos64",
    version,
    about = "GPU-accelerated certification and counterexample hunter for Erdős Problem #64 ($1,000 Falsifiable Prize: the Erdős–Gyárfás 2^k-Cycle Conjecture)"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Certifies Theorem 2.1 & 3.1 across the 18 canonical (3,g)-cages, generalized polygons, and snarks.
    CertifyCages {
        /// Maximum vertex count `n` for full cycle-spectrum enumeration (`C_l(G)` for all `3 <= l <= n`).
        #[arg(long, default_value_t = 24)]
        max_full_n: usize,
        /// Optional path to export JSON certification records.
        #[arg(long)]
        json_out: Option<String>,
    },
    /// Certifies Theorem 2.2 on the Open Frontier (`50 <= n <= 62`) cyclic voltage-lift and LCF-code cubic graphs.
    CertifyFrontier {
        /// Optional path to export JSON certification records.
        #[arg(long)]
        json_out: Option<String>,
    },
    /// Certifies Phase 3 (Theorem 2.3 Universal Abelian Commutator Floor + Non-Abelian / Periodic-LCF Breaking
    /// and Theorem 3.2 Bass-Ihara Hashimoto Spectral Trace + Hard-(C_4=C_8=0) C_16-Minimizing Tabu Search).
    CertifyPhase3 {
        /// Number of parallel Tabu restarts per seed graph.
        #[arg(long, default_value_t = 8)]
        restarts: usize,
        /// Number of 2-edge-switch Tabu steps per restart.
        #[arg(long, default_value_t = 250)]
        steps: usize,
        /// Optional path to export JSON Phase 3 bundle.
        #[arg(long)]
        json_out: Option<String>,
    },
    /// Runs 2-edge-switch simulated annealing on cubic graphs to minimize dyadic cycle energy.
    SearchCubic {
        /// Minimum even vertex count `n`.
        #[arg(long, default_value_t = 10)]
        min_n: usize,
        /// Maximum even vertex count `n`.
        #[arg(long, default_value_t = 24)]
        max_n: usize,
        /// Number of parallel Markov chain restarts per `n`.
        #[arg(long, default_value_t = 16)]
        restarts: usize,
        /// Number of 2-edge-switch annealing steps per restart.
        #[arg(long, default_value_t = 1500)]
        steps: usize,
        /// Base RNG seed.
        #[arg(long, default_value_t = 20260417)]
        seed: u64,
    },
    /// Verifies bit-exact CPU vs CUDA GPU parity on `k_eval_dyadic_cycles`.
    GpuParity {
        /// CUDA device ordinal.
        #[arg(long, default_value_t = 0)]
        device: usize,
        /// NVRTC compute architecture target (e.g., `compute_86`).
        #[arg(long, default_value = "compute_86")]
        arch: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::CertifyCages {
            max_full_n,
            json_out,
        } => {
            let records = certify_canonical_library(max_full_n);
            print_cage_certification_report(&records);
            if let Some(path) = json_out {
                let json = serde_json::to_string_pretty(&records)?;
                std::fs::write(&path, json)?;
                println!("Saved JSON certification report to {path}");
            }
        }
        Commands::CertifyFrontier { json_out } => {
            let records = certify_open_frontier_library();
            print_open_frontier_report(&records);
            if let Some(path) = json_out {
                let json = serde_json::to_string_pretty(&records)?;
                std::fs::write(&path, json)?;
                println!("Saved JSON frontier report to {path}");
            }
        }
        Commands::CertifyPhase3 {
            restarts,
            steps,
            json_out,
        } => {
            let bundle = certify_phase3_bundle(restarts, steps);
            print_phase3_report(&bundle);
            if let Some(path) = json_out {
                let json = serde_json::to_string_pretty(&bundle)?;
                std::fs::write(&path, json)?;
                println!("Saved JSON Phase 3 report to {path}");
            }
        }
        Commands::SearchCubic {
            min_n,
            max_n,
            restarts,
            steps,
            seed,
        } => {
            println!("====================================================================================");
            println!("  Theorem 3.1: 2-Edge-Switch Simulated Annealing on 3-Regular Graphs (min_deg = 3)");
            println!("====================================================================================");
            println!(
                "{:<6} {:>10} {:>10} {:>8} {:>8} {:>10} {:>10} {:>12} {:>8}",
                "n", "Restarts", "Steps", "Girth", "C_4", "Min C_8", "C_16", "Min 2^k", "Bip?"
            );
            println!("------------------------------------------------------------------------------------");
            let mut n = if min_n.is_multiple_of(2) {
                min_n
            } else {
                min_n + 1
            };
            while n <= max_n {
                let s = run_dyadic_annealing_search(n, restarts, steps, seed);
                println!(
                    "{:<6} {:>10} {:>10} {:>8} {:>8} {:>10} {:>10} {:>12} {:>8}",
                    s.n,
                    s.restarts,
                    s.steps_per_restart,
                    s.best_girth,
                    s.best_c4,
                    s.best_c8,
                    s.best_c16,
                    s.best_dyadic_total,
                    if s.is_bipartite { "yes" } else { "no" }
                );
                n += 2;
            }
            println!("====================================================================================");
        }
        Commands::GpuParity { device, arch } => {
            let gpu = Gpu::init(device, &arch)?;
            println!(
                "Initialized CUDA GPU: {} ({} SMs, arch={})",
                gpu.device_name, gpu.sm_count, arch
            );
            let lib = build_canonical_cage_and_snark_library();
            for g in &lib {
                if g.n > 24 {
                    continue;
                }
                let padded = g.padded_adj_64();
                let cpu_res = eval_dyadic_cycles_cpu(&padded, g.n, 1)[0];
                let gpu_res = gpu.eval_dyadic_cycles_batch(&padded, g.n, 1)?[0];
                assert_eq!(cpu_res, gpu_res, "Parity mismatch on {}", g.name);
                let cpu_hash = eval_hashimoto_b4_b8_cpu(&padded, g.n, 1)[0];
                let gpu_hash = gpu.eval_hashimoto_b4_b8_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_hash, gpu_hash,
                    "Hashimoto parity mismatch on {}",
                    g.name
                );
                let cpu_n16 = eval_n16_bicyclic_walks_cpu(&padded, g.n, 1)[0];
                let gpu_n16 = gpu.eval_n16_bicyclic_walks_batch(&padded, g.n, 1)?[0];
                assert_eq!(cpu_n16, gpu_n16, "N_16 walk parity mismatch on {}", g.name);
                let cpu_b10 = eval_b10_hankel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b10 = gpu.eval_b10_hankel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b10, gpu_b10,
                    "B^10 Hankel parity mismatch on {}",
                    g.name
                );
                let cpu_b12 = eval_b12_deg2_bessel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b12 = gpu.eval_b12_deg2_bessel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b12, gpu_b12,
                    "B^12 Bessel parity mismatch on {}",
                    g.name
                );
                let cpu_b14 = eval_b14_deg3_bessel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b14 = gpu.eval_b14_deg3_bessel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b14, gpu_b14,
                    "B^14 Bessel parity mismatch on {}",
                    g.name
                );
                let cpu_p9 = eval_multi_octave_discharging_cpu(&padded, g.n, 1)[0];
                let gpu_p9 = gpu.eval_multi_octave_discharging_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_p9, gpu_p9,
                    "Phase 9 multi-octave discharging parity mismatch on {}",
                    g.name
                );
                println!(
                    "  [OK] {:<26} (n={:>2}): C_4={:>3}, C_8={:>4}, C_14={:>4}, Tr(B^14)={:>5}, 4DBessel={:>6}, P9Closure={}",
                    g.name, g.n, gpu_res.c4, gpu_res.c8, gpu_b14.c14, gpu_b14.tr_b14, gpu_b14.tr_b16_4d_bessel_floor, gpu_p9.dyadic_closure_certified
                );
            }
            let frontier = build_open_frontier_voltage_and_lcf_library();
            for g in &frontier {
                let padded = g.padded_adj_64();
                let cpu_res = eval_dyadic_cycles_cpu(&padded, g.n, 1)[0];
                let gpu_res = gpu.eval_dyadic_cycles_batch(&padded, g.n, 1)?[0];
                assert_eq!(cpu_res, gpu_res, "Parity mismatch on frontier {}", g.name);
                let cpu_hash = eval_hashimoto_b4_b8_cpu(&padded, g.n, 1)[0];
                let gpu_hash = gpu.eval_hashimoto_b4_b8_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_hash, gpu_hash,
                    "Hashimoto parity mismatch on frontier {}",
                    g.name
                );
                let cpu_n16 = eval_n16_bicyclic_walks_cpu(&padded, g.n, 1)[0];
                let gpu_n16 = gpu.eval_n16_bicyclic_walks_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_n16, gpu_n16,
                    "N_16 walk parity mismatch on frontier {}",
                    g.name
                );
                let cpu_b10 = eval_b10_hankel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b10 = gpu.eval_b10_hankel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b10, gpu_b10,
                    "B^10 Hankel parity mismatch on frontier {}",
                    g.name
                );
                let cpu_b12 = eval_b12_deg2_bessel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b12 = gpu.eval_b12_deg2_bessel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b12, gpu_b12,
                    "B^12 Bessel parity mismatch on frontier {}",
                    g.name
                );
                let cpu_b14 = eval_b14_deg3_bessel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b14 = gpu.eval_b14_deg3_bessel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b14, gpu_b14,
                    "B^14 Bessel parity mismatch on frontier {}",
                    g.name
                );
                let cpu_p9 = eval_multi_octave_discharging_cpu(&padded, g.n, 1)[0];
                let gpu_p9 = gpu.eval_multi_octave_discharging_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_p9, gpu_p9,
                    "Phase 9 multi-octave discharging parity mismatch on frontier {}",
                    g.name
                );
                println!(
                    "  [OK] {:<30} (n={:>2}): C_4={:>3}, C_8={:>4}, C_14={:>4}, Tr(B^14)={:>5}, 4DBessel={:>6}, P9Closure={}",
                    g.name, g.n, gpu_res.c4, gpu_res.c8, gpu_b14.c14, gpu_b14.tr_b14, gpu_b14.tr_b16_4d_bessel_floor, gpu_p9.dyadic_closure_certified
                );
            }
            let phase3 = build_phase3_commutator_frontier_library();
            for g in &phase3 {
                let padded = g.padded_adj_64();
                let cpu_hash = eval_hashimoto_b4_b8_cpu(&padded, g.n, 1)[0];
                let gpu_hash = gpu.eval_hashimoto_b4_b8_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_hash, gpu_hash,
                    "Hashimoto parity mismatch on phase3 {}",
                    g.name
                );
                let cpu_n16 = eval_n16_bicyclic_walks_cpu(&padded, g.n, 1)[0];
                let gpu_n16 = gpu.eval_n16_bicyclic_walks_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_n16, gpu_n16,
                    "N_16 walk parity mismatch on phase3 {}",
                    g.name
                );
                let cpu_b10 = eval_b10_hankel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b10 = gpu.eval_b10_hankel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b10, gpu_b10,
                    "B^10 Hankel parity mismatch on phase3 {}",
                    g.name
                );
                let cpu_b12 = eval_b12_deg2_bessel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b12 = gpu.eval_b12_deg2_bessel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b12, gpu_b12,
                    "B^12 Bessel parity mismatch on phase3 {}",
                    g.name
                );
                let cpu_b14 = eval_b14_deg3_bessel_pinch_cpu(&padded, g.n, 1)[0];
                let gpu_b14 = gpu.eval_b14_deg3_bessel_pinch_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_b14, gpu_b14,
                    "B^14 Bessel parity mismatch on phase3 {}",
                    g.name
                );
                let cpu_p9 = eval_multi_octave_discharging_cpu(&padded, g.n, 1)[0];
                let gpu_p9 = gpu.eval_multi_octave_discharging_batch(&padded, g.n, 1)?[0];
                assert_eq!(
                    cpu_p9, gpu_p9,
                    "Phase 9 multi-octave discharging parity mismatch on phase3 {}",
                    g.name
                );
                println!(
                    "  [OK] {:<36} (n={:>2}): C_6={:>2}, C_14={:>4}, Tr(B^14)={:>5}, 4DBessel={:>6}, N16_ub={:>5}, P9Closure={}",
                    g.name,
                    g.n,
                    gpu_b14.c6,
                    gpu_b14.c14,
                    gpu_b14.tr_b14,
                    gpu_b14.tr_b16_4d_bessel_floor,
                    gpu_p9.n16_discharged_upper_bound,
                    gpu_p9.dyadic_closure_certified
                );
            }
            println!("All CPU vs CUDA GPU kernel parity checks passed bit-for-bit!");
        }
    }
    Ok(())
}
