//! ============================================================================
//! `gpu.rs` — CUDA driver context and runtime NVRTC compilation for `erdos64`.
//! ============================================================================
//!
//! Uses `cudarc` with `dynamic-loading`, so `libcuda.so` and `libnvrtc.so` are
//! loaded via `dlopen` only when a GPU path is invoked. CPU-only commands and
//! the default CI test suite run on hosts without an Nvidia GPU.

use anyhow::{Context, Result};
use cudarc::driver::{
    CudaContext, CudaFunction, CudaModule, CudaStream, LaunchConfig, PushKernelArg,
};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};
use std::sync::Arc;

use crate::records::MAX_VERTICES;
use crate::verify_cpu::{
    B10HankelPinchSummary, B12Deg2BesselSummary, B14Deg3BesselSummary, DyadicEvalSummary,
    HashimotoB4B8Summary, MultiOctaveDischargingSummary, N16WalkDecompositionSummary,
};

pub const KERNEL_SRC: &str = include_str!("kernels.cu");

pub struct Gpu {
    pub ctx: Arc<CudaContext>,
    pub stream: Arc<CudaStream>,
    pub module: Arc<CudaModule>,
    pub fn_eval_dyadic_cycles: CudaFunction,
    pub fn_eval_hashimoto_trace_b4_b8: CudaFunction,
    pub fn_enumerate_n16_bicyclic_walks: CudaFunction,
    pub fn_eval_b10_hankel_pinch: CudaFunction,
    pub fn_eval_b12_deg2_bessel_pinch: CudaFunction,
    pub fn_eval_b14_deg3_bessel_pinch: CudaFunction,
    pub fn_eval_multi_octave_discharging: CudaFunction,
    pub device_name: String,
    pub sm_count: i32,
}

impl Gpu {
    pub fn init(device_ordinal: usize, arch: &str) -> Result<Self> {
        let ctx = CudaContext::new(device_ordinal)
            .with_context(|| format!("failed to initialize CUDA device {device_ordinal}"))?;
        let stream = ctx.default_stream();
        let device_name = ctx
            .name()
            .unwrap_or_else(|_| "Unknown CUDA GPU".to_string());
        let sm_count = ctx
            .attribute(
                cudarc::driver::sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT,
            )
            .unwrap_or(82);

        let arch_static: &'static str = Box::leak(arch.to_owned().into_boxed_str());
        let opts = CompileOptions {
            arch: Some(arch_static),
            options: vec!["--device-int128".to_string()],
            ..Default::default()
        };
        let ptx = compile_ptx_with_opts(KERNEL_SRC, opts)
            .context("NVRTC compilation of src/kernels.cu failed")?;
        let module = ctx
            .load_module(ptx)
            .context("failed to load compiled PTX module")?;

        let fn_eval_dyadic_cycles = module
            .load_function("k_eval_dyadic_cycles")
            .context("missing kernel k_eval_dyadic_cycles")?;
        let fn_eval_hashimoto_trace_b4_b8 = module
            .load_function("k_eval_hashimoto_trace_b4_b8")
            .context("missing kernel k_eval_hashimoto_trace_b4_b8")?;
        let fn_enumerate_n16_bicyclic_walks = module
            .load_function("k_enumerate_n16_bicyclic_walks")
            .context("missing kernel k_enumerate_n16_bicyclic_walks")?;
        let fn_eval_b10_hankel_pinch = module
            .load_function("k_eval_b10_hankel_pinch")
            .context("missing kernel k_eval_b10_hankel_pinch")?;
        let fn_eval_b12_deg2_bessel_pinch = module
            .load_function("k_eval_b12_deg2_bessel_pinch")
            .context("missing kernel k_eval_b12_deg2_bessel_pinch")?;
        let fn_eval_b14_deg3_bessel_pinch = module
            .load_function("k_eval_b14_deg3_bessel_pinch")
            .context("missing kernel k_eval_b14_deg3_bessel_pinch")?;
        let fn_eval_multi_octave_discharging = module
            .load_function("k_eval_multi_octave_discharging")
            .context("missing kernel k_eval_multi_octave_discharging")?;

        Ok(Self {
            ctx,
            stream,
            module,
            fn_eval_dyadic_cycles,
            fn_eval_hashimoto_trace_b4_b8,
            fn_enumerate_n16_bicyclic_walks,
            fn_eval_b10_hankel_pinch,
            fn_eval_b12_deg2_bessel_pinch,
            fn_eval_b14_deg3_bessel_pinch,
            fn_eval_multi_octave_discharging,
            device_name,
            sm_count,
        })
    }

    pub fn launch_cfg(&self, n_items: u32) -> LaunchConfig {
        let block_dim = 128u32;
        let grid_dim = n_items.div_ceil(block_dim).max(1);
        LaunchConfig {
            grid_dim: (grid_dim, 1, 1),
            block_dim: (block_dim, 1, 1),
            shared_mem_bytes: 0,
        }
    }

    /// Evaluates `k_eval_dyadic_cycles` on a batch of `batch_size` graphs (`flat_adj.len() == batch_size * 64`).
    pub fn eval_dyadic_cycles_batch(
        &self,
        flat_adj: &[u64],
        n: usize,
        batch_size: usize,
    ) -> Result<Vec<DyadicEvalSummary>> {
        assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
        let d_adj = self.stream.clone_htod(flat_adj)?;
        let mut d_c4 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c8 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c16_0 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_min_deg = self.stream.alloc_zeros::<u32>(batch_size)?;

        let n_u32 = n as u32;
        let b_u32 = batch_size as u32;
        let cfg = self.launch_cfg(b_u32);

        unsafe {
            self.stream
                .launch_builder(&self.fn_eval_dyadic_cycles)
                .arg(&d_adj)
                .arg(&n_u32)
                .arg(&b_u32)
                .arg(&mut d_c4)
                .arg(&mut d_c8)
                .arg(&mut d_c16_0)
                .arg(&mut d_min_deg)
                .launch(cfg)?;
        }

        let h_c4 = self.stream.clone_dtoh(&d_c4)?;
        let h_c8 = self.stream.clone_dtoh(&d_c8)?;
        let h_c16_0 = self.stream.clone_dtoh(&d_c16_0)?;
        let h_min_deg = self.stream.clone_dtoh(&d_min_deg)?;

        let mut out = Vec::with_capacity(batch_size);
        for i in 0..batch_size {
            out.push(DyadicEvalSummary {
                c4: h_c4[i],
                c8: h_c8[i],
                c16_anchored0: h_c16_0[i],
                min_deg: h_min_deg[i],
            });
        }
        Ok(out)
    }

    /// Evaluates `k_eval_hashimoto_trace_b4_b8` (Theorem 3.2a) on a batch of `batch_size` graphs.
    pub fn eval_hashimoto_b4_b8_batch(
        &self,
        flat_adj: &[u64],
        n: usize,
        batch_size: usize,
    ) -> Result<Vec<HashimotoB4B8Summary>> {
        assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
        let d_adj = self.stream.clone_htod(flat_adj)?;
        let mut d_tr_a4 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_a6 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_a8 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b4 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b8 = self.stream.alloc_zeros::<u64>(batch_size)?;

        let n_u32 = n as u32;
        let b_u32 = batch_size as u32;
        let cfg = self.launch_cfg(b_u32);

        unsafe {
            self.stream
                .launch_builder(&self.fn_eval_hashimoto_trace_b4_b8)
                .arg(&d_adj)
                .arg(&n_u32)
                .arg(&b_u32)
                .arg(&mut d_tr_a4)
                .arg(&mut d_tr_a6)
                .arg(&mut d_tr_a8)
                .arg(&mut d_tr_b4)
                .arg(&mut d_tr_b8)
                .launch(cfg)?;
        }

        let h_tr_a4 = self.stream.clone_dtoh(&d_tr_a4)?;
        let h_tr_a6 = self.stream.clone_dtoh(&d_tr_a6)?;
        let h_tr_a8 = self.stream.clone_dtoh(&d_tr_a8)?;
        let h_tr_b4 = self.stream.clone_dtoh(&d_tr_b4)?;
        let h_tr_b8 = self.stream.clone_dtoh(&d_tr_b8)?;

        let mut out = Vec::with_capacity(batch_size);
        for i in 0..batch_size {
            out.push(HashimotoB4B8Summary {
                tr_a4: h_tr_a4[i],
                tr_a6: h_tr_a6[i],
                tr_a8: h_tr_a8[i],
                tr_b4: h_tr_b4[i],
                tr_b8: h_tr_b8[i],
            });
        }
        Ok(out)
    }

    /// Evaluates `k_enumerate_n16_bicyclic_walks` (Theorem 2.5a) on a batch of `batch_size` graphs.
    pub fn eval_n16_bicyclic_walks_batch(
        &self,
        flat_adj: &[u64],
        n: usize,
        batch_size: usize,
    ) -> Result<Vec<N16WalkDecompositionSummary>> {
        assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
        let d_adj = self.stream.clone_htod(flat_adj)?;
        let mut d_walks_c16 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_walks_c8 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_walks_c4 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_walks_dumbbell = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_walks_theta = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b16_total = self.stream.alloc_zeros::<u64>(batch_size)?;

        let n_u32 = n as u32;
        let b_u32 = batch_size as u32;
        let cfg = self.launch_cfg(b_u32);

        unsafe {
            self.stream
                .launch_builder(&self.fn_enumerate_n16_bicyclic_walks)
                .arg(&d_adj)
                .arg(&n_u32)
                .arg(&b_u32)
                .arg(&mut d_walks_c16)
                .arg(&mut d_walks_c8)
                .arg(&mut d_walks_c4)
                .arg(&mut d_walks_dumbbell)
                .arg(&mut d_walks_theta)
                .arg(&mut d_tr_b16_total)
                .launch(cfg)?;
        }

        let h_walks_c16 = self.stream.clone_dtoh(&d_walks_c16)?;
        let h_walks_c8 = self.stream.clone_dtoh(&d_walks_c8)?;
        let h_walks_c4 = self.stream.clone_dtoh(&d_walks_c4)?;
        let h_walks_dumbbell = self.stream.clone_dtoh(&d_walks_dumbbell)?;
        let h_walks_theta = self.stream.clone_dtoh(&d_walks_theta)?;
        let h_tr_b16_total = self.stream.clone_dtoh(&d_tr_b16_total)?;

        let mut out = Vec::with_capacity(batch_size);
        for i in 0..batch_size {
            out.push(N16WalkDecompositionSummary {
                walks_c16: h_walks_c16[i],
                walks_c8: h_walks_c8[i],
                walks_c4: h_walks_c4[i],
                walks_dumbbell: h_walks_dumbbell[i],
                walks_theta: h_walks_theta[i],
                tr_b16_total: h_tr_b16_total[i],
            });
        }
        Ok(out)
    }

    /// Evaluates `k_eval_b10_hankel_pinch` (Theorem 2.6) on a batch of `batch_size` graphs.
    pub fn eval_b10_hankel_pinch_batch(
        &self,
        flat_adj: &[u64],
        n: usize,
        batch_size: usize,
    ) -> Result<Vec<B10HankelPinchSummary>> {
        assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
        let d_adj = self.stream.clone_htod(flat_adj)?;
        let mut d_c6 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c10 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_a10 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b10 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c10_lo = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c10_hi = self.stream.alloc_zeros::<u64>(batch_size)?;

        let n_u32 = n as u32;
        let b_u32 = batch_size as u32;
        let cfg = self.launch_cfg(b_u32);

        unsafe {
            self.stream
                .launch_builder(&self.fn_eval_b10_hankel_pinch)
                .arg(&d_adj)
                .arg(&n_u32)
                .arg(&b_u32)
                .arg(&mut d_c6)
                .arg(&mut d_c10)
                .arg(&mut d_tr_a10)
                .arg(&mut d_tr_b10)
                .arg(&mut d_c10_lo)
                .arg(&mut d_c10_hi)
                .launch(cfg)?;
        }

        let h_c6 = self.stream.clone_dtoh(&d_c6)?;
        let h_c10 = self.stream.clone_dtoh(&d_c10)?;
        let h_tr_a10 = self.stream.clone_dtoh(&d_tr_a10)?;
        let h_tr_b10 = self.stream.clone_dtoh(&d_tr_b10)?;
        let h_c10_lo = self.stream.clone_dtoh(&d_c10_lo)?;
        let h_c10_hi = self.stream.clone_dtoh(&d_c10_hi)?;

        let mut out = Vec::with_capacity(batch_size);
        for i in 0..batch_size {
            out.push(B10HankelPinchSummary {
                c6: h_c6[i],
                c10: h_c10[i],
                tr_a10: h_tr_a10[i],
                tr_b10: h_tr_b10[i],
                c10_hankel_lower_bound: h_c10_lo[i],
                c10_hankel_upper_bound: h_c10_hi[i],
            });
        }
        Ok(out)
    }

    /// Evaluates `k_eval_b12_deg2_bessel_pinch` (Theorem 2.7) on a batch of `batch_size` graphs.
    pub fn eval_b12_deg2_bessel_pinch_batch(
        &self,
        flat_adj: &[u64],
        n: usize,
        batch_size: usize,
    ) -> Result<Vec<B12Deg2BesselSummary>> {
        assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
        let d_adj = self.stream.clone_htod(flat_adj)?;
        let mut d_c6 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c10 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c12 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tau1 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tau3 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_a12 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b12 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b16_3d = self.stream.alloc_zeros::<u64>(batch_size)?;

        let n_u32 = n as u32;
        let b_u32 = batch_size as u32;
        let cfg = self.launch_cfg(b_u32);

        unsafe {
            self.stream
                .launch_builder(&self.fn_eval_b12_deg2_bessel_pinch)
                .arg(&d_adj)
                .arg(&n_u32)
                .arg(&b_u32)
                .arg(&mut d_c6)
                .arg(&mut d_c10)
                .arg(&mut d_c12)
                .arg(&mut d_tau1)
                .arg(&mut d_tau3)
                .arg(&mut d_tr_a12)
                .arg(&mut d_tr_b12)
                .arg(&mut d_tr_b16_3d)
                .launch(cfg)?;
        }

        let h_c6 = self.stream.clone_dtoh(&d_c6)?;
        let h_c10 = self.stream.clone_dtoh(&d_c10)?;
        let h_c12 = self.stream.clone_dtoh(&d_c12)?;
        let h_tau1 = self.stream.clone_dtoh(&d_tau1)?;
        let h_tau3 = self.stream.clone_dtoh(&d_tau3)?;
        let h_tr_a12 = self.stream.clone_dtoh(&d_tr_a12)?;
        let h_tr_b12 = self.stream.clone_dtoh(&d_tr_b12)?;
        let h_tr_b16_3d = self.stream.clone_dtoh(&d_tr_b16_3d)?;

        let mut out = Vec::with_capacity(batch_size);
        for i in 0..batch_size {
            out.push(B12Deg2BesselSummary {
                c6: h_c6[i],
                c10: h_c10[i],
                c12: h_c12[i],
                tau1: h_tau1[i],
                tau3: h_tau3[i],
                tr_a12: h_tr_a12[i],
                tr_b12: h_tr_b12[i],
                tr_b16_3d_bessel_floor: h_tr_b16_3d[i],
            });
        }
        Ok(out)
    }

    /// Evaluates `k_eval_b14_deg3_bessel_pinch` (Theorem 2.8) on a batch of `batch_size` graphs.
    pub fn eval_b14_deg3_bessel_pinch_batch(
        &self,
        flat_adj: &[u64],
        n: usize,
        batch_size: usize,
    ) -> Result<Vec<B14Deg3BesselSummary>> {
        assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
        let d_adj = self.stream.clone_htod(flat_adj)?;
        let mut d_c6 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c10 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c12 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_c14 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_beta1 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_a14 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b14 = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b16_4d = self.stream.alloc_zeros::<u64>(batch_size)?;

        let n_u32 = n as u32;
        let b_u32 = batch_size as u32;
        let cfg = self.launch_cfg(b_u32);

        unsafe {
            self.stream
                .launch_builder(&self.fn_eval_b14_deg3_bessel_pinch)
                .arg(&d_adj)
                .arg(&n_u32)
                .arg(&b_u32)
                .arg(&mut d_c6)
                .arg(&mut d_c10)
                .arg(&mut d_c12)
                .arg(&mut d_c14)
                .arg(&mut d_beta1)
                .arg(&mut d_tr_a14)
                .arg(&mut d_tr_b14)
                .arg(&mut d_tr_b16_4d)
                .launch(cfg)?;
        }

        let h_c6 = self.stream.clone_dtoh(&d_c6)?;
        let h_c10 = self.stream.clone_dtoh(&d_c10)?;
        let h_c12 = self.stream.clone_dtoh(&d_c12)?;
        let h_c14 = self.stream.clone_dtoh(&d_c14)?;
        let h_beta1 = self.stream.clone_dtoh(&d_beta1)?;
        let h_tr_a14 = self.stream.clone_dtoh(&d_tr_a14)?;
        let h_tr_b14 = self.stream.clone_dtoh(&d_tr_b14)?;
        let h_tr_b16_4d = self.stream.clone_dtoh(&d_tr_b16_4d)?;

        let mut out = Vec::with_capacity(batch_size);
        for i in 0..batch_size {
            out.push(B14Deg3BesselSummary {
                c6: h_c6[i],
                c10: h_c10[i],
                c12: h_c12[i],
                c14: h_c14[i],
                beta1: h_beta1[i],
                tr_a14: h_tr_a14[i],
                tr_b14: h_tr_b14[i],
                tr_b16_4d_bessel_floor: h_tr_b16_4d[i],
            });
        }
        Ok(out)
    }

    /// Evaluates `k_eval_multi_octave_discharging` (Theorems 2.9 & 2.10 + Corollary 2.11) on a batch of `batch_size` graphs.
    pub fn eval_multi_octave_discharging_batch(
        &self,
        flat_adj: &[u64],
        n: usize,
        batch_size: usize,
    ) -> Result<Vec<MultiOctaveDischargingSummary>> {
        assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
        let d_adj = self.stream.clone_htod(flat_adj)?;
        let mut d_theta_ok = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_prism_ok = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_perron_k = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_perron_m = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_perron_margin = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_tr_b16_4d = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_n16_ub = self.stream.alloc_zeros::<u64>(batch_size)?;
        let mut d_closure_ok = self.stream.alloc_zeros::<u64>(batch_size)?;

        let n_u32 = n as u32;
        let b_u32 = batch_size as u32;
        let cfg = self.launch_cfg(b_u32);

        unsafe {
            self.stream
                .launch_builder(&self.fn_eval_multi_octave_discharging)
                .arg(&d_adj)
                .arg(&n_u32)
                .arg(&b_u32)
                .arg(&mut d_theta_ok)
                .arg(&mut d_prism_ok)
                .arg(&mut d_perron_k)
                .arg(&mut d_perron_m)
                .arg(&mut d_perron_margin)
                .arg(&mut d_tr_b16_4d)
                .arg(&mut d_n16_ub)
                .arg(&mut d_closure_ok)
                .launch(cfg)?;
        }

        let h_theta_ok = self.stream.clone_dtoh(&d_theta_ok)?;
        let h_prism_ok = self.stream.clone_dtoh(&d_prism_ok)?;
        let h_perron_k = self.stream.clone_dtoh(&d_perron_k)?;
        let h_perron_m = self.stream.clone_dtoh(&d_perron_m)?;
        let h_perron_margin = self.stream.clone_dtoh(&d_perron_margin)?;
        let h_tr_b16_4d = self.stream.clone_dtoh(&d_tr_b16_4d)?;
        let h_n16_ub = self.stream.clone_dtoh(&d_n16_ub)?;
        let h_closure_ok = self.stream.clone_dtoh(&d_closure_ok)?;

        let mut out = Vec::with_capacity(batch_size);
        for i in 0..batch_size {
            out.push(MultiOctaveDischargingSummary {
                theta_conservation_verified: h_theta_ok[i],
                prism_conservation_verified: h_prism_ok[i],
                perron_dominance_k: h_perron_k[i],
                perron_dominance_m: h_perron_m[i],
                perron_dominance_margin: h_perron_margin[i],
                tr_b16_4d_bessel_floor: h_tr_b16_4d[i],
                n16_discharged_upper_bound: h_n16_ub[i],
                dyadic_closure_certified: h_closure_ok[i],
            });
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::{
        build_canonical_cage_and_snark_library, build_phase3_commutator_frontier_library,
    };
    use crate::verify_cpu::{
        eval_b10_hankel_pinch_cpu, eval_b12_deg2_bessel_pinch_cpu, eval_b14_deg3_bessel_pinch_cpu,
        eval_dyadic_cycles_cpu, eval_hashimoto_b4_b8_cpu, eval_multi_octave_discharging_cpu,
        eval_n16_bicyclic_walks_cpu,
    };

    #[test]
    #[ignore = "requires CUDA GPU"]
    fn test_gpu_cpu_dyadic_cycles_parity_on_cages_and_snarks() {
        let gpu = Gpu::init(0, "compute_86").expect("failed to init CUDA GPU");
        let lib = build_canonical_cage_and_snark_library();

        for g in &lib {
            if g.n > 20 {
                continue;
            }
            let padded = g.padded_adj_64();
            let cpu_res = eval_dyadic_cycles_cpu(&padded, g.n, 1);
            let gpu_res = gpu
                .eval_dyadic_cycles_batch(&padded, g.n, 1)
                .expect("GPU batch evaluation failed");
            assert_eq!(
                cpu_res[0], gpu_res[0],
                "CPU vs GPU dyadic cycle parity mismatch on {}",
                g.name
            );
            let cpu_hash = eval_hashimoto_b4_b8_cpu(&padded, g.n, 1);
            let gpu_hash = gpu
                .eval_hashimoto_b4_b8_batch(&padded, g.n, 1)
                .expect("GPU Hashimoto evaluation failed");
            assert_eq!(
                cpu_hash[0], gpu_hash[0],
                "CPU vs GPU Hashimoto trace parity mismatch on {}",
                g.name
            );
            let cpu_n16 = eval_n16_bicyclic_walks_cpu(&padded, g.n, 1);
            let gpu_n16 = gpu
                .eval_n16_bicyclic_walks_batch(&padded, g.n, 1)
                .expect("GPU N_16 walk decomposition failed");
            assert_eq!(
                cpu_n16[0], gpu_n16[0],
                "CPU vs GPU N_16 topological walk decomposition parity mismatch on {}",
                g.name
            );
            let cpu_b10 = eval_b10_hankel_pinch_cpu(&padded, g.n, 1);
            let gpu_b10 = gpu
                .eval_b10_hankel_pinch_batch(&padded, g.n, 1)
                .expect("GPU B^10 Hankel pinch failed");
            assert_eq!(
                cpu_b10[0], gpu_b10[0],
                "CPU vs GPU B^10 Hankel pinch parity mismatch on {}",
                g.name
            );
            let cpu_b12 = eval_b12_deg2_bessel_pinch_cpu(&padded, g.n, 1);
            let gpu_b12 = gpu
                .eval_b12_deg2_bessel_pinch_batch(&padded, g.n, 1)
                .expect("GPU B^12 Bessel pinch failed");
            assert_eq!(
                cpu_b12[0], gpu_b12[0],
                "CPU vs GPU B^12 Bessel pinch parity mismatch on {}",
                g.name
            );
            let cpu_b14 = eval_b14_deg3_bessel_pinch_cpu(&padded, g.n, 1);
            let gpu_b14 = gpu
                .eval_b14_deg3_bessel_pinch_batch(&padded, g.n, 1)
                .expect("GPU B^14 Bessel pinch failed");
            assert_eq!(
                cpu_b14[0], gpu_b14[0],
                "CPU vs GPU B^14 Bessel pinch parity mismatch on {}",
                g.name
            );
            let cpu_p9 = eval_multi_octave_discharging_cpu(&padded, g.n, 1);
            let gpu_p9 = gpu
                .eval_multi_octave_discharging_batch(&padded, g.n, 1)
                .expect("GPU Phase 9 multi-octave discharging failed");
            assert_eq!(
                cpu_p9[0], gpu_p9[0],
                "CPU vs GPU Phase 9 multi-octave discharging parity mismatch on {}",
                g.name
            );
        }

        let phase3 = build_phase3_commutator_frontier_library();
        for g in &phase3 {
            let padded = g.padded_adj_64();
            let cpu_hash = eval_hashimoto_b4_b8_cpu(&padded, g.n, 1);
            let gpu_hash = gpu
                .eval_hashimoto_b4_b8_batch(&padded, g.n, 1)
                .expect("GPU Hashimoto evaluation failed on phase3");
            assert_eq!(
                cpu_hash[0], gpu_hash[0],
                "CPU vs GPU Hashimoto trace parity mismatch on {}",
                g.name
            );
            let cpu_n16 = eval_n16_bicyclic_walks_cpu(&padded, g.n, 1);
            let gpu_n16 = gpu
                .eval_n16_bicyclic_walks_batch(&padded, g.n, 1)
                .expect("GPU N_16 walk decomposition failed on phase3");
            assert_eq!(
                cpu_n16[0], gpu_n16[0],
                "CPU vs GPU N_16 topological walk decomposition parity mismatch on {}",
                g.name
            );
            let cpu_b10 = eval_b10_hankel_pinch_cpu(&padded, g.n, 1);
            let gpu_b10 = gpu
                .eval_b10_hankel_pinch_batch(&padded, g.n, 1)
                .expect("GPU B^10 Hankel pinch failed on phase3");
            assert_eq!(
                cpu_b10[0], gpu_b10[0],
                "CPU vs GPU B^10 Hankel pinch parity mismatch on {}",
                g.name
            );
            let cpu_b12 = eval_b12_deg2_bessel_pinch_cpu(&padded, g.n, 1);
            let gpu_b12 = gpu
                .eval_b12_deg2_bessel_pinch_batch(&padded, g.n, 1)
                .expect("GPU B^12 Bessel pinch failed on phase3");
            assert_eq!(
                cpu_b12[0], gpu_b12[0],
                "CPU vs GPU B^12 Bessel pinch parity mismatch on {}",
                g.name
            );
            let cpu_b14 = eval_b14_deg3_bessel_pinch_cpu(&padded, g.n, 1);
            let gpu_b14 = gpu
                .eval_b14_deg3_bessel_pinch_batch(&padded, g.n, 1)
                .expect("GPU B^14 Bessel pinch failed on phase3");
            assert_eq!(
                cpu_b14[0], gpu_b14[0],
                "CPU vs GPU B^14 Bessel pinch parity mismatch on {}",
                g.name
            );
            let cpu_p9 = eval_multi_octave_discharging_cpu(&padded, g.n, 1);
            let gpu_p9 = gpu
                .eval_multi_octave_discharging_batch(&padded, g.n, 1)
                .expect("GPU Phase 9 multi-octave discharging failed on phase3");
            assert_eq!(
                cpu_p9[0], gpu_p9[0],
                "CPU vs GPU Phase 9 multi-octave discharging parity mismatch on {}",
                g.name
            );
        }
    }
}
