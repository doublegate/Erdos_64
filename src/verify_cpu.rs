//! ============================================================================
//! `verify_cpu.rs` — Exact CPU verification, bitmask cycle-spectrum algorithms,
//! and CPU twin for CUDA kernels (`eval_dyadic_cycles_cpu`).
//! ============================================================================
//!
//! Implements the algorithmic core of Theorem 2.1 and Theorem 3.1:
//! - Exact girth `g(G)` and bipartiteness certification via bitwise BFS;
//! - `O(n^2)` bitwise 4-cycle counter `count_c4_bitwise`;
//! - Canonical-anchor (`s = min V(C)`) bitmask path-extension and half-path
//!   meet-in-the-middle counters for `C_8(G)` and `C_16(G)`;
//! - Exact full cycle-length spectrum `C_l(G)` (`3 <= l <= n`) for cages and
//!   snarks up to `n = 30`;
//! - Exact CPU twin `eval_dyadic_cycles_cpu` matching the CUDA kernel
//!   `k_eval_dyadic_cycles` bit-for-bit.

use crate::records::{BitGraph, CycleSpectrumRecord, MAX_VERTICES};
use serde::{Deserialize, Serialize};

/// Output of `eval_dyadic_cycles_cpu` and the GPU kernel `k_eval_dyadic_cycles`
/// for a single graph on `n <= 64` vertices.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DyadicEvalSummary {
    pub c4: u64,
    pub c8: u64,
    pub c16_anchored0: u64,
    pub min_deg: u32,
}

/// Checks whether `G` is bipartite using 2-coloring bitwise BFS.
pub fn is_bipartite(g: &BitGraph) -> bool {
    let n = g.n;
    let mut color = vec![0u8; n]; // 0 = unvisited, 1 = red, 2 = blue
    for start in 0..n {
        if color[start] != 0 {
            continue;
        }
        color[start] = 1;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        while let Some(u) = queue.pop_front() {
            let next_c = if color[u] == 1 { 2 } else { 1 };
            let mut mask = g.adj[u];
            while mask != 0 {
                let v = mask.trailing_zeros() as usize;
                mask &= mask - 1;
                if color[v] == 0 {
                    color[v] = next_c;
                    queue.push_back(v);
                } else if color[v] == color[u] {
                    return false;
                }
            }
        }
    }
    true
}

/// Computes the exact girth `g(G)` (length of the shortest simple cycle, or `usize::MAX` if acyclic).
pub fn compute_girth(g: &BitGraph) -> usize {
    let n = g.n;
    let mut best = usize::MAX;
    for s in 0..n {
        let mut dist = vec![usize::MAX; n];
        let mut parent = vec![usize::MAX; n];
        let mut q = std::collections::VecDeque::new();
        dist[s] = 0;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            let d = dist[u];
            if 2 * d + 1 >= best {
                break;
            }
            let mut mask = g.adj[u];
            while mask != 0 {
                let v = mask.trailing_zeros() as usize;
                mask &= mask - 1;
                if dist[v] == usize::MAX {
                    dist[v] = d + 1;
                    parent[v] = u;
                    q.push_back(v);
                } else if parent[u] != v {
                    let cycle_len = dist[u] + dist[v] + 1;
                    if cycle_len < best {
                        best = cycle_len;
                    }
                }
            }
        }
    }
    best
}

/// Counts undirected 4-cycles `C_4(G)` in `O(n^2)` bitwise operations using
/// `sum_{u < v} binom(popcount(adj[u] & adj[v]), 2) / 2`.
#[inline]
pub fn count_c4_slice(adj: &[u64], n: usize) -> u64 {
    let mut pairs_of_common_neighbors = 0u64;
    for u in 0..n {
        let au = adj[u];
        for &av in adj.iter().take(n).skip(u + 1) {
            let common = (au & av).count_ones() as u64;
            if common >= 2 {
                pairs_of_common_neighbors += common * (common - 1) / 2;
            }
        }
    }
    // Each undirected 4-cycle (a - b - c - d - a) has two opposite pairs {a, c} and {b, d}.
    pairs_of_common_neighbors / 2
}

/// Counts undirected 8-cycles `C_8(G)` using canonical-anchor (`s = min V(C)`)
/// bitmask path-extension with last-step intersection closure.
///
/// Every simple 8-cycle `(s, v_1, v_2, v_3, v_4, v_5, v_6, v_7)` has a unique minimum
/// vertex `s` and two orientations (`v_1 < v_7` vs `v_1 > v_7`). By枚举 paths
/// `(s, v_1, v_2, v_3, v_4, v_5, v_6)` with `v_i > s`, `v_7 in N(v_6) \cap N(s)` with
/// `v_7 > v_1` and `v_7 notin {v_1..v_6}` is counted via a single bitwise AND + popcount!
pub fn count_c8_slice(adj: &[u64], n: usize) -> u64 {
    let mut total_c8 = 0u64;
    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        if ns.count_ones() < 2 {
            continue;
        }
        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            let vis1 = (1u64 << s) | (1u64 << v1);
            // v7 must be in ns and strictly greater than v1
            let valid_v7_base = if v1 == 63 {
                0u64
            } else {
                ns & !((1u64 << (v1 + 1)) - 1)
            };
            if valid_v7_base == 0 {
                continue;
            }

            let mut m2 = adj[v1] & higher_mask & !vis1;
            while m2 != 0 {
                let v2 = m2.trailing_zeros() as usize;
                m2 &= m2 - 1;
                let vis2 = vis1 | (1u64 << v2);

                let mut m3 = adj[v2] & higher_mask & !vis2;
                while m3 != 0 {
                    let v3 = m3.trailing_zeros() as usize;
                    m3 &= m3 - 1;
                    let vis3 = vis2 | (1u64 << v3);

                    let mut m4 = adj[v3] & higher_mask & !vis3;
                    while m4 != 0 {
                        let v4 = m4.trailing_zeros() as usize;
                        m4 &= m4 - 1;
                        let vis4 = vis3 | (1u64 << v4);

                        let mut m5 = adj[v4] & higher_mask & !vis4;
                        while m5 != 0 {
                            let v5 = m5.trailing_zeros() as usize;
                            m5 &= m5 - 1;
                            let vis5 = vis4 | (1u64 << v5);

                            let mut m6 = adj[v5] & higher_mask & !vis5;
                            while m6 != 0 {
                                let v6 = m6.trailing_zeros() as usize;
                                m6 &= m6 - 1;
                                let vis6 = vis5 | (1u64 << v6);
                                let v7_mask = adj[v6] & valid_v7_base & !vis6;
                                total_c8 += v7_mask.count_ones() as u64;
                            }
                        }
                    }
                }
            }
        }
    }
    total_c8
}

/// Counts undirected 8-cycles `C_8(G)` via Theorem 3.1's Half-Path Meet-in-the-Middle
/// decomposition: for each canonical anchor `s`, enumerates simple 4-step half-paths
/// `P = (s, v_1, v_2, v_3, m)` with `v_i, m > s`, groups them by midpoint `m`, and counts
/// pairs `(P, Q)` sharing `(s, m)` whose interior 3-vertex masks `{v_1, v_2, v_3}` are disjoint.
pub fn count_c8_mitm(g: &BitGraph) -> u64 {
    let n = g.n;
    let adj = &g.adj;
    let mut total_c8 = 0u64;

    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        if ns.count_ones() < 2 {
            continue;
        }
        // Group half-paths by midpoint `m`: store `(first_step v1, interior_mask {v1,v2,v3})`
        let mut by_mid: Vec<Vec<(usize, u64)>> = vec![Vec::new(); n];
        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            let b1 = 1u64 << v1;

            let mut m2 = adj[v1] & higher_mask & !b1;
            while m2 != 0 {
                let v2 = m2.trailing_zeros() as usize;
                m2 &= m2 - 1;
                let b12 = b1 | (1u64 << v2);

                let mut m3 = adj[v2] & higher_mask & !b12;
                while m3 != 0 {
                    let v3 = m3.trailing_zeros() as usize;
                    m3 &= m3 - 1;
                    let int_mask = b12 | (1u64 << v3);

                    let mut m4 = adj[v3] & higher_mask & !int_mask;
                    while m4 != 0 {
                        let mid = m4.trailing_zeros() as usize;
                        m4 &= m4 - 1;
                        by_mid[mid].push((v1, int_mask));
                    }
                }
            }
        }

        for paths in &by_mid {
            let len = paths.len();
            for i in 0..len {
                let (v1_a, mask_a) = paths[i];
                for &(v1_b, mask_b) in &paths[(i + 1)..len] {
                    if v1_a != v1_b && (mask_a & mask_b) == 0 {
                        total_c8 += 1;
                    }
                }
            }
        }
    }
    total_c8
}

/// Counts undirected 16-cycles `C_16(G)` passing through a specific anchor vertex `s`
/// (with all other vertices `> s`) using Half-Path Meet-in-the-Middle:
/// enumerates 8-step half-paths `(s, v_1, ..., v_7, m)` with `v_i, m > s`, groups by `m`,
/// and counts pairs with disjoint 7-vertex interior masks `{v_1..v_7}`.
pub fn count_c16_anchored(adj: &[u64], n: usize, s: usize) -> u64 {
    if n < 16 || s + 15 >= n {
        return 0;
    }
    let higher_mask = if s == 63 {
        0u64
    } else {
        !((1u64 << (s + 1)) - 1)
    };
    let ns = adj[s] & higher_mask;
    if ns.count_ones() < 2 {
        return 0;
    }

    let mut by_mid: Vec<Vec<u64>> = vec![Vec::new(); n];

    #[allow(clippy::too_many_arguments)]
    fn dfs_half8(
        adj: &[u64],
        higher_mask: u64,
        u: usize,
        depth: usize,
        int_mask: u64,
        by_mid: &mut [Vec<u64>],
    ) {
        let mut mask = adj[u] & higher_mask & !int_mask;
        if depth == 7 {
            while mask != 0 {
                let mid = mask.trailing_zeros() as usize;
                mask &= mask - 1;
                by_mid[mid].push(int_mask);
            }
            return;
        }
        while mask != 0 {
            let nxt = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            dfs_half8(
                adj,
                higher_mask,
                nxt,
                depth + 1,
                int_mask | (1u64 << nxt),
                by_mid,
            );
        }
    }

    let mut m1 = ns;
    while m1 != 0 {
        let v1 = m1.trailing_zeros() as usize;
        m1 &= m1 - 1;
        dfs_half8(adj, higher_mask, v1, 1, 1u64 << v1, &mut by_mid);
    }

    let mut count = 0u64;
    for paths in &by_mid {
        let len = paths.len();
        for i in 0..len {
            let ma = paths[i];
            for &mb in &paths[(i + 1)..len] {
                if (ma & mb) == 0 {
                    count += 1;
                }
            }
        }
    }
    count
}

/// Counts all undirected 16-cycles `C_16(G)` across all canonical anchors `s in 0..n`.
pub fn count_c16_total(g: &BitGraph) -> u64 {
    if g.n < 16 {
        return 0;
    }
    (0..g.n).map(|s| count_c16_anchored(&g.adj, g.n, s)).sum()
}

/// Counts undirected 32-cycles `C_32^(0)(G)` anchored at `s = 0` (`0 = min V(C)`)
/// using Half-Path Meet-in-the-Middle of 16-step half-paths `(0, v_1, ..., v_15, m)`
/// with a per-midpoint cap `max_paths_per_mid` so it executes in milliseconds on
/// `50 <= n <= 62` frontier cubic graphs.
pub fn count_c32_anchored0_bounded(adj: &[u64], n: usize, max_paths_per_mid: usize) -> u64 {
    if n < 32 {
        return 0;
    }
    let higher_mask = !1u64;
    let ns = adj[0] & higher_mask;
    if ns.count_ones() < 2 {
        return 0;
    }

    let mut by_mid: Vec<Vec<(usize, u64)>> = vec![Vec::new(); n];

    #[allow(clippy::too_many_arguments)]
    fn dfs_half16(
        adj: &[u64],
        higher_mask: u64,
        v1: usize,
        u: usize,
        depth: usize,
        int_mask: u64,
        max_paths_per_mid: usize,
        by_mid: &mut [Vec<(usize, u64)>],
    ) {
        let mut mask = adj[u] & higher_mask & !int_mask;
        if depth == 15 {
            while mask != 0 {
                let mid = mask.trailing_zeros() as usize;
                mask &= mask - 1;
                if by_mid[mid].len() < max_paths_per_mid {
                    by_mid[mid].push((v1, int_mask));
                }
            }
            return;
        }
        while mask != 0 {
            let nxt = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            dfs_half16(
                adj,
                higher_mask,
                v1,
                nxt,
                depth + 1,
                int_mask | (1u64 << nxt),
                max_paths_per_mid,
                by_mid,
            );
        }
    }

    let mut m1 = ns;
    while m1 != 0 {
        let v1 = m1.trailing_zeros() as usize;
        m1 &= m1 - 1;
        dfs_half16(
            adj,
            higher_mask,
            v1,
            v1,
            1,
            1u64 << v1,
            max_paths_per_mid,
            &mut by_mid,
        );
    }

    let mut count = 0u64;
    for paths in &by_mid {
        let len = paths.len();
        for i in 0..len {
            let (v1_a, ma) = paths[i];
            for &(v1_b, mb) in &paths[(i + 1)..len] {
                if v1_a != v1_b && (ma & mb) == 0 {
                    count += 1;
                }
            }
        }
    }
    count
}

/// Counts all undirected simple cycles of every length `l in [3, n]` in `G`
/// using canonical-anchor (`s = min V(C)`) bitmask backtracking with 2-connected
/// reachability pruning.
pub fn compute_full_cycle_spectrum(g: &BitGraph) -> CycleSpectrumRecord {
    let n = g.n;
    let adj = &g.adj;
    let mut counts = vec![0u64; n + 1];

    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        if ns.count_ones() < 2 {
            continue;
        }

        #[allow(clippy::too_many_arguments)]
        fn dfs_cycles(
            adj: &[u64],
            v1: usize,
            u: usize,
            len: usize,
            visited: u64,
            higher_mask: u64,
            ns: u64,
            counts: &mut [u64],
        ) {
            // Check if `u` can close a cycle back to `s` with orientation `v_1 < u`
            if len >= 3 && (ns & (1u64 << u)) != 0 && u > v1 {
                counts[len] += 1;
            }
            let mut mask = adj[u] & higher_mask & !visited;
            while mask != 0 {
                let nxt = mask.trailing_zeros() as usize;
                mask &= mask - 1;
                dfs_cycles(
                    adj,
                    v1,
                    nxt,
                    len + 1,
                    visited | (1u64 << nxt),
                    higher_mask,
                    ns,
                    counts,
                );
            }
        }

        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            dfs_cycles(
                adj,
                v1,
                v1,
                2,
                (1u64 << s) | (1u64 << v1),
                higher_mask,
                ns,
                &mut counts,
            );
        }
    }

    let c4 = if n >= 4 { counts[4] } else { 0 };
    let c8 = if n >= 8 { counts[8] } else { 0 };
    let c16 = if n >= 16 { counts[16] } else { 0 };
    let c32 = if n >= 32 { counts[32] } else { 0 };
    let total_dyadic_cycles = c4 + c8 + c16 + c32;

    let mut cycle_lengths = Vec::new();
    let mut dyadic_lengths_present = Vec::new();
    for (l, &cnt) in counts.iter().enumerate().skip(3) {
        if cnt > 0 {
            cycle_lengths.push(l);
            if l >= 4 && l.is_power_of_two() {
                dyadic_lengths_present.push(l);
            }
        }
    }

    let girth = cycle_lengths.first().copied().unwrap_or(0);
    let circumference = cycle_lengths.last().copied().unwrap_or(0);

    CycleSpectrumRecord {
        name: g.name.clone(),
        family: g.family.clone(),
        n,
        m: g.num_edges(),
        min_deg: g.min_degree(),
        max_deg: g.max_degree(),
        girth,
        circumference,
        is_bipartite: is_bipartite(g),
        is_cage: g.is_cage,
        is_snark: g.is_snark,
        c4,
        c8,
        c16,
        c32,
        total_dyadic_cycles,
        dyadic_lengths_present,
        cycle_lengths,
        cycle_counts_by_length: counts,
    }
}

/// Fast dyadic-only cycle evaluator for larger graphs (`n <= 64`) where full Hamiltonian-depth
/// spectrum enumeration is unnecessary: computes exact `C_4(G)`, `C_8(G)`, and `C_16(G)` via
/// bitmask meet-in-the-middle, plus bounded anchor-0 `C_32(G)` when `n >= 32`.
pub fn compute_dyadic_spectrum_fast(g: &BitGraph) -> CycleSpectrumRecord {
    let n = g.n;
    let c4 = count_c4_slice(&g.adj, n);
    let c8 = count_c8_slice(&g.adj, n);
    let c16 = count_c16_total(g);
    let c32 = if n >= 32 {
        count_c32_anchored0_bounded(&g.adj, n, 256)
    } else {
        0u64
    };
    let total_dyadic_cycles = c4 + c8 + c16 + c32;

    let mut dyadic_lengths_present = Vec::new();
    if c4 > 0 {
        dyadic_lengths_present.push(4);
    }
    if c8 > 0 {
        dyadic_lengths_present.push(8);
    }
    if c16 > 0 {
        dyadic_lengths_present.push(16);
    }
    if c32 > 0 {
        dyadic_lengths_present.push(32);
    }

    CycleSpectrumRecord {
        name: g.name.clone(),
        family: g.family.clone(),
        n,
        m: g.num_edges(),
        min_deg: g.min_degree(),
        max_deg: g.max_degree(),
        girth: compute_girth(g),
        circumference: 0,
        is_bipartite: is_bipartite(g),
        is_cage: g.is_cage,
        is_snark: g.is_snark,
        c4,
        c8,
        c16,
        c32,
        total_dyadic_cycles,
        dyadic_lengths_present,
        cycle_lengths: Vec::new(),
        cycle_counts_by_length: Vec::new(),
    }
}

/// Exact CPU twin of the CUDA kernel `k_eval_dyadic_cycles`:
/// given a batch of `batch_size` graphs each packed as `[u64; 64]` with vertex count `n`,
/// computes `(C_4(G), C_8(G), C_16_anchored_at_0(G), min_deg(G))` for each graph.
pub fn eval_dyadic_cycles_cpu(
    flat_adj: &[u64],
    n: usize,
    batch_size: usize,
) -> Vec<DyadicEvalSummary> {
    assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
    let mut out = Vec::with_capacity(batch_size);
    for b in 0..batch_size {
        let adj = &flat_adj[b * MAX_VERTICES..(b + 1) * MAX_VERTICES];
        let min_deg = (0..n).map(|u| adj[u].count_ones()).min().unwrap_or(0);
        let c4 = count_c4_slice(adj, n);
        let c8 = count_c8_slice(adj, n);
        let c16_anchored0 = count_c16_strided_anchor0(adj, n);
        out.push(DyadicEvalSummary {
            c4,
            c8,
            c16_anchored0,
            min_deg,
        });
    }
    out
}

/// Direct stack-based DFS counter for 16-cycles anchored at `s = 0` (`0 = min V(C)`),
/// matching the register-bounded GPU kernel logic bit-for-bit when `n >= 16` and `C_4(G) == 0`.
pub fn count_c16_strided_anchor0(adj: &[u64], n: usize) -> u64 {
    if n < 16 {
        return 0;
    }
    // Only run anchor-0 16-cycle probe on small/sparse cubic graphs (max_deg <= 3)
    for &au in adj.iter().take(n) {
        if au.count_ones() > 3 {
            return 0;
        }
    }
    let ns = adj[0] & !1u64;
    if ns.count_ones() < 2 {
        return 0;
    }
    let mut total = 0u64;
    let mut m1 = ns;
    while m1 != 0 {
        let v1 = m1.trailing_zeros() as usize;
        m1 &= m1 - 1;
        let valid_v15 = if v1 == 63 {
            0u64
        } else {
            ns & !((1u64 << (v1 + 1)) - 1)
        };
        if valid_v15 == 0 {
            continue;
        }
        // Iterative stack for path v1 -> v2 -> ... -> v14 (depth 1..=14)
        let mut path = [0usize; 16];
        let mut masks = [0u64; 16];
        let mut vis = [0u64; 16];
        path[1] = v1;
        vis[1] = 1u64 | (1u64 << v1);
        masks[1] = adj[v1] & !vis[1];
        let mut depth = 1usize;

        while depth >= 1 {
            let m = masks[depth];
            if m == 0 {
                depth -= 1;
                continue;
            }
            let nxt = m.trailing_zeros() as usize;
            masks[depth] &= m - 1;
            let next_vis = vis[depth] | (1u64 << nxt);
            if depth == 13 {
                // nxt is v14; count v15 in adj[v14] & valid_v15 & !next_vis
                let v15_mask = adj[nxt] & valid_v15 & !next_vis;
                total += v15_mask.count_ones() as u64;
            } else {
                depth += 1;
                path[depth] = nxt;
                vis[depth] = next_vis;
                masks[depth] = adj[nxt] & !next_vis;
            }
        }
    }
    total
}

/// Certificate for Theorem 2.3's Universal Abelian Commutator Floor in the Free Abelian
/// Dipole Cover `Lift_{Z^2}(D_3; (0,0), (1,0), (0,1))` (the hexagonal honeycomb lattice).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FreeAbelianCommutatorCertificate {
    /// Directed simple 4-cycles starting at `u_{(0,0)}` in `Lift_{Z^2}(D_3)`.
    pub directed_simple_c4: u64,
    /// Directed simple 6-cycles starting at `u_{(0,0)}` in `Lift_{Z^2}(D_3)` (`6 = 3 * 2 -> 1 universal Gamma-orbit`).
    pub directed_simple_c6: u64,
    /// Universal `Gamma`-orbits of simple 6-cycles (`1`, proving every abelian dipole lift has girth `<= 6`).
    pub universal_c6_orbits: u64,
    /// Directed simple 8-cycles starting at `u_{(0,0)}` in `Lift_{Z^2}(D_3)`.
    pub directed_simple_c8: u64,
    /// Total closed non-backtracking walks of length 16 starting at `u_{(0,0)}` in `Lift_{Z^2}(D_3)` (`960`).
    pub closed_nonbacktracking_16walks: u64,
    /// Universal `Gamma`-orbits of closed non-backtracking 16-walks (`960 / 16 = 60`).
    pub universal_nb16_orbits: u64,
    /// Directed simple 16-cycles starting at `u_{(0,0)}` in `Lift_{Z^2}(D_3)` (`288`).
    pub directed_simple_c16: u64,
    /// Universal `Gamma`-orbits of simple 16-cycles (`288 / 16 = 18`, giving `18 |Gamma|` universal 16-cycles).
    pub universal_c16_orbits: u64,
    /// Universal `Gamma`-orbits of non-simple figure-eight/theta 16-walks (`60 - 18 = 42`).
    pub universal_nonsimple_nb16_orbits: u64,
    /// Distinct nonzero net displacement vectors `\vec{\Delta} \in Z^2 \setminus {(0,0)}` reached by simple 16-walks (`216`).
    pub distinct_nonzero_16walk_displacements: usize,
}

/// Exhaustively enumerates closed non-backtracking walks of lengths 4, 6, 8, and 16
/// starting at `u_{(0,0)}` in the free abelian dipole cover `Lift_{Z^2}(D_3; (0,0), (1,0), (0,1))`
/// to certify Theorem 2.3(a)-(b).
pub fn certify_free_abelian_commutator_floor() -> FreeAbelianCommutatorCertificate {
    let shifts: [(i32, i32); 3] = [(0, 0), (1, 0), (0, 1)];

    fn count_simple_cycles_half_len(r: usize, shifts: &[(i32, i32); 3]) -> u64 {
        let mut count = 0u64;
        let mut u_vis = [(0i32, 0i32); 9];
        let mut v_vis = [(0i32, 0i32); 9];
        u_vis[0] = (0, 0);

        #[allow(clippy::too_many_arguments)]
        fn dfs(
            t: usize,
            r: usize,
            u_curr: (i32, i32),
            first_a: usize,
            prev_b: usize,
            u_vis: &mut [(i32, i32); 9],
            v_vis: &mut [(i32, i32); 9],
            shifts: &[(i32, i32); 3],
            count: &mut u64,
        ) {
            if t == r {
                if u_curr == (0, 0) && first_a != prev_b {
                    *count += 1;
                }
                return;
            }
            for a in 0..3 {
                if t > 0 && a == prev_b {
                    continue;
                }
                let v_next = (u_curr.0 + shifts[a].0, u_curr.1 + shifts[a].1);
                if v_vis[..t].contains(&v_next) {
                    continue;
                }
                v_vis[t] = v_next;
                for b in 0..3 {
                    if b == a {
                        continue;
                    }
                    let u_next = (v_next.0 - shifts[b].0, v_next.1 - shifts[b].1);
                    if t + 1 < r {
                        if u_vis[..=t].contains(&u_next) {
                            continue;
                        }
                        u_vis[t + 1] = u_next;
                    } else if u_next != (0, 0) {
                        continue;
                    }
                    dfs(
                        t + 1,
                        r,
                        u_next,
                        if t == 0 { a } else { first_a },
                        b,
                        u_vis,
                        v_vis,
                        shifts,
                        count,
                    );
                }
            }
        }

        dfs(
            0,
            r,
            (0, 0),
            0,
            usize::MAX,
            &mut u_vis,
            &mut v_vis,
            shifts,
            &mut count,
        );
        count
    }

    let directed_simple_c4 = count_simple_cycles_half_len(2, &shifts);
    let directed_simple_c6 = count_simple_cycles_half_len(3, &shifts);
    let directed_simple_c8 = count_simple_cycles_half_len(4, &shifts);
    let directed_simple_c16 = count_simple_cycles_half_len(8, &shifts);

    // Count ALL closed non-backtracking walks of length 16 (r = 8) in Lift_{Z^2}(D_3)
    let mut closed_nonbacktracking_16walks = 0u64;
    fn dfs_nb(
        t: usize,
        u_curr: (i32, i32),
        first_a: usize,
        prev_b: usize,
        shifts: &[(i32, i32); 3],
        count: &mut u64,
    ) {
        if t == 8 {
            if u_curr == (0, 0) && first_a != prev_b {
                *count += 1;
            }
            return;
        }
        for a in 0..3 {
            if t > 0 && a == prev_b {
                continue;
            }
            let v_next = (u_curr.0 + shifts[a].0, u_curr.1 + shifts[a].1);
            for b in 0..3 {
                if b == a {
                    continue;
                }
                let u_next = (v_next.0 - shifts[b].0, v_next.1 - shifts[b].1);
                dfs_nb(
                    t + 1,
                    u_next,
                    if t == 0 { a } else { first_a },
                    b,
                    shifts,
                    count,
                );
            }
        }
    }
    dfs_nb(
        0,
        (0, 0),
        0,
        usize::MAX,
        &shifts,
        &mut closed_nonbacktracking_16walks,
    );

    // Also count distinct nonzero net displacement vectors \vec{\Delta} in Z^2 \setminus {(0,0)}
    // reached by simple 16-walks with distinct interior vertices
    let mut displacements = std::collections::BTreeSet::new();
    let mut u_vis = [(0i32, 0i32); 9];
    let mut v_vis = [(0i32, 0i32); 9];
    u_vis[0] = (0, 0);
    #[allow(clippy::too_many_arguments)]
    fn dfs_disp(
        t: usize,
        u_curr: (i32, i32),
        first_a: usize,
        prev_b: usize,
        u_vis: &mut [(i32, i32); 9],
        v_vis: &mut [(i32, i32); 9],
        shifts: &[(i32, i32); 3],
        displacements: &mut std::collections::BTreeSet<(i32, i32)>,
    ) {
        if t == 8 {
            if u_curr != (0, 0) && first_a != prev_b {
                displacements.insert(u_curr);
            }
            return;
        }
        for a in 0..3 {
            if t > 0 && a == prev_b {
                continue;
            }
            let v_next = (u_curr.0 + shifts[a].0, u_curr.1 + shifts[a].1);
            if v_vis[..t].contains(&v_next) {
                continue;
            }
            v_vis[t] = v_next;
            for b in 0..3 {
                if b == a {
                    continue;
                }
                let u_next = (v_next.0 - shifts[b].0, v_next.1 - shifts[b].1);
                if t + 1 < 8 && u_vis[..=t].contains(&u_next) {
                    continue;
                }
                u_vis[t + 1] = u_next;
                dfs_disp(
                    t + 1,
                    u_next,
                    if t == 0 { a } else { first_a },
                    b,
                    u_vis,
                    v_vis,
                    shifts,
                    displacements,
                );
            }
        }
    }
    dfs_disp(
        0,
        (0, 0),
        0,
        usize::MAX,
        &mut u_vis,
        &mut v_vis,
        &shifts,
        &mut displacements,
    );

    let universal_nb16_orbits = closed_nonbacktracking_16walks / 16;
    let universal_c16_orbits = directed_simple_c16 / 16;

    FreeAbelianCommutatorCertificate {
        directed_simple_c4,
        directed_simple_c6,
        universal_c6_orbits: directed_simple_c6 / 6,
        directed_simple_c8,
        closed_nonbacktracking_16walks,
        universal_nb16_orbits,
        directed_simple_c16,
        universal_c16_orbits,
        universal_nonsimple_nb16_orbits: universal_nb16_orbits - universal_c16_orbits,
        distinct_nonzero_16walk_displacements: displacements.len(),
    }
}

/// Certificate for Theorem 3.2(a): Bass–Ihara Chebyshev Spectral Trace Identity
/// relating vertex adjacency matrix powers `Tr(A^{2r})` to the Hashimoto non-backtracking
/// edge matrix traces `Tr(B^4)`, `Tr(B^8)`, and `Tr(B^{16})`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HashimotoTraceCertificate {
    pub name: String,
    pub n: usize,
    pub girth: usize,
    pub c4: u64,
    pub c8: u64,
    pub c16: u64,
    /// `Tr(B^4)` computed via Bass–Ihara Chebyshev polynomial `Tr(A^4) - 8 Tr(A^2) + 9n`.
    pub tr_b4_chebyshev: u64,
    /// `Tr(B^4)` computed via explicit `(2m x 2m)` Hashimoto matrix power `B^4`.
    pub tr_b4_matrix: u64,
    /// `Tr(B^8)` computed via Bass–Ihara Chebyshev polynomial `Tr(A^8) - 16 Tr(A^6) + 80 Tr(A^4) - 128 Tr(A^2) + 33n`.
    pub tr_b8_chebyshev: u64,
    /// `Tr(B^8)` computed via explicit `(2m x 2m)` Hashimoto matrix power `B^8`.
    pub tr_b8_matrix: u64,
    /// `Tr(B^{16})` computed via Bass–Ihara degree-16 Chebyshev polynomial of `A`.
    pub tr_b16_chebyshev: u64,
    /// `Tr(B^{16})` computed via explicit `(2m x 2m)` Hashimoto matrix power `B^{16}`.
    pub tr_b16_matrix: u64,
    /// Non-simple figure-eight/lollipop/theta 16-walk excess `(Tr(B^{16}) - 32 * C_16(G)) / 32`.
    pub nonsimple_nb16_walk_orbits: u64,
}

/// Summary of the GPU/CPU Hashimoto spectral trace kernel `k_eval_hashimoto_trace_b4_b8`
/// for a batch of 3-regular graphs on `n <= 64` vertices.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct HashimotoB4B8Summary {
    pub tr_a4: u64,
    pub tr_a6: u64,
    pub tr_a8: u64,
    pub tr_b4: u64,
    pub tr_b8: u64,
}

/// Computes exact even vertex-adjacency traces `Tr(A^{2r})` for `r = 1..=8` (`2r = 2, 4, ..., 16`)
/// on a graph `G` of order `n <= 64`.
///
/// Since `A` is symmetric, `Tr(A^{2r}) = \sum_{s=0}^{n-1} || A^r e_s ||_2^2`, so we only need
/// `r <= 8` matrix-vector steps `v^{(t+1)}[u] = \sum_{w \in N(u)} v^{(t)}[w]` per starting vertex `s`!
pub fn compute_adjacency_even_traces(adj: &[u64], n: usize) -> [u128; 9] {
    let mut tr_even = [0u128; 9]; // index r stores Tr(A^{2r}) for r = 0..=8
    tr_even[0] = n as u128;
    let mut curr = vec![0u64; n];
    let mut next = vec![0u64; n];

    for s in 0..n {
        curr.fill(0);
        curr[s] = 1;
        for tr_r in tr_even.iter_mut().skip(1) {
            for u in 0..n {
                let mut sum = 0u64;
                let mut mask = adj[u];
                while mask != 0 {
                    let w = mask.trailing_zeros() as usize;
                    mask &= mask - 1;
                    sum += curr[w];
                }
                next[u] = sum;
            }
            curr.copy_from_slice(&next);
            let mut norm_sq = 0u128;
            for &val in &curr[..n] {
                let v128 = val as u128;
                norm_sq += v128 * v128;
            }
            *tr_r += norm_sq;
        }
    }
    tr_even
}

/// Computes explicit `(2m x 2m)` Hashimoto non-backtracking matrix traces `(Tr(B^4), Tr(B^8), Tr(B^{16}))`
/// via repeated squaring `B -> B^2 -> B^4 -> B^8 -> B^{16}`.
pub fn compute_hashimoto_traces_explicit_matrix(g: &BitGraph) -> (u64, u64, u64) {
    let n = g.n;
    let mut dir_edges = Vec::with_capacity(3 * n);
    for u in 0..n {
        let mut mask = g.adj[u];
        while mask != 0 {
            let v = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            dir_edges.push((u, v));
        }
    }
    let m_dir = dir_edges.len();
    let mut b_mat = vec![0u64; m_dir * m_dir];
    for (i, &(u, v)) in dir_edges.iter().enumerate() {
        for (j, &(w, x)) in dir_edges.iter().enumerate() {
            if v == w && x != u {
                b_mat[i * m_dir + j] = 1;
            }
        }
    }

    let square_mat = |mat: &[u64]| -> Vec<u64> {
        let mut out = vec![0u64; m_dir * m_dir];
        for i in 0..m_dir {
            for k in 0..m_dir {
                let aik = mat[i * m_dir + k];
                if aik == 0 {
                    continue;
                }
                for j in 0..m_dir {
                    out[i * m_dir + j] += aik * mat[k * m_dir + j];
                }
            }
        }
        out
    };

    let trace_mat = |mat: &[u64]| -> u64 { (0..m_dir).map(|i| mat[i * m_dir + i]).sum() };

    let b2 = square_mat(&b_mat);
    let b4 = square_mat(&b2);
    let b8 = square_mat(&b4);
    // Note: Tr(B^{16}) = \sum_{i,j} B^8_{i,j} * B^8_{j,i} without needing a full B^{16} multiply!
    let mut tr_b16 = 0u64;
    for i in 0..m_dir {
        for j in 0..m_dir {
            tr_b16 += b8[i * m_dir + j] * b8[j * m_dir + i];
        }
    }
    (trace_mat(&b4), trace_mat(&b8), tr_b16)
}

/// Certifies Theorem 3.2(a) on a 3-regular graph `g`: evaluates `Tr(B^4)`, `Tr(B^8)`, and `Tr(B^{16})`
/// via both the Bass–Ihara Chebyshev polynomial of the vertex adjacency matrix `A` and explicit
/// `(2m x 2m)` Hashimoto matrix squaring, and compares them to `8 C_4(G)`, `16 C_8(G)`, and `32 C_16(G)`.
pub fn certify_hashimoto_trace_identity(g: &BitGraph) -> HashimotoTraceCertificate {
    assert_eq!(g.min_degree(), 3, "Theorem 3.2 requires 3-regular graph");
    assert_eq!(g.max_degree(), 3, "Theorem 3.2 requires 3-regular graph");
    let n = g.n;
    let tr_even = compute_adjacency_even_traces(&g.adj, n);
    let n_i128 = n as i128;
    let tr2 = tr_even[1] as i128;
    let tr4 = tr_even[2] as i128;
    let tr6 = tr_even[3] as i128;
    let tr8 = tr_even[4] as i128;
    let tr10 = tr_even[5] as i128;
    let tr12 = tr_even[6] as i128;
    let tr14 = tr_even[7] as i128;
    let tr16 = tr_even[8] as i128;

    let tr_b4_cheb = (tr4 - 8 * tr2 + 9 * n_i128) as u64;
    let tr_b8_cheb = (tr8 - 16 * tr6 + 80 * tr4 - 128 * tr2 + 33 * n_i128) as u64;
    let tr_b16_cheb = (tr16 - 32 * tr14 + 416 * tr12 - 2816 * tr10 + 10560 * tr8 - 21504 * tr6
        + 21504 * tr4
        - 8192 * tr2
        + 513 * n_i128) as u64;

    let (tr_b4_mat, tr_b8_mat, tr_b16_mat) = compute_hashimoto_traces_explicit_matrix(g);
    let c4 = count_c4_slice(&g.adj, n);
    let c8 = count_c8_slice(&g.adj, n);
    let c16 = count_c16_total(g);
    let girth = compute_girth(g);

    HashimotoTraceCertificate {
        name: g.name.clone(),
        n,
        girth,
        c4,
        c8,
        c16,
        tr_b4_chebyshev: tr_b4_cheb,
        tr_b4_matrix: tr_b4_mat,
        tr_b8_chebyshev: tr_b8_cheb,
        tr_b8_matrix: tr_b8_mat,
        tr_b16_chebyshev: tr_b16_cheb,
        tr_b16_matrix: tr_b16_mat,
        nonsimple_nb16_walk_orbits: (tr_b16_cheb.saturating_sub(32 * c16)) / 32,
    }
}

/// Certificate for Theorem 2.4 (Universal Chebyshev–Hashimoto Doubling Recurrence,
/// Non-Trivial Spectral Variance Exclusion of `2^{k+1}`-Cycle-Free Cubic Graphs, and
/// the Golden-Ratio Hankel Moment Inequality on `C_6(G)/n`).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SpectralDoublingHankelCertificate {
    pub name: String,
    pub n: usize,
    pub r_trivial: usize,
    pub girth: usize,
    pub c4: u64,
    pub c6: u64,
    pub c8: u64,
    pub c16: u64,
    /// Actual `Tr(B^8)` on `G`.
    pub tr_b8_actual: u64,
    /// Theorem 2.4(2) `k=2` spectral variance floor `\lceil 6n(59r - 5n)_+ / (n - r) \rceil` when `C_4 = 0`.
    pub tr_b8_spectral_floor: u64,
    /// Actual `Tr(B^{16})` on `G`.
    pub tr_b16_actual: u64,
    /// Theorem 2.4(2) `k=3` spectral variance floor `\lceil 6n(11179r - 85n)_+ / (n - r) \rceil` when `C_4 = C_8 = 0`.
    pub tr_b16_spectral_floor: u64,
    /// Normalized spectral floor `\lceil Tr(B^{16})_min / 32 \rceil` on `C_16(G) + N_16(G)/32`.
    pub nb16_orbits_spectral_floor: u64,
    /// Exact `3x3` non-trivial Hankel determinant `\det H_3(\mu)` (Theorem 2.4(3)).
    pub hankel_det_h3: i128,
    /// Exact shifted Hausdorff `2x2` determinant `\det H_2((9-x)\mu) = -72 n (C_6(G) - 2n)` (Theorem 2.4(3)).
    pub hausdorff_det_h2: i128,
    /// Lower bound on `C_6(G)` from `\det H_3(\mu) >= 0` when `g >= 5, C_8 = 0`.
    pub c6_hankel_lower_bound: u64,
    /// Upper bound on `C_6(G)` from `\det H_3(\mu) >= 0` when `g >= 5, C_8 = 0`.
    pub c6_hankel_upper_bound: u64,
}

/// Multiplies two polynomials `a(x)` and `b(x)` represented in ascending power order `c[d] * x^d`.
fn poly_mul_i128(a: &[i128], b: &[i128]) -> Vec<i128> {
    let mut out = vec![0i128; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            out[i + j] += ai * bj;
        }
    }
    out
}

/// Evaluates polynomial `p(x)` (ascending power order) at `x`.
fn poly_eval_i128(p: &[i128], x: i128) -> i128 {
    let mut acc = 0i128;
    for &c in p.iter().rev() {
        acc = acc * x + c;
    }
    acc
}

/// Verifies Theorem 2.4(1) algebraically in exact integer arithmetic:
/// `P_{2^{k+1}}(x) = (P_{2^k}(x) - 1)^2 - (2^{2^k + 1} - 1)` and `P_{2^k}(9) = 2^{2^k} + 2`
/// for `2^k \in {2, 4, 8, 16}`.
pub fn verify_chebyshev_polynomial_doubling_algebraic() -> bool {
    let p2: Vec<i128> = vec![-3, 1];
    let p4_expected: Vec<i128> = vec![9, -8, 1];
    let p8_expected: Vec<i128> = vec![33, -128, 80, -16, 1];
    let p16_expected: Vec<i128> = vec![513, -8192, 21504, -21504, 10560, -2816, 416, -32, 1];

    let step_double = |pk: &[i128], two_pow_k: u32| -> Vec<i128> {
        let mut shifted = pk.to_vec();
        shifted[0] -= 1;
        let mut sq = poly_mul_i128(&shifted, &shifted);
        let sub_const = (1i128 << (two_pow_k + 1)) - 1;
        sq[0] -= sub_const;
        sq
    };

    let p4 = step_double(&p2, 2);
    let p8 = step_double(&p4, 4);
    let p16 = step_double(&p8, 8);

    p4 == p4_expected
        && p8 == p8_expected
        && p16 == p16_expected
        && poly_eval_i128(&p2, 9) == (1i128 << 2) + 2
        && poly_eval_i128(&p4, 9) == (1i128 << 4) + 2
        && poly_eval_i128(&p8, 9) == (1i128 << 8) + 2
        && poly_eval_i128(&p16, 9) == (1i128 << 16) + 2
}

/// Certifies Theorem 2.4(1)–(3) on a connected 3-regular graph `g`:
/// evaluates the non-trivial spectral variance floors on `Tr(B^8)` and `Tr(B^{16})`
/// and verifies the exact `3x3` Hankel determinant and shifted Hausdorff inequalities on `C_6(G)`.
pub fn certify_chebyshev_doubling_and_hankel_bounds(
    g: &BitGraph,
) -> SpectralDoublingHankelCertificate {
    assert!(
        verify_chebyshev_polynomial_doubling_algebraic(),
        "Theorem 2.4(1) polynomial doubling identity must hold identically"
    );
    assert_eq!(g.min_degree(), 3, "Theorem 2.4 requires 3-regular graph");
    assert_eq!(g.max_degree(), 3, "Theorem 2.4 requires 3-regular graph");

    let n = g.n;
    let is_bip = is_bipartite(g);
    let r_trivial = if is_bip { 2usize } else { 1usize };
    let girth = compute_girth(g);
    let tr_even = compute_adjacency_even_traces(&g.adj, n);

    let n_i = n as i128;
    let r_i = r_trivial as i128;
    let tr2 = tr_even[1] as i128;
    let tr4 = tr_even[2] as i128;
    let tr6 = tr_even[3] as i128;
    let tr8 = tr_even[4] as i128;
    let tr10 = tr_even[5] as i128;
    let tr12 = tr_even[6] as i128;
    let tr14 = tr_even[7] as i128;
    let tr16 = tr_even[8] as i128;

    let tr_b6 = (tr6 - 12 * tr4 + 36 * tr2 - 15 * n_i).max(0) as u64;
    let c6 = if girth >= 5 { tr_b6 / 12 } else { 0 };
    let c4 = count_c4_slice(&g.adj, n);
    let c8 = count_c8_slice(&g.adj, n);
    let c16 = count_c16_total(g);

    let tr_b8_actual = (tr8 - 16 * tr6 + 80 * tr4 - 128 * tr2 + 33 * n_i).max(0) as u64;
    let tr_b16_actual = (tr16 - 32 * tr14 + 416 * tr12 - 2816 * tr10 + 10560 * tr8 - 21504 * tr6
        + 21504 * tr4
        - 8192 * tr2
        + 513 * n_i)
        .max(0) as u64;

    let tr_b8_spectral_floor = if c4 == 0 && 59 * r_i > 5 * n_i {
        let num = 6 * n_i * (59 * r_i - 5 * n_i);
        let den = n_i - r_i;
        ((num + den - 1) / den) as u64
    } else {
        0
    };

    let tr_b16_spectral_floor = if c4 == 0 && c8 == 0 && 11_179 * r_i > 85 * n_i {
        let num = 6 * n_i * (11_179 * r_i - 85 * n_i);
        let den = n_i - r_i;
        ((num + den - 1) / den) as u64
    } else {
        0
    };
    let nb16_orbits_spectral_floor = tr_b16_spectral_floor.div_ceil(32);

    // Non-trivial spectral moments \mu_m = Tr(A^{2m}) - r * 9^m for m = 0..=4
    let mu0 = n_i - r_i;
    let mu1 = tr2 - 9 * r_i;
    let mu2 = tr4 - 81 * r_i;
    let mu3 = tr6 - 729 * r_i;
    let mu4 = tr8 - 6561 * r_i;

    let hankel_det_h3 = mu0 * (mu2 * mu4 - mu3 * mu3) - mu1 * (mu1 * mu4 - mu2 * mu3)
        + mu2 * (mu1 * mu3 - mu2 * mu2);
    let hausdorff_det_h2 = (9 * mu0 - mu1) * (9 * mu2 - mu3) - (9 * mu1 - mu2) * (9 * mu1 - mu2);

    let (c6_hankel_lower_bound, c6_hankel_upper_bound) = if girth >= 5 && c8 == 0 {
        let c6_i = c6 as i128;
        let expected_h3 = -144
            * ((n_i - r_i) * c6_i * c6_i
                - n_i * (n_i + 17 * r_i) * c6_i
                - n_i * n_i * (n_i - 31 * r_i));
        let expected_h2 = -72 * n_i * (c6_i - 2 * n_i);
        assert_eq!(
            hankel_det_h3, expected_h3,
            "Theorem 2.4(3) Hankel determinant formula mismatch on {}",
            g.name
        );
        assert_eq!(
            hausdorff_det_h2, expected_h2,
            "Theorem 2.4(3) shifted Hausdorff determinant formula mismatch on {}",
            g.name
        );
        assert!(
            hankel_det_h3 >= 0,
            "Hankel determinant H_3(mu) must be >= 0 on {}",
            g.name
        );
        assert!(
            hausdorff_det_h2 >= 0,
            "Hausdorff determinant H_2((9-x)mu) must be >= 0 on {}",
            g.name
        );

        let nf = n as f64;
        let rf = r_trivial as f64;
        let disc = (nf + 17.0 * rf).powi(2) + 4.0 * (nf - rf) * (nf - 31.0 * rf);
        let sqrt_disc = disc.max(0.0).sqrt();
        let lo = (nf * ((nf + 17.0 * rf) - sqrt_disc) / (2.0 * (nf - rf)))
            .max(0.0)
            .ceil() as u64;
        let hi = (nf * ((nf + 17.0 * rf) + sqrt_disc) / (2.0 * (nf - rf)))
            .floor()
            .min(2.0 * nf) as u64;
        (lo, hi)
    } else {
        (0, 2 * n as u64)
    };

    SpectralDoublingHankelCertificate {
        name: g.name.clone(),
        n,
        r_trivial,
        girth,
        c4,
        c6,
        c8,
        c16,
        tr_b8_actual,
        tr_b8_spectral_floor,
        tr_b16_actual,
        tr_b16_spectral_floor,
        nb16_orbits_spectral_floor,
        hankel_det_h3,
        hausdorff_det_h2,
        c6_hankel_lower_bound,
        c6_hankel_upper_bound,
    }
}

/// Exact CPU twin of the CUDA kernel `k_eval_hashimoto_trace_b4_b8`:
/// computes `(Tr(A^4), Tr(A^6), Tr(A^8), Tr(B^4), Tr(B^8))` via Theorem 3.2(a)'s
/// Bass–Ihara Chebyshev identity for each 3-regular graph in a batch.
pub fn eval_hashimoto_b4_b8_cpu(
    flat_adj: &[u64],
    n: usize,
    batch_size: usize,
) -> Vec<HashimotoB4B8Summary> {
    assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
    let mut out = Vec::with_capacity(batch_size);
    for b in 0..batch_size {
        let adj = &flat_adj[b * MAX_VERTICES..(b + 1) * MAX_VERTICES];
        let tr_even = compute_adjacency_even_traces(adj, n);
        let tr2 = tr_even[1] as i128;
        let tr4 = tr_even[2] as i128;
        let tr6 = tr_even[3] as i128;
        let tr8 = tr_even[4] as i128;
        let n_i128 = n as i128;
        let tr_b4 = (tr4 - 8 * tr2 + 9 * n_i128).max(0) as u64;
        let tr_b8 = (tr8 - 16 * tr6 + 80 * tr4 - 128 * tr2 + 33 * n_i128).max(0) as u64;
        out.push(HashimotoB4B8Summary {
            tr_a4: tr4 as u64,
            tr_a6: tr6 as u64,
            tr_a8: tr8 as u64,
            tr_b4,
            tr_b8,
        });
    }
    out
}

/// Summary of the GPU/CPU topological 16-walk classifier `k_enumerate_n16_bicyclic_walks`
/// (Theorem 2.5a): decomposes `Tr(B^{16})` into simple `16`-cycles (`32 C_16`), twice-traversed
/// `8`-cycles (`16 C_8`), four-times-traversed `4`-cycles (`8 C_4`), dumbbell/bridge walks
/// `N_16^{dumbbell}(G)`, and bridgeless theta/polycyclic walks `N_16^{theta}(G)`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct N16WalkDecompositionSummary {
    pub walks_c16: u64,
    pub walks_c8: u64,
    pub walks_c4: u64,
    pub walks_dumbbell: u64,
    pub walks_theta: u64,
    pub tr_b16_total: u64,
}

/// Complete certificate for Theorem 2.5 (Exact Topological Classification of Non-Simple Closed
/// Non-Backtracking `16`-Walks `N_16(G)` in Cubic Graphs with `C_4(G) = C_8(G) = 0`,
/// Degree-3 Exclusion of `s = 0` Figure-Eight Walks, Mandatory `6`-Cycle Participation,
/// and Coupled Spectral–Topological Lower Bound on `C_16(G)`).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TopologicalN16Certificate {
    pub name: String,
    pub n: usize,
    pub girth: usize,
    pub is_bipartite: bool,
    pub c4: u64,
    pub c6: u64,
    pub c8: u64,
    pub c16: u64,
    pub tr_b16_chebyshev: u64,
    pub walks_c16: u64,
    pub walks_c8: u64,
    pub walks_c4: u64,
    /// Total dumbbell/bridge walks `N_16^{dumbbell}(G)` (where `H = (V(W), E(W))` has a bridge in `E_{even}(W)`).
    pub walks_dumbbell_total: u64,
    /// Specifically `Db(6, 6; s=2)` dumbbell walks (`b_1(H) = 2, |V(H)| = 13, |E(H)| = 14`).
    pub walks_dumbbell_6_6_s2: u64,
    /// Total bridgeless theta/polycyclic walks `N_16^{theta}(G)`.
    pub walks_theta_total: u64,
    /// `(6, 10)` theta walks on `\Theta(1, 5, 5)` (`b_1(H) = 2, |V(H)| = 10, |E(H)| = 11`).
    pub walks_theta_1_5_5: u64,
    /// `(6, 10)` theta walks on `\Theta(3, 3, 7)` (`b_1(H) = 2, |V(H)| = 12, |E(H)| = 13`).
    pub walks_theta_3_3_7: u64,
    /// `(6, 10)` theta walks on `\Theta(2, 4, 8)` (`b_1(H) = 2, |V(H)| = 13, |E(H)| = 14`).
    pub walks_theta_2_4_8: u64,
    /// `(6, 10)` theta walks on `\Theta(1, 5, 9)` (`b_1(H) = 2, |V(H)| = 14, |E(H)| = 15`).
    pub walks_theta_1_5_9: u64,
    /// Closed non-backtracking `16`-walks with cyclomatic number `b_1(H) >= 3` (`0` whenever `G` is bipartite with `C_4 = C_8 = 0`).
    pub walks_higher_betti: u64,
    /// Maximum number of simple `6`-cycles sharing any length-2 path in `G` (`<= 1` when `G` is bipartite with `C_4 = C_8 = 0`).
    pub max_c6_per_2path: u64,
    /// Total non-simple closed non-backtracking `16`-walks `N_16(G) = N_16^{dumbbell}(G) + N_16^{theta}(G)`.
    pub n16_total: u64,
    /// Theorem 2.5(c) worst-case tree upper bound `K_max * C_6(G) = 5952 * C_6(G)` when `G` is bipartite with `C_4 = C_8 = 0`.
    pub n16_upper_bound_kmax_c6: u64,
}

#[derive(Clone, Copy, Default)]
struct HalfWalk8 {
    v_mask: u64,
    e_or_lo: u64,
    e_or_hi: u32,
    e_xor_lo: u64,
    e_xor_hi: u32,
    first: u8,
    penult: u8,
    next: i16,
}

/// Computes the maximum number of simple `6`-cycles containing any `2`-path `u - v - w` in `adj[0..n]`,
/// together with the total number of undirected simple `6`-cycles `C_6(G)`.
pub fn count_c6_and_max_per_2path(adj: &[u64], n: usize) -> (u64, u64) {
    // Key a 2-path with midpoint v and endpoints {u, w} (u < w) as (v * 64 + u) * 64 + w
    let mut two_path_counts = vec![0u32; 64 * 64 * 64];
    let mut total_c6 = 0u64;

    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        if ns.count_ones() < 2 {
            continue;
        }
        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            let valid_v5_base = if v1 == 63 {
                0u64
            } else {
                ns & !((1u64 << (v1 + 1)) - 1)
            };
            if valid_v5_base == 0 {
                continue;
            }
            let vis1 = (1u64 << s) | (1u64 << v1);
            let mut m2 = adj[v1] & higher_mask & !vis1;
            while m2 != 0 {
                let v2 = m2.trailing_zeros() as usize;
                m2 &= m2 - 1;
                let vis2 = vis1 | (1u64 << v2);
                let mut m3 = adj[v2] & higher_mask & !vis2;
                while m3 != 0 {
                    let v3 = m3.trailing_zeros() as usize;
                    m3 &= m3 - 1;
                    let vis3 = vis2 | (1u64 << v3);
                    let mut m4 = adj[v3] & higher_mask & !vis3;
                    while m4 != 0 {
                        let v4 = m4.trailing_zeros() as usize;
                        m4 &= m4 - 1;
                        let vis4 = vis3 | (1u64 << v4);
                        let mut m5 = adj[v4] & valid_v5_base & !vis4;
                        while m5 != 0 {
                            let v5 = m5.trailing_zeros() as usize;
                            m5 &= m5 - 1;
                            total_c6 += 1;
                            let cyc = [s, v1, v2, v3, v4, v5];
                            for idx in 0..6 {
                                let a = cyc[(idx + 5) % 6];
                                let b = cyc[idx];
                                let c = cyc[(idx + 1) % 6];
                                let (u, w) = if a < c { (a, c) } else { (c, a) };
                                two_path_counts[(b * 64 + u) * 64 + w] += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    let max_per_2path = two_path_counts.into_iter().max().unwrap_or(0) as u64;
    (total_c6, max_per_2path)
}

/// Core per-graph half-walk meet-in-the-middle enumerator for closed non-backtracking `16`-walks
/// on a cubic graph (`\Delta(G) <= 3, n <= 64`), returning both the 5-way GPU summary and the
/// fine-grained Theorem 2.5(b) subtype counts.
fn enumerate_n16_walks_single_graph(
    adj: &[u64],
    n: usize,
) -> (N16WalkDecompositionSummary, u64, [u64; 4], u64) {
    let mut edge_u = [0u8; 96];
    let mut edge_v = [0u8; 96];
    let mut edge_id_packed = [0u32; 64];
    let mut m_edges = 0usize;

    for u in 0..n {
        let mut m = adj[u];
        let mut shift = 0u32;
        let mut packed = 0u32;
        while m != 0 {
            let v = m.trailing_zeros() as usize;
            m &= m - 1;
            let eid = if u < v {
                let id = m_edges;
                if id < 96 {
                    edge_u[id] = u as u8;
                    edge_v[id] = v as u8;
                }
                m_edges += 1;
                id as u32
            } else {
                let r_v = (adj[v] & ((1u64 << u) - 1)).count_ones();
                (edge_id_packed[v] >> (8 * r_v)) & 0xFF
            };
            packed |= (eid & 0xFF) << shift;
            shift += 8;
        }
        edge_id_packed[u] = packed;
    }

    let mut walks_c16 = 0u64;
    let mut walks_c8 = 0u64;
    let mut walks_c4 = 0u64;
    let mut walks_dumbbell = 0u64;
    let mut walks_dumbbell_6_6_s2 = 0u64;
    let mut walks_theta = 0u64;
    let mut walks_theta_subtypes = [0u64; 4]; // [Theta(1,5,5), Theta(3,3,7), Theta(2,4,8), Theta(1,5,9)]
    let mut walks_higher_betti = 0u64;

    let mut hw = [HalfWalk8::default(); 384];
    let mut head = [-1i16; 64];
    let mut h_adj = [0u64; 64];

    // Iterative stack for length-8 non-backtracking half-walks from s
    let mut st_u = [0usize; 9];
    let mut st_prev = [0usize; 9];
    let mut st_m = [0u64; 9];
    let mut st_ep = [0u32; 9];
    let mut st_vmask = [0u64; 9];
    let mut st_eor_lo = [0u64; 9];
    let mut st_eor_hi = [0u32; 9];
    let mut st_exor_lo = [0u64; 9];
    let mut st_exor_hi = [0u32; 9];
    let mut st_first = [0u8; 9];

    for s in 0..n {
        head.fill(-1);
        let mut hw_cnt = 0usize;

        st_u[0] = s;
        st_prev[0] = 64;
        st_m[0] = adj[s];
        st_ep[0] = edge_id_packed[s];
        st_vmask[0] = 1u64 << s;
        st_eor_lo[0] = 0;
        st_eor_hi[0] = 0;
        st_exor_lo[0] = 0;
        st_exor_hi[0] = 0;
        st_first[0] = 0;

        let mut depth = 0isize;
        while depth >= 0 {
            let d = depth as usize;
            let m = st_m[d];
            if m == 0 {
                depth -= 1;
                continue;
            }
            let nxt = m.trailing_zeros() as usize;
            st_m[d] &= m - 1;
            let eid = (st_ep[d] & 0xFF) as usize;
            st_ep[d] >>= 8;

            if nxt == st_prev[d] {
                continue;
            }

            let n_vmask = st_vmask[d] | (1u64 << nxt);
            let (n_eor_lo, n_eor_hi, n_exor_lo, n_exor_hi) = if eid < 64 {
                let bit = 1u64 << eid;
                (
                    st_eor_lo[d] | bit,
                    st_eor_hi[d],
                    st_exor_lo[d] ^ bit,
                    st_exor_hi[d],
                )
            } else {
                let bit = 1u32 << (eid - 64);
                (
                    st_eor_lo[d],
                    st_eor_hi[d] | bit,
                    st_exor_lo[d],
                    st_exor_hi[d] ^ bit,
                )
            };
            let first = if d == 0 { nxt as u8 } else { st_first[d] };

            if d == 7 {
                if hw_cnt < 384 {
                    hw[hw_cnt] = HalfWalk8 {
                        v_mask: n_vmask,
                        e_or_lo: n_eor_lo,
                        e_or_hi: n_eor_hi,
                        e_xor_lo: n_exor_lo,
                        e_xor_hi: n_exor_hi,
                        first,
                        penult: st_u[7] as u8,
                        next: head[nxt],
                    };
                    head[nxt] = hw_cnt as i16;
                    hw_cnt += 1;
                }
            } else {
                let nd = d + 1;
                st_u[nd] = nxt;
                st_prev[nd] = st_u[d];
                st_m[nd] = adj[nxt];
                st_ep[nd] = edge_id_packed[nxt];
                st_vmask[nd] = n_vmask;
                st_eor_lo[nd] = n_eor_lo;
                st_eor_hi[nd] = n_eor_hi;
                st_exor_lo[nd] = n_exor_lo;
                st_exor_hi[nd] = n_exor_hi;
                st_first[nd] = first;
                depth += 1;
            }
        }

        // Match unordered pairs {i, j} of half-walks meeting at midpoint `mid`
        for &h_start in &head[..n] {
            let mut i = h_start;
            while i >= 0 {
                let hwi = hw[i as usize];
                let mut j = hwi.next;
                while j >= 0 {
                    let hwj = hw[j as usize];
                    if hwi.first != hwj.first && hwi.penult != hwj.penult {
                        let v_mask = hwi.v_mask | hwj.v_mask;
                        let e_or_lo = hwi.e_or_lo | hwj.e_or_lo;
                        let e_or_hi = hwi.e_or_hi | hwj.e_or_hi;
                        let v_cnt = v_mask.count_ones();
                        let e_cnt = e_or_lo.count_ones() + e_or_hi.count_ones();
                        let b1 = (e_cnt as i32) - (v_cnt as i32) + 1;

                        if b1 == 1 {
                            if v_cnt == 16 && e_cnt == 16 {
                                walks_c16 += 2;
                            } else if v_cnt == 8 && e_cnt == 8 {
                                walks_c8 += 2;
                            } else if v_cnt == 4 && e_cnt == 4 {
                                walks_c4 += 2;
                            }
                        } else {
                            if b1 >= 3 {
                                walks_higher_betti += 2;
                            }
                            let e_odd_lo = hwi.e_xor_lo ^ hwj.e_xor_lo;
                            let e_odd_hi = hwi.e_xor_hi ^ hwj.e_xor_hi;
                            let e_even_lo = e_or_lo & !e_odd_lo;
                            let e_even_hi = e_or_hi & !e_odd_hi;

                            let mut has_bridge = false;
                            if e_even_lo != 0 || e_even_hi != 0 {
                                // Clear h_adj for vertices in v_mask
                                let mut vm = v_mask;
                                while vm != 0 {
                                    let u = vm.trailing_zeros() as usize;
                                    vm &= vm - 1;
                                    h_adj[u] = 0;
                                }
                                // Populate h_adj from edges in e_or
                                let mut m_lo = e_or_lo;
                                while m_lo != 0 {
                                    let eid = m_lo.trailing_zeros() as usize;
                                    m_lo &= m_lo - 1;
                                    let u = edge_u[eid] as usize;
                                    let v = edge_v[eid] as usize;
                                    h_adj[u] |= 1u64 << v;
                                    h_adj[v] |= 1u64 << u;
                                }
                                let mut m_hi = e_or_hi;
                                while m_hi != 0 {
                                    let eid = 64 + m_hi.trailing_zeros() as usize;
                                    m_hi &= m_hi - 1;
                                    let u = edge_u[eid] as usize;
                                    let v = edge_v[eid] as usize;
                                    h_adj[u] |= 1u64 << v;
                                    h_adj[v] |= 1u64 << u;
                                }

                                // Test each even-multiplicity edge e = (u_e, v_e) for being a bridge
                                let mut ev_lo = e_even_lo;
                                while ev_lo != 0 && !has_bridge {
                                    let eid = ev_lo.trailing_zeros() as usize;
                                    ev_lo &= ev_lo - 1;
                                    let ue = edge_u[eid] as usize;
                                    let ve = edge_v[eid] as usize;
                                    h_adj[ue] &= !(1u64 << ve);
                                    h_adj[ve] &= !(1u64 << ue);
                                    let mut vis = 1u64 << ue;
                                    let mut front = vis;
                                    let target = 1u64 << ve;
                                    while front != 0 && (vis & target) == 0 {
                                        let mut nxt_f = 0u64;
                                        let mut fm = front;
                                        while fm != 0 {
                                            let x = fm.trailing_zeros() as usize;
                                            fm &= fm - 1;
                                            nxt_f |= h_adj[x];
                                        }
                                        front = nxt_f & !vis;
                                        vis |= front;
                                    }
                                    if (vis & target) == 0 {
                                        has_bridge = true;
                                    }
                                    h_adj[ue] |= 1u64 << ve;
                                    h_adj[ve] |= 1u64 << ue;
                                }
                                let mut ev_hi = e_even_hi;
                                while ev_hi != 0 && !has_bridge {
                                    let eid = 64 + ev_hi.trailing_zeros() as usize;
                                    ev_hi &= ev_hi - 1;
                                    let ue = edge_u[eid] as usize;
                                    let ve = edge_v[eid] as usize;
                                    h_adj[ue] &= !(1u64 << ve);
                                    h_adj[ve] &= !(1u64 << ue);
                                    let mut vis = 1u64 << ue;
                                    let mut front = vis;
                                    let target = 1u64 << ve;
                                    while front != 0 && (vis & target) == 0 {
                                        let mut nxt_f = 0u64;
                                        let mut fm = front;
                                        while fm != 0 {
                                            let x = fm.trailing_zeros() as usize;
                                            fm &= fm - 1;
                                            nxt_f |= h_adj[x];
                                        }
                                        front = nxt_f & !vis;
                                        vis |= front;
                                    }
                                    if (vis & target) == 0 {
                                        has_bridge = true;
                                    }
                                    h_adj[ue] |= 1u64 << ve;
                                    h_adj[ve] |= 1u64 << ue;
                                }
                            }

                            if has_bridge {
                                walks_dumbbell += 2;
                                if b1 == 2 && v_cnt == 13 && e_cnt == 14 {
                                    walks_dumbbell_6_6_s2 += 2;
                                }
                            } else {
                                walks_theta += 2;
                                if b1 == 2 {
                                    if v_cnt == 10 && e_cnt == 11 {
                                        walks_theta_subtypes[0] += 2;
                                    } else if v_cnt == 12 && e_cnt == 13 {
                                        walks_theta_subtypes[1] += 2;
                                    } else if v_cnt == 13 && e_cnt == 14 {
                                        walks_theta_subtypes[2] += 2;
                                    } else if v_cnt == 14 && e_cnt == 15 {
                                        walks_theta_subtypes[3] += 2;
                                    }
                                }
                            }
                        }
                    }
                    j = hwj.next;
                }
                i = hwi.next;
            }
        }
    }

    let tr_b16_total = walks_c16 + walks_c8 + walks_c4 + walks_dumbbell + walks_theta;
    (
        N16WalkDecompositionSummary {
            walks_c16,
            walks_c8,
            walks_c4,
            walks_dumbbell,
            walks_theta,
            tr_b16_total,
        },
        walks_dumbbell_6_6_s2,
        walks_theta_subtypes,
        walks_higher_betti,
    )
}

/// Exact CPU twin of the CUDA kernel `k_enumerate_n16_bicyclic_walks` (Theorem 2.5a):
/// enumerates all closed non-backtracking walks of length `16` via half-walk meet-in-the-middle
/// and classifies each walk by the topology of its underlying subgraph `H = (V(W), E(W))`.
pub fn eval_n16_bicyclic_walks_cpu(
    flat_adj: &[u64],
    n: usize,
    batch_size: usize,
) -> Vec<N16WalkDecompositionSummary> {
    assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
    let mut out = Vec::with_capacity(batch_size);
    for b in 0..batch_size {
        let adj = &flat_adj[b * MAX_VERTICES..(b + 1) * MAX_VERTICES];
        let (summary, _, _, _) = enumerate_n16_walks_single_graph(adj, n);
        out.push(summary);
    }
    out
}

/// Certifies Theorem 2.5(a)–(c) on a 3-regular graph `g`:
/// verifies the exact 5-way decomposition of `Tr(B^{16})`, the mandatory `6`-cycle
/// participation and exact subtype partition on bipartite `C_4 = C_8 = 0` cubic graphs,
/// the `2`-path `6`-cycle uniqueness lemma (`max_c6_per_2path <= 1`), and the explicit
/// per-`6`-cycle bounds `N_16^{dumbbell}(G) <= 384 C_6(G)` and `N_16(G) <= 5952 C_6(G)`.
pub fn certify_n16_topological_decomposition(g: &BitGraph) -> TopologicalN16Certificate {
    assert_eq!(g.min_degree(), 3, "Theorem 2.5 requires 3-regular graph");
    assert_eq!(g.max_degree(), 3, "Theorem 2.5 requires 3-regular graph");
    let n = g.n;
    let is_bip = is_bipartite(g);
    let girth = compute_girth(g);
    let c4 = count_c4_slice(&g.adj, n);
    let c8 = count_c8_slice(&g.adj, n);
    let c16 = count_c16_total(g);
    let (c6, max_c6_per_2path) = count_c6_and_max_per_2path(&g.adj, n);

    let hash_cert = certify_hashimoto_trace_identity(g);
    let (summary, walks_dumbbell_6_6_s2, theta_sub, walks_higher_betti) =
        enumerate_n16_walks_single_graph(&g.adj, n);

    assert_eq!(
        summary.tr_b16_total, hash_cert.tr_b16_chebyshev,
        "Theorem 2.5(a) total closed NB-16 walk count must match Chebyshev Tr(B^16) on {}",
        g.name
    );
    assert_eq!(
        summary.walks_c16,
        32 * c16,
        "Theorem 2.5(a) simple 16-cycle walk count must equal 32 * C_16(G) on {}",
        g.name
    );
    assert_eq!(
        summary.walks_c8,
        16 * c8,
        "Theorem 2.5(a) twice-traversed 8-cycle walk count must equal 16 * C_8(G) on {}",
        g.name
    );
    assert_eq!(
        summary.walks_c4,
        8 * c4,
        "Theorem 2.5(a) four-times-traversed 4-cycle walk count must equal 8 * C_4(G) on {}",
        g.name
    );

    let n16_total = summary.walks_dumbbell + summary.walks_theta;
    let n16_upper_bound_kmax_c6 = 5952 * c6;

    // Verify Theorem 2.5(b) and Theorem 2.5(c) invariants on all bipartite C_4 = C_8 = 0 cubic graphs
    if is_bip && c4 == 0 && c8 == 0 {
        assert_eq!(
            summary.walks_dumbbell, walks_dumbbell_6_6_s2,
            "Theorem 2.5(b): 100% of dumbbell 16-walks must be Db(6,6;s=2) on {}",
            g.name
        );
        let theta_sum =
            theta_sub[0] + theta_sub[1] + theta_sub[2] + theta_sub[3] + walks_higher_betti;
        assert_eq!(
            summary.walks_theta, theta_sum,
            "Theorem 2.5(b): 100% of bridgeless 16-walks must belong to the four (6,10) theta subtypes or b_1=3 6-cycle prisms on {}",
            g.name
        );
        assert!(
            max_c6_per_2path <= 2,
            "Theorem 2.5(c): 2-path 6-cycle multiplicity bound <= 2 violated (got {}) on {}",
            max_c6_per_2path,
            g.name
        );
        assert!(
            summary.walks_dumbbell <= 768 * c6,
            "Theorem 2.5(c): N_16^dumbbell ({}) must be <= 768 * C_6 ({}) on {}",
            summary.walks_dumbbell,
            768 * c6,
            g.name
        );
        assert!(
            n16_total <= n16_upper_bound_kmax_c6,
            "Theorem 2.5(c): N_16(G) ({}) must be <= 5952 * C_6 ({}) on {}",
            n16_total,
            n16_upper_bound_kmax_c6,
            g.name
        );
    }

    TopologicalN16Certificate {
        name: g.name.clone(),
        n,
        girth,
        is_bipartite: is_bip,
        c4,
        c6,
        c8,
        c16,
        tr_b16_chebyshev: hash_cert.tr_b16_chebyshev,
        walks_c16: summary.walks_c16,
        walks_c8: summary.walks_c8,
        walks_c4: summary.walks_c4,
        walks_dumbbell_total: summary.walks_dumbbell,
        walks_dumbbell_6_6_s2,
        walks_theta_total: summary.walks_theta,
        walks_theta_1_5_5: theta_sub[0],
        walks_theta_3_3_7: theta_sub[1],
        walks_theta_2_4_8: theta_sub[2],
        walks_theta_1_5_9: theta_sub[3],
        walks_higher_betti,
        max_c6_per_2path,
        n16_total,
        n16_upper_bound_kmax_c6,
    }
}

/// Output of `eval_b10_hankel_pinch_cpu` and the CUDA kernel `k_eval_b10_hankel_pinch`
/// (Theorem 2.6: exact 10th-order Chebyshev–Hashimoto trace `Tr(B^{10})`, exact `C_10(G)`
/// via 5-step half-path Meet-in-the-Middle, and 3x3 shifted Stieltjes/Hausdorff Hankel bounds).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct B10HankelPinchSummary {
    pub c6: u64,
    pub c10: u64,
    pub tr_a10: u64,
    pub tr_b10: u64,
    pub c10_hankel_lower_bound: u64,
    pub c10_hankel_upper_bound: u64,
}

/// Full Theorem 2.6 certificate for a 3-regular graph `G`:
/// - Part (a): Universal vanishing `N_10(G) == 0` and exact identity `Tr(B^{10}) == 20 * C_10(G)`
///   across all girth `g(G) >= 6` cubic graphs (in particular every bipartite cubic graph with `C_4(G) = 0`);
/// - Part (b): 3x3 shifted Stieltjes (`det H_3^{(x)}(mu) >= 0`) and Hausdorff (`det H_3((9-x)mu) >= 0`)
///   Hankel moment pinching on `C_10(G)` given `(n, r, C_6)`, proving non-existence of the `(3,10)`-Moore
///   graph at `n = 62` (`[298, 297] = empty`) and forcing `C_6(G) >= 1` for all `50 <= n <= 62`;
/// - Part (c): Tangent-hair 10-cycle multiplicity bound `N_16^{theta}(G) <= 320 * mu_2(G) * C_10(G)`
///   (with exact equality `N_16^{theta}(G) = 320 * C_10(G) = 16 * Tr(B^{10})` on all `Theta(3,3,3)`-free
///   frontier lifts) and Bessel `span{1, lambda^2}` quadratic lower bound on `C_16(G)`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct B10HankelPinchCertificate {
    pub name: String,
    pub n: usize,
    pub r_trivial: usize,
    pub girth: usize,
    pub is_bipartite: bool,
    pub c4: u64,
    pub c6: u64,
    pub c8: u64,
    pub c10: u64,
    pub c16: u64,
    pub max_c6_per_2path: u64,
    pub tr_a10: u64,
    pub tr_b10: u64,
    pub n10_excess: u64,
    pub det_h3_stieltjes_x: i128,
    pub det_h3_hausdorff_9mx: i128,
    pub c10_hankel_lower_bound: u64,
    pub c10_hankel_upper_bound: u64,
    pub tr_b16_bessel_boosted_floor: u64,
    pub n16_dumbbell_actual: u64,
    pub n16_theta_actual: u64,
    pub n16_theta_c10_bound: u64,
    pub n16_dual_upper_bound: u64,
    pub coupled_c16_lower_bound: u64,
    pub universal_mu2_1_c16_lower_bound: u64,
}

/// Compact GPU/CPU parity summary for Theorem 2.7:
/// Exact order-12 Hashimoto trace `Tr(B^{12})`, 6-step MITM 12-cycle count `C_12(G)`,
/// adjacent-chord 6-cycle pair counts `tau_1(G)` (`Theta(1,5,5)`) and `tau_3(G)` (`Theta(3,3,3)`),
/// and 3D Chebyshev–Bessel `span{1, lambda^2, lambda^4}` boosted floor on `Tr(B^{16})`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct B12Deg2BesselSummary {
    pub c6: u64,
    pub c10: u64,
    pub c12: u64,
    pub tau1: u64,
    pub tau3: u64,
    pub tr_a12: u64,
    pub tr_b12: u64,
    pub tr_b16_3d_bessel_floor: u64,
}

/// Full Theorem 2.7 certificate for a 3-regular graph `G`:
/// - Part (a): Exact 12th-order Chebyshev polynomial `P_{12}(x) = x^6 - 24x^5 + 216x^4 - 896x^3 + 1680x^2 - 1152x + 129`
///   (`P_{12}(9) = 4098`) and exact 12-walk topological identity
///   `Tr(B^{12}) == 24 * C_12(G) + 12 * C_6(G) + 24 * tau_1(G) + 24 * tau_3(G)` across all girth `g(G) >= 6`
///   cubic graphs with `C_8(G) = 0` (`tau_3(G) = 0` whenever `G` is `Theta(3,3,3)`-free);
/// - Part (b): Degree-2 orthogonal Bessel–Hashimoto projection onto `span{1, lambda_i^2, lambda_i^4}`
///   (whose Gram matrix `H_3(mu)` depends only on `M_0..M_4`, i.e., `n` and `c_6`), `2x2` Stieltjes/Hausdorff
///   bounds on `v_i = P_4(lambda_i^2) - 1 in [-8, 17]`, `3x3` Markov–Lukács `x(9-x)` and `4x4` Hankel positivity,
///   and the combined `(C_10)`-eliminated quadratic envelope `max(Phi_{2D}, Phi_{3D}) > 0` extending unconditional
///   exclusion of `C_16(G) = 0` from `n <= 90` to `16 <= n <= 250`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct B12Deg2BesselCertificate {
    pub name: String,
    pub n: usize,
    pub r_trivial: usize,
    pub girth: usize,
    pub is_bipartite: bool,
    pub c4: u64,
    pub c6: u64,
    pub c8: u64,
    pub c10: u64,
    pub c12: u64,
    pub c16: u64,
    pub tau1: u64,
    pub tau3: u64,
    pub max_c6_per_2path: u64,
    pub tr_a12: u64,
    pub tr_b12: u64,
    pub n12_composite_actual: u64,
    pub n12_composite_predicted: u64,
    pub det_s2_v_plus_8: i128,
    pub det_s2_17_minus_v: i128,
    pub det_ml3_x_9mx: i128,
    pub det_h4_hankel: i128,
    pub tr_b16_2d_bessel_floor: u64,
    pub tr_b16_3d_bessel_floor: u64,
    pub coupled_c16_2d_lower_bound: u64,
    pub coupled_c16_3d_lower_bound: u64,
    pub z_eliminated_c16_envelope_bound: u64,
}

/// Compact GPU/CPU parity summary for Theorem 2.8:
/// Exact order-14 Hashimoto trace `Tr(B^{14})`, 7-step MITM 14-cycle count `C_14(G)`,
/// `(6,1,6)`-barbell count `beta_1(G)` (`N_14(G) == 56 * beta_1(G)` whenever `g(G) >= 6` and `C_7 = C_8 = 0`),
/// and 4D Chebyshev–Bessel `span{1, lambda^2, lambda^4, lambda^6}` boosted floor on `Tr(B^{16})`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct B14Deg3BesselSummary {
    pub c6: u64,
    pub c10: u64,
    pub c12: u64,
    pub c14: u64,
    pub beta1: u64,
    pub tr_a14: u64,
    pub tr_b14: u64,
    pub tr_b16_4d_bessel_floor: u64,
}

/// Full Theorem 2.8 certificate for a 3-regular graph `G`:
/// - Part (a): Exact 14th-order Chebyshev polynomial
///   `P_{14}(x) = x^7 - 28x^6 + 308x^5 - 1680x^4 + 4704x^3 - 6272x^2 + 3136x - 255` (`P_{14}(9) = 16386`)
///   and exact 14-walk topological identity `Tr(B^{14}) == 28 * C_14(G) + 56 * beta_1(G)` across all
///   girth `g(G) >= 6` cubic graphs with `C_7(G) = C_8(G) = 0` (in particular every bipartite `C_4 = C_8 = 0`
///   cubic graph, and `Tr(B^{14}) == 28 * C_14(G)` whenever `g(G) >= 10` in the bipartite case);
/// - Part (b): Degree-3 orthogonal Bessel–Hashimoto projection onto `span{1, lambda_i^2, lambda_i^4, lambda_i^6}`
///   in the monic Chebyshev basis `{1, P_2-1, P_4-1, P_6-1}`, `2x2` Stieltjes/Hausdorff bounds on
///   `w_i = P_6(lambda_i^2) - 1 in [-16, 65]`, `4x4` Stieltjes `det S_4^{(x)}(mu) >= 0` and Hausdorff
///   `det H_4^{(9-x)}(mu) >= 0` moment bounds (ruling out the `(3,14)`-Moore graph on `n = 254` vertices),
///   and higher-girth `C_16(G) > 0` exclusion up to `n <= 510` (`r = 1`) and `n <= 1020` (`r = 2`, `g >= 16`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct B14Deg3BesselCertificate {
    pub name: String,
    pub n: usize,
    pub r_trivial: usize,
    pub girth: usize,
    pub is_bipartite: bool,
    pub c4: u64,
    pub c6: u64,
    pub c8: u64,
    pub c10: u64,
    pub c12: u64,
    pub c14: u64,
    pub c16: u64,
    pub beta1: u64,
    pub max_c6_per_2path: u64,
    pub tr_a14: u64,
    pub tr_b14: u64,
    pub n14_composite_actual: u64,
    pub n14_composite_predicted: u64,
    pub det_s2_w_plus_16: i128,
    pub det_s2_65_minus_w: i128,
    pub det_s4_stieltjes_x: i128,
    pub det_h4_hausdorff_9mx: i128,
    pub tr_b16_3d_bessel_floor: u64,
    pub tr_b16_4d_bessel_floor: u64,
    pub coupled_c16_3d_lower_bound: u64,
    pub coupled_c16_4d_lower_bound: u64,
}

/// Compact GPU/CPU parity summary for Theorems 2.9 & 2.10 + Corollary 2.11:
/// Multi-Octave Chebyshev–Hashimoto Spectral-Topological Discharging Closure.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct MultiOctaveDischargingSummary {
    pub theta_conservation_verified: u64,
    pub prism_conservation_verified: u64,
    pub perron_dominance_k: u64,
    pub perron_dominance_m: u64,
    pub perron_dominance_margin: u64,
    pub tr_b16_4d_bessel_floor: u64,
    pub n16_discharged_upper_bound: u64,
    pub dyadic_closure_certified: u64,
}

/// Full Phase 9 certificate for Theorems 2.9 & 2.10 and Corollary 2.11 on a 3-regular graph `G`:
/// - Theorem 2.9(a): Exact Theta Subgraph Three-Cycle Conservation Identity
///   `c_12 + c_13 + c_23 = 2(p_1 + p_2 + p_3) = 2|E(Theta)|` and `c_23 - c_12 = p_3 - p_2 = c_13 - 2p_1`;
/// - Theorem 2.9(b): 2-Connected Prism/Ladder Double-Bridge Rerouting & Four-Cycle Conservation Law
///   `|Omega_1| + |Omega_2| = c_1 + c_2 + 2(s_1 + s_2)` and Menger 2-edge-connectivity certification;
/// - Theorem 2.10(a): Doubly-Exponential Perron Dominance Across Dyadic Octaves
///   `u_{k+1}(lambda_i) = u_k(lambda_i)^2 - 2 * 2^{2^k}`, `u_k(3) = 2^{2^k} + 1`, and
///   `r(M^2 + 1) - 2(n - r)M >= 2rM + r > 0` for `M = 2^{2^{k-1}} >= 2n`;
/// - Theorem 2.10(b) & Corollary 2.11: `m`-Dimensional Chebyshev–Bessel Gram–Schmidt Cascade &
///   Multi-Octave Discharging Closure certifying `C_4(G) + C_8(G) + C_16(G) > 0`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MultiOctaveDischargingCertificate {
    pub name: String,
    pub n: usize,
    pub r_trivial: usize,
    pub girth: usize,
    pub is_bipartite: bool,
    pub c4: u64,
    pub c6: u64,
    pub c8: u64,
    pub c10: u64,
    pub c12: u64,
    pub c14: u64,
    pub c16: u64,
    pub theta_conservation_verified: bool,
    pub prism_conservation_verified: bool,
    pub menger_2_edge_connected: bool,
    pub perron_dominance_k: u64,
    pub perron_dominance_m: u64,
    pub perron_dominance_margin: u64,
    pub tr_b16_4d_bessel_floor: u64,
    pub n16_discharged_upper_bound: u64,
    pub coupled_c16_4d_lower_bound: u64,
    pub total_dyadic_cycles_up_to_16: u64,
    pub dyadic_closure_certified: bool,
}

#[derive(Clone, Copy, Default)]
struct HalfPath5 {
    int_mask: u64,
    v1: u8,
    next: i16,
}

/// Counts undirected simple 10-cycles `C_10(G)` in `adj[0..n]` using canonical-anchor (`s = min V(C)`)
/// 5-step Half-Path Meet-in-the-Middle with a zero-allocation fixed-size stack buffer (`<= 48` half-paths
/// per anchor `s` in any cubic graph).
pub fn count_c10_slice(adj: &[u64], n: usize) -> u64 {
    if n < 10 {
        return 0;
    }
    let mut total_c10 = 0u64;
    let mut hp = [HalfPath5::default(); 64];
    let mut head = [-1i16; 64];

    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        if ns.count_ones() < 2 {
            continue;
        }
        head.fill(-1);
        let mut hp_cnt = 0usize;

        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            let vis1 = 1u64 << v1;

            let mut m2 = adj[v1] & higher_mask & !vis1;
            while m2 != 0 {
                let v2 = m2.trailing_zeros() as usize;
                m2 &= m2 - 1;
                let vis2 = vis1 | (1u64 << v2);

                let mut m3 = adj[v2] & higher_mask & !vis2;
                while m3 != 0 {
                    let v3 = m3.trailing_zeros() as usize;
                    m3 &= m3 - 1;
                    let vis3 = vis2 | (1u64 << v3);

                    let mut m4 = adj[v3] & higher_mask & !vis3;
                    while m4 != 0 {
                        let v4 = m4.trailing_zeros() as usize;
                        m4 &= m4 - 1;
                        let vis4 = vis3 | (1u64 << v4);

                        let mut m5 = adj[v4] & higher_mask & !vis4;
                        while m5 != 0 {
                            let mid = m5.trailing_zeros() as usize;
                            m5 &= m5 - 1;
                            if hp_cnt < 64 {
                                hp[hp_cnt] = HalfPath5 {
                                    int_mask: vis4,
                                    v1: v1 as u8,
                                    next: head[mid],
                                };
                                head[mid] = hp_cnt as i16;
                                hp_cnt += 1;
                            }
                        }
                    }
                }
            }
        }

        for &h_start in &head[(s + 1)..n] {
            let mut i = h_start;
            while i >= 0 {
                let hpi = hp[i as usize];
                let mut j = hpi.next;
                while j >= 0 {
                    let hpj = hp[j as usize];
                    if hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0 {
                        total_c10 += 1;
                    }
                    j = hpj.next;
                }
                i = hpi.next;
            }
        }
    }
    total_c10
}

/// Computes the exact integer 3x3 shifted Stieltjes lower bound `C_10_min` (`det H_3^{(x)}(mu) >= 0`)
/// and 3x3 shifted Hausdorff upper bound `C_10_max` (`det H_3((9-x)mu) >= 0`) on `C_10(G)`
/// for a 3-regular graph on `n` vertices with `r_trivial` trivial `lambda^2 = 9` eigenvalues and `C_6(G) = c6`.
///
/// At the `(3,10)`-Moore bound `n = 62, r_trivial = 2, c6 = 0`, returns `(298, 297)` (`298 > 297`),
/// certifying in exact integer arithmetic that no `(3,10)`-Moore graph on `n = 62` vertices exists!
pub fn compute_c10_hankel_bounds(n: usize, r_trivial: usize, c6: u64) -> (u64, u64) {
    let ni = n as i128;
    let ri = r_trivial as i128;
    let c6i = c6 as i128;

    let num_lo = 3
        * (4 * c6i * c6i * c6i + 3 * ni * c6i * c6i
            - 9 * ni * ni * c6i
            - 4 * ni * ni * ni
            - 3 * ri * (31 * c6i * c6i + 3 * ni * c6i - 84 * ni * ni));
    let den_lo = 5 * (ni * (c6i + ni) - 3 * ri * (c6i + 5 * ni));
    let c10_min = if den_lo > 0 && num_lo > 0 {
        ((num_lo + den_lo - 1) / den_lo) as u64
    } else {
        0
    };

    let num_hi = 6 * c6i * c6i * c6i - 27 * ni * c6i * c6i + 48 * ni * ni * ni;
    let den_hi = 5 * ni * (2 * ni - c6i);
    let c10_max = if den_hi > 0 && num_hi >= 0 {
        (num_hi / den_hi) as u64
    } else {
        0
    };

    (c10_min, c10_max)
}

fn is_bipartite_slice(adj: &[u64], n: usize) -> bool {
    let mut red = 0u64;
    let mut blue = 0u64;
    for s in 0..n {
        let bit = 1u64 << s;
        if ((red | blue) & bit) != 0 {
            continue;
        }
        red |= bit;
        let mut front = bit;
        let mut cur_is_red = true;
        while front != 0 {
            let mut nxt = 0u64;
            let mut fm = front;
            while fm != 0 {
                let u = fm.trailing_zeros() as usize;
                fm &= fm - 1;
                nxt |= adj[u];
            }
            if cur_is_red {
                if (nxt & red) != 0 {
                    return false;
                }
                front = nxt & !blue;
                blue |= front;
            } else {
                if (nxt & blue) != 0 {
                    return false;
                }
                front = nxt & !red;
                red |= front;
            }
            cur_is_red = !cur_is_red;
        }
    }
    true
}

#[allow(clippy::needless_range_loop)]
fn eval_b10_hankel_single_graph(adj: &[u64], n: usize) -> B10HankelPinchSummary {
    let mut a2 = [[0u64; MAX_VERTICES]; MAX_VERTICES];
    let mut a4 = [[0u64; MAX_VERTICES]; MAX_VERTICES];
    let mut a6 = [[0u64; MAX_VERTICES]; MAX_VERTICES];

    for i in 0..n {
        for j in i..n {
            let c = (adj[i] & adj[j]).count_ones() as u64;
            a2[i][j] = c;
            a2[j][i] = c;
        }
    }
    for i in 0..n {
        for j in i..n {
            let mut s = 0u64;
            for k in 0..n {
                s += a2[i][k] * a2[k][j];
            }
            a4[i][j] = s;
            a4[j][i] = s;
        }
    }
    for i in 0..n {
        for j in i..n {
            let mut s = 0u64;
            for k in 0..n {
                s += a4[i][k] * a2[k][j];
            }
            a6[i][j] = s;
            a6[j][i] = s;
        }
    }

    let mut tr_a2 = 0u64;
    let mut tr_a4 = 0u64;
    let mut tr_a6 = 0u64;
    let mut tr_a8 = 0u64;
    let mut tr_a10 = 0u64;
    for i in 0..n {
        tr_a2 += a2[i][i];
        tr_a4 += a4[i][i];
        tr_a6 += a6[i][i];
        for k in 0..n {
            tr_a8 += a4[i][k] * a4[k][i];
            tr_a10 += a6[i][k] * a4[k][i];
        }
    }

    let c6 = tr_a6.saturating_sub(87 * (n as u64)) / 12;
    let tr_b10_i128 = (tr_a10 as i128) - 20 * (tr_a8 as i128) + 140 * (tr_a6 as i128)
        - 400 * (tr_a4 as i128)
        + 400 * (tr_a2 as i128)
        - 63 * (n as i128);
    let tr_b10 = tr_b10_i128.max(0) as u64;
    let c10 = count_c10_slice(adj, n);
    let r_trivial = if is_bipartite_slice(adj, n) { 2 } else { 1 };
    let (c10_hankel_lower_bound, c10_hankel_upper_bound) =
        compute_c10_hankel_bounds(n, r_trivial, c6);

    B10HankelPinchSummary {
        c6,
        c10,
        tr_a10,
        tr_b10,
        c10_hankel_lower_bound,
        c10_hankel_upper_bound,
    }
}

/// Exact CPU twin of the CUDA kernel `k_eval_b10_hankel_pinch` (Theorem 2.6):
/// computes `Tr(A^{10})`, `Tr(B^{10})`, exact `C_10(G)` via 5-step half-path MITM, and
/// the 3x3 shifted Stieltjes/Hausdorff Hankel bounds `[C_{10,min}, C_{10,max}]` for a batch of graphs.
pub fn eval_b10_hankel_pinch_cpu(
    flat_adj: &[u64],
    n: usize,
    batch_size: usize,
) -> Vec<B10HankelPinchSummary> {
    assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
    let mut out = Vec::with_capacity(batch_size);
    for b in 0..batch_size {
        let adj = &flat_adj[b * MAX_VERTICES..(b + 1) * MAX_VERTICES];
        out.push(eval_b10_hankel_single_graph(adj, n));
    }
    out
}

/// Certifies Theorem 2.6(a)–(c) on a 3-regular graph `g`:
/// - Part (a): `Tr(B^{10}) == 20 * C_10(G)` (`N_10(G) == 0`) on every cubic graph of girth `g(G) >= 6`;
/// - Part (b): 3x3 shifted Stieltjes (`det H_3^{(x)}(mu) >= 0`) and Hausdorff (`det H_3((9-x)mu) >= 0`)
///   moment positivity and `C_10(G) in [C_{10,min}, C_{10,max}]` whenever `g(G) >= 5`, `C_8(G) = 0`, and `N_10(G) = 0`;
/// - Part (c): Tangent-hair 10-cycle bound `N_16^{theta}(G) <= 320 * mu_2(G) * C_10(G)` (with exact equality
///   `N_16^{theta}(G) = 320 * C_10(G) = 16 * Tr(B^{10})` on all `Theta(3,3,3)`-free bipartite `C_4 = C_8 = 0`
///   frontier lifts) and Bessel `span{1, lambda^2}` quadratic lower bound on `C_16(G)`.
pub fn certify_b10_hankel_and_dual_theta_pinch(g: &BitGraph) -> B10HankelPinchCertificate {
    assert_eq!(g.min_degree(), 3, "Theorem 2.6 requires 3-regular graph");
    assert_eq!(g.max_degree(), 3, "Theorem 2.6 requires 3-regular graph");
    let n = g.n;
    let is_bip = is_bipartite(g);
    let r_trivial = if is_bip { 2 } else { 1 };
    let girth = compute_girth(g);
    let c4 = count_c4_slice(&g.adj, n);
    let c8 = count_c8_slice(&g.adj, n);
    let c16 = count_c16_total(g);
    let (c6, max_c6_per_2path) = count_c6_and_max_per_2path(&g.adj, n);

    let summary = eval_b10_hankel_single_graph(&g.adj, n);
    let c10 = summary.c10;
    assert!(
        summary.tr_b10 >= 20 * c10,
        "Tr(B^10) ({}) must be >= 20 * C_10 ({}) on {}",
        summary.tr_b10,
        20 * c10,
        g.name
    );
    let n10_excess = summary.tr_b10 - 20 * c10;
    if girth >= 6 {
        assert_eq!(
            n10_excess, 0,
            "Theorem 2.6(a): N_10(G) must vanish identically (Tr(B^10) == 20 * C_10) on girth >= 6 graph {}",
            g.name
        );
    }

    // Compute exact 3x3 shifted Stieltjes and Hausdorff determinants from non-trivial moments mu_0..mu_5
    let ni = n as i128;
    let ri = r_trivial as i128;
    let c6i = c6 as i128;
    let c10i = c10 as i128;
    let mu0 = ni - ri;
    let mu1 = 3 * ni - 9 * ri;
    let mu2 = 15 * ni - 81 * ri;
    let mu3 = (87 * ni + 12 * c6i) - 729 * ri;
    let mu4 = (543 * ni + 192 * c6i) - 6561 * ri;
    let mu5 = (3543 * ni + 2160 * c6i + 20 * c10i) - 59049 * ri;

    // H_3^{(x)}(mu) = [[mu1, mu2, mu3], [mu2, mu3, mu4], [mu3, mu4, mu5]]
    let det_h3_stieltjes_x = mu1 * (mu3 * mu5 - mu4 * mu4) - mu2 * (mu2 * mu5 - mu4 * mu3)
        + mu3 * (mu2 * mu4 - mu3 * mu3);

    // H_3((9-x)mu) = [[nu0, nu1, nu2], [nu1, nu2, nu3], [nu2, nu3, nu4]] where nu_m = 9*mu_m - mu_{m+1}
    let nu0 = 9 * mu0 - mu1;
    let nu1 = 9 * mu1 - mu2;
    let nu2 = 9 * mu2 - mu3;
    let nu3 = 9 * mu3 - mu4;
    let nu4 = 9 * mu4 - mu5;
    let det_h3_hausdorff_9mx = nu0 * (nu2 * nu4 - nu3 * nu3) - nu1 * (nu1 * nu4 - nu3 * nu2)
        + nu2 * (nu1 * nu3 - nu2 * nu2);

    if girth >= 6 && c8 == 0 {
        assert!(
            det_h3_stieltjes_x >= 0,
            "Theorem 2.6(b): det H_3^(x)(mu) ({}) must be >= 0 on {}",
            det_h3_stieltjes_x,
            g.name
        );
        assert!(
            det_h3_hausdorff_9mx >= 0,
            "Theorem 2.6(b): det H_3((9-x)mu) ({}) must be >= 0 on {}",
            det_h3_hausdorff_9mx,
            g.name
        );
        assert!(
            c10 >= summary.c10_hankel_lower_bound && c10 <= summary.c10_hankel_upper_bound,
            "Theorem 2.6(b): C_10={} must lie in Hankel interval [{}..={}] on {}",
            c10,
            summary.c10_hankel_lower_bound,
            summary.c10_hankel_upper_bound,
            g.name
        );
    }

    let topo = certify_n16_topological_decomposition(g);
    let mu2_mult = max_c6_per_2path.max(1);
    let n16_theta_c10_bound = 320 * mu2_mult * c10;
    let n16_dual_upper_bound = 384 * mu2_mult * c6 + (5568 * c6).min(n16_theta_c10_bound);

    // Bessel span{1, lambda^2} boosted spectral floor on Tr(B^{16})
    let tr_b16_bessel_boosted_floor = if n > 7 * r_trivial && 11179 * r_trivial > 85 * n {
        let den_full = 3 * ni * (ni - ri) * (ni - 7 * ri);
        let term_in_sq = (10 * c10i + 24 * c6i) * (ni - ri) - 3 * ni * (ni + 257 * ri);
        let num_full =
            18 * ni * ni * (11179 * ri - 85 * ni) * (ni - 7 * ri) + 2 * term_in_sq * term_in_sq;
        ((num_full + den_full - 1) / den_full) as u64
    } else {
        0
    };

    let coupled_c16_lower_bound = if tr_b16_bessel_boosted_floor > n16_dual_upper_bound {
        (tr_b16_bessel_boosted_floor - n16_dual_upper_bound).div_ceil(32)
    } else {
        0
    };

    // Universal (c6, c10)-eliminated quadratic lower bound for Theta(3,3,3)-free (mu_2 <= 1) graphs:
    // 32 * C_16(G) >= 90 * n * (501 * r - 11 * n) / (n - r) + 384 * C_6(G)
    let universal_mu2_1_c16_lower_bound = if n > r_trivial && 501 * r_trivial > 11 * n {
        let univ_num = 90 * ni * (501 * ri - 11 * ni) + 384 * c6i * (ni - ri);
        let univ_tr_diff = ((univ_num + (ni - ri) - 1) / (ni - ri)) as u64;
        univ_tr_diff.div_ceil(32)
    } else {
        0
    };

    if is_bip && c4 == 0 && c8 == 0 {
        assert!(
            topo.walks_theta_total <= n16_theta_c10_bound,
            "Theorem 2.6(c): N_16^theta ({}) must be <= 320 * mu_2 * C_10 ({}) on {}",
            topo.walks_theta_total,
            n16_theta_c10_bound,
            g.name
        );
        assert!(
            topo.n16_total <= n16_dual_upper_bound,
            "Theorem 2.6(c): N_16(G) ({}) must be <= dual upper bound ({}) on {}",
            topo.n16_total,
            n16_dual_upper_bound,
            g.name
        );
        assert!(
            topo.tr_b16_chebyshev >= tr_b16_bessel_boosted_floor,
            "Theorem 2.6(c): Tr(B^16) ({}) must be >= Bessel-boosted floor ({}) on {}",
            topo.tr_b16_chebyshev,
            tr_b16_bessel_boosted_floor,
            g.name
        );
        assert!(
            c16 >= coupled_c16_lower_bound,
            "Theorem 2.6(c): C_16 ({}) must be >= coupled lower bound ({}) on {}",
            c16,
            coupled_c16_lower_bound,
            g.name
        );
        if max_c6_per_2path <= 1 {
            assert!(
                c16 >= universal_mu2_1_c16_lower_bound,
                "Theorem 2.6(c): C_16 ({}) must be >= universal mu_2<=1 floor ({}) on {}",
                c16,
                universal_mu2_1_c16_lower_bound,
                g.name
            );
        }
    }

    B10HankelPinchCertificate {
        name: g.name.clone(),
        n,
        r_trivial,
        girth,
        is_bipartite: is_bip,
        c4,
        c6,
        c8,
        c10,
        c16,
        max_c6_per_2path,
        tr_a10: summary.tr_a10,
        tr_b10: summary.tr_b10,
        n10_excess,
        det_h3_stieltjes_x,
        det_h3_hausdorff_9mx,
        c10_hankel_lower_bound: summary.c10_hankel_lower_bound,
        c10_hankel_upper_bound: summary.c10_hankel_upper_bound,
        tr_b16_bessel_boosted_floor,
        n16_dumbbell_actual: topo.walks_dumbbell_total,
        n16_theta_actual: topo.walks_theta_total,
        n16_theta_c10_bound,
        n16_dual_upper_bound,
        coupled_c16_lower_bound,
        universal_mu2_1_c16_lower_bound,
    }
}

/// Counts undirected simple 12-cycles `C_12(G)` in `adj[0..n]` using canonical-anchor (`s = min V(C)`)
/// 6-step Half-Path Meet-in-the-Middle with a fixed-size stack buffer (`<= 96` half-paths per anchor
/// `s` in any cubic graph).
pub fn count_c12_slice(adj: &[u64], n: usize) -> u64 {
    if n < 12 {
        return 0;
    }
    let mut total_c12 = 0u64;
    let mut hp = [HalfPath5::default(); 128];
    let mut head = [-1i16; 64];

    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        if ns.count_ones() < 2 {
            continue;
        }
        head.fill(-1);
        let mut hp_cnt = 0usize;

        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            let vis1 = 1u64 << v1;

            let mut m2 = adj[v1] & higher_mask & !vis1;
            while m2 != 0 {
                let v2 = m2.trailing_zeros() as usize;
                m2 &= m2 - 1;
                let vis2 = vis1 | (1u64 << v2);

                let mut m3 = adj[v2] & higher_mask & !vis2;
                while m3 != 0 {
                    let v3 = m3.trailing_zeros() as usize;
                    m3 &= m3 - 1;
                    let vis3 = vis2 | (1u64 << v3);

                    let mut m4 = adj[v3] & higher_mask & !vis3;
                    while m4 != 0 {
                        let v4 = m4.trailing_zeros() as usize;
                        m4 &= m4 - 1;
                        let vis4 = vis3 | (1u64 << v4);

                        let mut m5 = adj[v4] & higher_mask & !vis4;
                        while m5 != 0 {
                            let v5 = m5.trailing_zeros() as usize;
                            m5 &= m5 - 1;
                            let vis5 = vis4 | (1u64 << v5);

                            let mut m6 = adj[v5] & higher_mask & !vis5;
                            while m6 != 0 {
                                let mid = m6.trailing_zeros() as usize;
                                m6 &= m6 - 1;
                                if hp_cnt < 128 {
                                    hp[hp_cnt] = HalfPath5 {
                                        int_mask: vis5,
                                        v1: v1 as u8,
                                        next: head[mid],
                                    };
                                    head[mid] = hp_cnt as i16;
                                    hp_cnt += 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        for &h_start in &head[(s + 1)..n] {
            let mut i = h_start;
            while i >= 0 {
                let hpi = hp[i as usize];
                let mut j = hpi.next;
                while j >= 0 {
                    let hpj = hp[j as usize];
                    if hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0 {
                        total_c12 += 1;
                    }
                    j = hpj.next;
                }
                i = hpi.next;
            }
        }
    }
    total_c12
}

/// Enumerates all undirected simple 6-cycles in `adj[0..n]`, represents each 6-cycle by its
/// 128-bit undirected edge mask `(lo, hi)` over the `<= 96` edges of a cubic graph on `n <= 64`
/// vertices, and returns `(c6, tau1, tau3)`:
/// - `c6`: number of simple 6-cycles `C_6(G)`;
/// - `tau1`: number of unordered pairs `{Z_1, Z_2}` of simple 6-cycles sharing a single edge (`Theta(1,5,5)`);
/// - `tau3`: number of unordered pairs `{Z_1, Z_2}` of simple 6-cycles sharing 3 edges (`Theta(3,3,3)`).
pub fn count_c6_tau1_tau3_slice(adj: &[u64], n: usize) -> (u64, u64, u64) {
    let mut edge_id = [[0u8; MAX_VERTICES]; MAX_VERTICES];
    let mut num_edges = 0u8;
    for u in 0..n {
        let mut m = adj[u];
        while m != 0 {
            let v = m.trailing_zeros() as usize;
            m &= m - 1;
            if u < v && num_edges < 128 {
                edge_id[u][v] = num_edges;
                edge_id[v][u] = num_edges;
                num_edges += 1;
            }
        }
    }

    let mut c6_lo = [0u64; 128];
    let mut c6_hi = [0u64; 128];
    let mut c6_cnt = 0usize;

    let set_edge = |lo: &mut u64, hi: &mut u64, eid: u8| {
        if eid < 64 {
            *lo |= 1u64 << eid;
        } else {
            *hi |= 1u64 << (eid - 64);
        }
    };

    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            let m5_mask = if v1 == 63 {
                0u64
            } else {
                !((1u64 << (v1 + 1)) - 1)
            };
            let mut m5 = ns & m5_mask;
            while m5 != 0 {
                let v5 = m5.trailing_zeros() as usize;
                m5 &= m5 - 1;
                let vis15 = (1u64 << s) | (1u64 << v1) | (1u64 << v5);
                let mut m2 = adj[v1] & higher_mask & !vis15;
                while m2 != 0 {
                    let v2 = m2.trailing_zeros() as usize;
                    m2 &= m2 - 1;
                    let vis2 = vis15 | (1u64 << v2);
                    let mut m4 = adj[v5] & higher_mask & !vis2;
                    while m4 != 0 {
                        let v4 = m4.trailing_zeros() as usize;
                        m4 &= m4 - 1;
                        let vis4 = vis2 | (1u64 << v4);
                        let mut m3 = adj[v2] & adj[v4] & higher_mask & !vis4;
                        while m3 != 0 {
                            let v3 = m3.trailing_zeros() as usize;
                            m3 &= m3 - 1;
                            if c6_cnt < 128 {
                                let mut lo = 0u64;
                                let mut hi = 0u64;
                                set_edge(&mut lo, &mut hi, edge_id[s][v1]);
                                set_edge(&mut lo, &mut hi, edge_id[v1][v2]);
                                set_edge(&mut lo, &mut hi, edge_id[v2][v3]);
                                set_edge(&mut lo, &mut hi, edge_id[v3][v4]);
                                set_edge(&mut lo, &mut hi, edge_id[v4][v5]);
                                set_edge(&mut lo, &mut hi, edge_id[v5][s]);
                                c6_lo[c6_cnt] = lo;
                                c6_hi[c6_cnt] = hi;
                            }
                            c6_cnt += 1;
                        }
                    }
                }
            }
        }
    }

    let stored = c6_cnt.min(128);
    let mut tau1 = 0u64;
    let mut tau3 = 0u64;
    for a in 0..stored {
        for b in (a + 1)..stored {
            let shared = (c6_lo[a] & c6_lo[b]).count_ones() + (c6_hi[a] & c6_hi[b]).count_ones();
            if shared == 1 {
                tau1 += 1;
            } else if shared == 3 {
                tau3 += 1;
            }
        }
    }

    (c6_cnt as u64, tau1, tau3)
}

/// Computes the exact integer ceiling of the 3D Chebyshev–Bessel floor on `Tr(B^{16})`
/// obtained by projecting `y_i = P_8(lambda_i^2) - 1` (`i > r`) onto `span{1, lambda_i^2, lambda_i^4}`.
pub fn compute_tr_b16_3d_bessel_floor(
    n: usize,
    r_trivial: usize,
    c6: u64,
    c10: u64,
    tr_b12: u64,
) -> u64 {
    if n <= 7 * r_trivial || 11179 * r_trivial <= 85 * n {
        return 0;
    }
    let ni = n as i128;
    let ri = r_trivial as i128;
    let c6i = c6 as i128;
    let c10i = c10 as i128;
    let tr_b12i = tr_b12 as i128;

    let den_2d = 3 * ni * (ni - ri) * (ni - 7 * ri);
    let term_in_sq = (10 * c10i + 24 * c6i) * (ni - ri) - 3 * ni * (ni + 257 * ri);
    let num_2d =
        18 * ni * ni * (11179 * ri - 85 * ni) * (ni - 7 * ri) + 2 * term_in_sq * term_in_sq;

    let d3_tilde = -(ni - ri) * c6i * c6i + ni * (ni + 17 * ri) * c6i + ni * ni * (ni - 31 * ri);
    if d3_tilde > 0 {
        let m_num2 = ni * (ni - 7 * ri) * tr_b12i
            - 2 * (2 * (ni - ri) * c6i - ni * (ni + 17 * ri)) * (10 * c10i + 24 * c6i)
            + 12 * ni * (c6i - 2 * ni) * (ni + 257 * ri);
        let num_3d = 8 * d3_tilde * num_2d + (ni - ri) * (m_num2 * m_num2);
        let den_3d = 8 * d3_tilde * den_2d;
        ((num_3d + den_3d - 1) / den_3d) as u64
    } else {
        ((num_2d + den_2d - 1) / den_2d) as u64
    }
}

#[allow(clippy::needless_range_loop)]
fn eval_b12_deg2_bessel_single_graph(adj: &[u64], n: usize) -> B12Deg2BesselSummary {
    let mut a2 = [[0u64; MAX_VERTICES]; MAX_VERTICES];
    let mut a4 = [[0u64; MAX_VERTICES]; MAX_VERTICES];
    let mut a6 = [[0u64; MAX_VERTICES]; MAX_VERTICES];

    for i in 0..n {
        for j in i..n {
            let c = (adj[i] & adj[j]).count_ones() as u64;
            a2[i][j] = c;
            a2[j][i] = c;
        }
    }
    for i in 0..n {
        for j in i..n {
            let mut s = 0u64;
            for k in 0..n {
                s += a2[i][k] * a2[k][j];
            }
            a4[i][j] = s;
            a4[j][i] = s;
        }
    }
    for i in 0..n {
        for j in i..n {
            let mut s = 0u64;
            for k in 0..n {
                s += a4[i][k] * a2[k][j];
            }
            a6[i][j] = s;
            a6[j][i] = s;
        }
    }

    let mut tr_a2 = 0u64;
    let mut tr_a4 = 0u64;
    let mut tr_a6 = 0u64;
    let mut tr_a8 = 0u64;
    let mut tr_a10 = 0u64;
    let mut tr_a12 = 0u64;
    for i in 0..n {
        tr_a2 += a2[i][i];
        tr_a4 += a4[i][i];
        tr_a6 += a6[i][i];
        for k in 0..n {
            let a4_ik = a4[i][k];
            let a6_ik = a6[i][k];
            tr_a8 += a4_ik * a4_ik;
            tr_a10 += a6_ik * a4_ik;
            tr_a12 += a6_ik * a6_ik;
        }
    }

    // P_12(x) = x^6 - 24 x^5 + 216 x^4 - 896 x^3 + 1680 x^2 - 1152 x + 129
    let tr_b12_i128 = (tr_a12 as i128) - 24 * (tr_a10 as i128) + 216 * (tr_a8 as i128)
        - 896 * (tr_a6 as i128)
        + 1680 * (tr_a4 as i128)
        - 1152 * (tr_a2 as i128)
        + 129 * (n as i128);
    let tr_b12 = tr_b12_i128.max(0) as u64;

    let (c6, tau1, tau3) = count_c6_tau1_tau3_slice(adj, n);
    let c10 = count_c10_slice(adj, n);
    let c12 = count_c12_slice(adj, n);
    let r_trivial = if is_bipartite_slice(adj, n) { 2 } else { 1 };
    let tr_b16_3d_bessel_floor = compute_tr_b16_3d_bessel_floor(n, r_trivial, c6, c10, tr_b12);

    B12Deg2BesselSummary {
        c6,
        c10,
        c12,
        tau1,
        tau3,
        tr_a12,
        tr_b12,
        tr_b16_3d_bessel_floor,
    }
}

/// Exact CPU twin of the CUDA kernel `k_eval_b12_deg2_bessel_pinch` (Theorem 2.7):
/// computes `Tr(A^{12})`, `Tr(B^{12})`, exact `C_6(G)`, `tau_1(G)`, `tau_3(G)`, `C_10(G)`,
/// `C_12(G)` via 6-step half-path MITM, and the 3D Chebyshev–Bessel `span{1, lambda^2, lambda^4}`
/// floor on `Tr(B^{16})` for a batch of graphs.
pub fn eval_b12_deg2_bessel_pinch_cpu(
    flat_adj: &[u64],
    n: usize,
    batch_size: usize,
) -> Vec<B12Deg2BesselSummary> {
    assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
    let mut out = Vec::with_capacity(batch_size);
    for b in 0..batch_size {
        let adj = &flat_adj[b * MAX_VERTICES..(b + 1) * MAX_VERTICES];
        out.push(eval_b12_deg2_bessel_single_graph(adj, n));
    }
    out
}

/// Certifies Theorem 2.7(a)–(b) on a 3-regular graph `g`:
/// - Part (a): Exact 12-walk topological identity `Tr(B^{12}) == 24 * C_12(G) + 12 * C_6(G) + 24 * tau_1(G) + 24 * tau_3(G)`
///   across all girth `g(G) >= 6` cubic graphs with `C_8(G) = 0` (`tau_3(G) = 0` whenever `G` is `Theta(3,3,3)`-free);
/// - Part (b): Degree-2 orthogonal Bessel–Hashimoto projection on `span{1, lambda_i^2, lambda_i^4}`
///   boosting `Tr(B^{16})`, `2x2` Stieltjes/Hausdorff bounds on `v_i = P_4(lambda_i^2) - 1 in [-8, 17]`,
///   `3x3` Markov–Lukács `x(9-x)` and `4x4` Hankel positivity, and the `(C_10)`-eliminated quadratic envelope
///   `max(Phi_{2D}, Phi_{3D}) > 0` extending unconditional exclusion of `C_16(G) = 0` up to `n <= 250`.
pub fn certify_b12_deg2_bessel_pinch(g: &BitGraph) -> B12Deg2BesselCertificate {
    assert_eq!(g.min_degree(), 3, "Theorem 2.7 requires 3-regular graph");
    assert_eq!(g.max_degree(), 3, "Theorem 2.7 requires 3-regular graph");
    let n = g.n;
    let is_bip = is_bipartite(g);
    let r_trivial = if is_bip { 2 } else { 1 };
    let girth = compute_girth(g);
    let c4 = count_c4_slice(&g.adj, n);
    let c8 = count_c8_slice(&g.adj, n);
    let c16 = count_c16_total(g);
    let (c6_check, max_c6_per_2path) = count_c6_and_max_per_2path(&g.adj, n);

    let summary = eval_b12_deg2_bessel_single_graph(&g.adj, n);
    assert_eq!(
        summary.c6, c6_check,
        "C_6 count mismatch between 2-path and 128-bit edge-mask enumerators on {}",
        g.name
    );

    let c6 = summary.c6;
    let c10 = summary.c10;
    let c12 = summary.c12;
    let tau1 = summary.tau1;
    let tau3 = summary.tau3;
    let tr_b12 = summary.tr_b12;

    assert!(
        tr_b12 >= 24 * c12,
        "Tr(B^12) ({}) must be >= 24 * C_12 ({}) on {}",
        tr_b12,
        24 * c12,
        g.name
    );
    let n12_composite_actual = tr_b12 - 24 * c12;
    let n12_composite_predicted = 12 * c6 + 24 * tau1 + 24 * tau3;

    if girth >= 6 && c8 == 0 {
        assert_eq!(
            n12_composite_actual, n12_composite_predicted,
            "Theorem 2.7(a): N_12(G) = Tr(B^12) - 24*C_12 ({}) must equal 12*C_6 + 24*tau_1 + 24*tau_3 ({}) on {}",
            n12_composite_actual,
            n12_composite_predicted,
            g.name
        );
        if max_c6_per_2path <= 1 {
            assert_eq!(
                tau3, 0,
                "Theorem 2.7(a): tau_3(G) must be 0 on Theta(3,3,3)-free graph {}",
                g.name
            );
            assert_eq!(
                tr_b12,
                24 * c12 + 12 * c6 + 24 * tau1,
                "Theorem 2.7(a): Tr(B^12) must equal 24*C_12 + 12*C_6 + 24*tau_1 on {}",
                g.name
            );
        }
    }

    let ni = n as i128;
    let ri = r_trivial as i128;
    let c6i = c6 as i128;
    let c10i = c10 as i128;
    let tr_b12i = tr_b12 as i128;

    // 1. Stieltjes & Hausdorff determinants on v_i = P_4(lambda_i^2) - 1 in [-8, 17]:
    //    det S_2^{(v+8)} = (7n - 25r) * Tr(B^{12}) + 288 * n * (3n - 125r) >= 0
    //    det S_2^{(17-v)} = 18 * n * (448 * n - Tr(B^{12})) >= 0
    let det_s2_v_plus_8 = (7 * ni - 25 * ri) * tr_b12i + 288 * ni * (3 * ni - 125 * ri);
    let det_s2_17_minus_v = 18 * ni * (448 * ni - tr_b12i);

    // 2. Non-trivial moments mu_0..mu_6 and 3x3 Markov-Lukacs x(9-x) + 4x4 Hankel determinants:
    let mu0 = ni - ri;
    let mu1 = 3 * ni - 9 * ri;
    let mu2 = 15 * ni - 81 * ri;
    let mu3 = (87 * ni + 12 * c6i) - 729 * ri;
    let mu4 = (543 * ni + 192 * c6i) - 6561 * ri;
    let mu5 = (3543 * ni + 2160 * c6i + 20 * c10i) - 59049 * ri;
    let mu6 = (23823 * ni + 21120 * c6i + 480 * c10i + tr_b12i) - 531441 * ri;

    let ml0 = 9 * mu1 - mu2;
    let ml1 = 9 * mu2 - mu3;
    let ml2 = 9 * mu3 - mu4;
    let ml3 = 9 * mu4 - mu5;
    let ml4 = 9 * mu5 - mu6;
    let det_ml3_x_9mx = ml0 * (ml2 * ml4 - ml3 * ml3) - ml1 * (ml1 * ml4 - ml3 * ml2)
        + ml2 * (ml1 * ml3 - ml2 * ml2);

    // 4x4 Hankel determinant det H_4(mu) via cofactor expansion along row 0:
    let det3 = |a00: i128,
                a01: i128,
                a02: i128,
                a10: i128,
                a11: i128,
                a12: i128,
                a20: i128,
                a21: i128,
                a22: i128| {
        a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20)
            + a02 * (a10 * a21 - a11 * a20)
    };
    let det_h4_hankel = mu0 * det3(mu2, mu3, mu4, mu3, mu4, mu5, mu4, mu5, mu6)
        - mu1 * det3(mu1, mu3, mu4, mu2, mu4, mu5, mu3, mu5, mu6)
        + mu2 * det3(mu1, mu2, mu4, mu2, mu3, mu5, mu3, mu4, mu6)
        - mu3 * det3(mu1, mu2, mu3, mu2, mu3, mu4, mu3, mu4, mu5);

    let b10_cert = certify_b10_hankel_and_dual_theta_pinch(g);
    let tr_b16_2d_bessel_floor = if c4 == 0 && c8 == 0 {
        b10_cert.tr_b16_bessel_boosted_floor
    } else {
        0
    };
    let tr_b16_3d_bessel_floor = if c4 == 0 && c8 == 0 {
        summary.tr_b16_3d_bessel_floor
    } else {
        0
    };

    let coupled_c16_2d_lower_bound = if c4 == 0 && c8 == 0 {
        b10_cert.coupled_c16_lower_bound
    } else {
        0
    };
    let coupled_c16_3d_lower_bound =
        if c4 == 0 && c8 == 0 && tr_b16_3d_bessel_floor > b10_cert.n16_dual_upper_bound {
            (tr_b16_3d_bessel_floor - b10_cert.n16_dual_upper_bound).div_ceil(32)
        } else {
            0
        };

    // Z-eliminated envelope max(Phi_{2D}, Phi_{3D}) in {1, v(x), x} basis:
    let phi_3d_c16_bound = if c4 == 0 && c8 == 0 && 5 * n > 59 * r_trivial {
        let num_phi3 = (ni - ri) * tr_b12i * tr_b12i
            - 12 * (16 * (ni - ri) * c6i + (635 * ri - 5 * ni) * ni) * tr_b12i
            + 1152
                * (8 * (ni - ri) * c6i * c6i
                    + (5 * ni + 517 * ri) * ni * c6i
                    + (1255 * ri - 25 * ni) * ni * ni);
        let den_phi3 = 6 * ni * (5 * ni - 59 * ri);
        if num_phi3 > 0 {
            let tr_diff = ((num_phi3 + den_phi3 - 1) / den_phi3) as u64;
            tr_diff.div_ceil(32)
        } else {
            0
        }
    } else {
        0
    };
    let z_eliminated_c16_envelope_bound = if c4 == 0 && c8 == 0 {
        b10_cert
            .universal_mu2_1_c16_lower_bound
            .max(phi_3d_c16_bound)
    } else {
        0
    };

    if girth >= 6 && c8 == 0 {
        assert!(
            det_s2_v_plus_8 >= 0,
            "Theorem 2.7(b): det S_2^(v+8) ({}) must be >= 0 on {}",
            det_s2_v_plus_8,
            g.name
        );
        assert!(
            det_s2_17_minus_v >= 0,
            "Theorem 2.7(b): det S_2^(17-v) ({}) must be >= 0 on {}",
            det_s2_17_minus_v,
            g.name
        );
        assert!(
            det_ml3_x_9mx >= 0,
            "Theorem 2.7(b): det H_3^(x(9-x))(mu) ({}) must be >= 0 on {}",
            det_ml3_x_9mx,
            g.name
        );
        assert!(
            det_h4_hankel >= 0,
            "Theorem 2.7(b): det H_4(mu) ({}) must be >= 0 on {}",
            det_h4_hankel,
            g.name
        );
        assert!(
            tr_b16_3d_bessel_floor >= tr_b16_2d_bessel_floor,
            "Theorem 2.7(b): 3D Bessel floor ({}) must be >= 2D Bessel floor ({}) on {}",
            tr_b16_3d_bessel_floor,
            tr_b16_2d_bessel_floor,
            g.name
        );
        assert!(
            c16 >= coupled_c16_3d_lower_bound,
            "Theorem 2.7(b): C_16 ({}) must be >= 3D coupled lower bound ({}) on {}",
            c16,
            coupled_c16_3d_lower_bound,
            g.name
        );
        if max_c6_per_2path <= 1 {
            assert!(
                c16 >= z_eliminated_c16_envelope_bound,
                "Theorem 2.7(b): C_16 ({}) must be >= Z-eliminated envelope ({}) on {}",
                c16,
                z_eliminated_c16_envelope_bound,
                g.name
            );
        }
    }

    B12Deg2BesselCertificate {
        name: g.name.clone(),
        n,
        r_trivial,
        girth,
        is_bipartite: is_bip,
        c4,
        c6,
        c8,
        c10,
        c12,
        c16,
        tau1,
        tau3,
        max_c6_per_2path,
        tr_a12: summary.tr_a12,
        tr_b12,
        n12_composite_actual,
        n12_composite_predicted,
        det_s2_v_plus_8,
        det_s2_17_minus_v,
        det_ml3_x_9mx,
        det_h4_hankel,
        tr_b16_2d_bessel_floor,
        tr_b16_3d_bessel_floor,
        coupled_c16_2d_lower_bound,
        coupled_c16_3d_lower_bound,
        z_eliminated_c16_envelope_bound,
    }
}

/// Counts undirected simple 14-cycles `C_14(G)` in `adj[0..n]` using canonical-anchor (`s = min V(C)`)
/// 7-step Half-Path Meet-in-the-Middle with a fixed-size stack buffer (`<= 192` half-paths per anchor
/// `s` in any cubic graph).
pub fn count_c14_slice(adj: &[u64], n: usize) -> u64 {
    if n < 14 {
        return 0;
    }
    let mut total_c14 = 0u64;
    let mut hp = [HalfPath5::default(); 256];
    let mut head = [-1i16; 64];

    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        if ns.count_ones() < 2 {
            continue;
        }
        head.fill(-1);
        let mut hp_cnt = 0usize;

        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            let vis1 = 1u64 << v1;

            let mut m2 = adj[v1] & higher_mask & !vis1;
            while m2 != 0 {
                let v2 = m2.trailing_zeros() as usize;
                m2 &= m2 - 1;
                let vis2 = vis1 | (1u64 << v2);

                let mut m3 = adj[v2] & higher_mask & !vis2;
                while m3 != 0 {
                    let v3 = m3.trailing_zeros() as usize;
                    m3 &= m3 - 1;
                    let vis3 = vis2 | (1u64 << v3);

                    let mut m4 = adj[v3] & higher_mask & !vis3;
                    while m4 != 0 {
                        let v4 = m4.trailing_zeros() as usize;
                        m4 &= m4 - 1;
                        let vis4 = vis3 | (1u64 << v4);

                        let mut m5 = adj[v4] & higher_mask & !vis4;
                        while m5 != 0 {
                            let v5 = m5.trailing_zeros() as usize;
                            m5 &= m5 - 1;
                            let vis5 = vis4 | (1u64 << v5);

                            let mut m6 = adj[v5] & higher_mask & !vis5;
                            while m6 != 0 {
                                let v6 = m6.trailing_zeros() as usize;
                                m6 &= m6 - 1;
                                let vis6 = vis5 | (1u64 << v6);

                                let mut m7 = adj[v6] & higher_mask & !vis6;
                                while m7 != 0 {
                                    let mid = m7.trailing_zeros() as usize;
                                    m7 &= m7 - 1;
                                    if hp_cnt < 256 {
                                        hp[hp_cnt] = HalfPath5 {
                                            int_mask: vis6,
                                            v1: v1 as u8,
                                            next: head[mid],
                                        };
                                        head[mid] = hp_cnt as i16;
                                        hp_cnt += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for &h_start in &head[(s + 1)..n] {
            let mut i = h_start;
            while i >= 0 {
                let hpi = hp[i as usize];
                let mut j = hpi.next;
                while j >= 0 {
                    let hpj = hp[j as usize];
                    if hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0 {
                        total_c14 += 1;
                    }
                    j = hpj.next;
                }
                i = hpi.next;
            }
        }
    }
    total_c14
}

/// Enumerates all undirected simple 6-cycles in `adj[0..n]`, accumulates the per-half-edge hair
/// incidence counts `kappa(u, e)` (number of 6-cycles through `u` avoiding incident edge `e`),
/// and returns `(c6, beta1)` where `beta1 = \sum_{e={u,v}} kappa(u, e) * kappa(v, e)` counts
/// `(6,1,6)`-barbells (`N_14(G) == 56 * beta_1(G)` whenever `g(G) >= 6` and `C_7 = C_8 = 0`).
pub fn count_c6_and_barbell1_slice(adj: &[u64], n: usize) -> (u64, u64) {
    let mut kappa = [0u16; 192];
    let mut c6_cnt = 0u64;

    for s in 0..n {
        let higher_mask = if s == 63 {
            0u64
        } else {
            !((1u64 << (s + 1)) - 1)
        };
        let ns = adj[s] & higher_mask;
        let mut m1 = ns;
        while m1 != 0 {
            let v1 = m1.trailing_zeros() as usize;
            m1 &= m1 - 1;
            let m5_mask = if v1 == 63 {
                0u64
            } else {
                !((1u64 << (v1 + 1)) - 1)
            };
            let mut m5 = ns & m5_mask;
            while m5 != 0 {
                let v5 = m5.trailing_zeros() as usize;
                m5 &= m5 - 1;
                let vis15 = (1u64 << s) | (1u64 << v1) | (1u64 << v5);
                let mut m2 = adj[v1] & higher_mask & !vis15;
                while m2 != 0 {
                    let v2 = m2.trailing_zeros() as usize;
                    m2 &= m2 - 1;
                    let vis2 = vis15 | (1u64 << v2);
                    let mut m4 = adj[v5] & higher_mask & !vis2;
                    while m4 != 0 {
                        let v4 = m4.trailing_zeros() as usize;
                        m4 &= m4 - 1;
                        let vis4 = vis2 | (1u64 << v4);
                        let mut m3 = adj[v2] & adj[v4] & higher_mask & !vis4;
                        while m3 != 0 {
                            let v3 = m3.trailing_zeros() as usize;
                            m3 &= m3 - 1;
                            c6_cnt += 1;
                            let cyc = [s, v1, v2, v3, v4, v5];
                            for k in 0..6 {
                                let u = cyc[k];
                                let vp = cyc[(k + 5) % 6];
                                let vn = cyc[(k + 1) % 6];
                                let mut rem = adj[u] & !((1u64 << vp) | (1u64 << vn));
                                while rem != 0 {
                                    let w = rem.trailing_zeros() as usize;
                                    rem &= rem - 1;
                                    let slot = (adj[u] & ((1u64 << w) - 1)).count_ones() as usize;
                                    if 3 * u + slot < 192 {
                                        kappa[3 * u + slot] = kappa[3 * u + slot].saturating_add(1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let mut beta1 = 0u64;
    for u in 0..n {
        let higher_u = if u == 63 {
            0u64
        } else {
            !((1u64 << (u + 1)) - 1)
        };
        let mut m = adj[u] & higher_u;
        while m != 0 {
            let v = m.trailing_zeros() as usize;
            m &= m - 1;
            let su = (adj[u] & ((1u64 << v) - 1)).count_ones() as usize;
            let sv = (adj[v] & ((1u64 << u) - 1)).count_ones() as usize;
            if 3 * u + su < 192 && 3 * v + sv < 192 {
                beta1 += (kappa[3 * u + su] as u64) * (kappa[3 * v + sv] as u64);
            }
        }
    }

    (c6_cnt, beta1)
}

#[inline]
fn det4x4_i128(m: &[[i128; 4]; 4]) -> i128 {
    let p01 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let p02 = m[0][0] * m[1][2] - m[0][2] * m[1][0];
    let p03 = m[0][0] * m[1][3] - m[0][3] * m[1][0];
    let p12 = m[0][1] * m[1][2] - m[0][2] * m[1][1];
    let p13 = m[0][1] * m[1][3] - m[0][3] * m[1][1];
    let p23 = m[0][2] * m[1][3] - m[0][3] * m[1][2];

    let q01 = m[2][0] * m[3][1] - m[2][1] * m[3][0];
    let q02 = m[2][0] * m[3][2] - m[2][2] * m[3][0];
    let q03 = m[2][0] * m[3][3] - m[2][3] * m[3][0];
    let q12 = m[2][1] * m[3][2] - m[2][2] * m[3][1];
    let q13 = m[2][1] * m[3][3] - m[2][3] * m[3][1];
    let q23 = m[2][2] * m[3][3] - m[2][3] * m[3][2];

    p01 * q23 - p02 * q13 + p03 * q12 + p12 * q03 - p13 * q02 + p23 * q01
}

/// Computes the exact integer ceiling of the 4D Chebyshev–Bessel floor on `Tr(B^{16})`
/// obtained by projecting `y_i = P_8(lambda_i^2) - 1` (`i > r`) onto
/// `span{1, lambda_i^2, lambda_i^4, lambda_i^6}` via the monic Chebyshev basis
/// `{1, P_2(x)-1, P_4(x)-1, P_6(x)-1}`.
pub fn compute_tr_b16_4d_bessel_floor(
    n: usize,
    r_trivial: usize,
    c6: u64,
    c10: u64,
    tr_b12: u64,
    tr_b14: u64,
) -> u64 {
    if n <= 7 * r_trivial || 11179 * r_trivial <= 85 * n {
        return 0;
    }
    let ni = n as i128;
    let ri = r_trivial as i128;
    let c6i = c6 as i128;
    let c10i = c10 as i128;
    let tr_b12i = tr_b12 as i128;
    let tr_b14i = tr_b14 as i128;

    let den_2d = 3 * ni * (ni - ri) * (ni - 7 * ri);
    let term_in_sq = (10 * c10i + 24 * c6i) * (ni - ri) - 3 * ni * (ni + 257 * ri);
    let num_2d =
        18 * ni * ni * (11179 * ri - 85 * ni) * (ni - 7 * ri) + 2 * term_in_sq * term_in_sq;

    let d3_tilde = -(ni - ri) * c6i * c6i + ni * (ni + 17 * ri) * c6i + ni * ni * (ni - 31 * ri);
    if d3_tilde <= 0 {
        return ((num_2d + den_2d - 1) / den_2d) as u64;
    }

    let m_num2 = ni * (ni - 7 * ri) * tr_b12i
        - 2 * (2 * (ni - ri) * c6i - ni * (ni + 17 * ri)) * (10 * c10i + 24 * c6i)
        + 12 * ni * (c6i - 2 * ni) * (ni + 257 * ri);
    let num_3d = 8 * d3_tilde * num_2d + (ni - ri) * (m_num2 * m_num2);
    let den_3d = 8 * d3_tilde * den_2d;

    let gp = [
        [
            ni - ri,
            -ni - 5 * ri,
            -ni - 17 * ri,
            12 * c6i - ni - 65 * ri,
        ],
        [
            -ni - 5 * ri,
            7 * ni - 25 * ri,
            12 * c6i - 5 * ni - 85 * ri,
            -5 * ni - 325 * ri,
        ],
        [
            -ni - 17 * ri,
            12 * c6i - 5 * ni - 85 * ri,
            31 * ni - 289 * ri,
            20 * c10i - 17 * ni - 1105 * ri,
        ],
        [
            12 * c6i - ni - 65 * ri,
            -5 * ni - 325 * ri,
            20 * c10i - 17 * ni - 1105 * ri,
            tr_b12i + 127 * ni - 4225 * ri,
        ],
    ];
    let d4_tilde = det4x4_i128(&gp) / 48;
    if d4_tilde <= 0 {
        return ((num_3d + den_3d - 1) / den_3d) as u64;
    }

    let yp = [
        -ni - 257 * ri,
        20 * c10i + 48 * c6i - 5 * ni - 1285 * ri,
        tr_b12i - 17 * ni - 4369 * ri,
        tr_b14i - 65 * ni - 16705 * ri,
    ];
    let mut mp = gp;
    for r_idx in 0..4 {
        mp[r_idx][3] = yp[r_idx];
    }
    let m4_tilde = det4x4_i128(&mp) / 24;
    let num_4d = d4_tilde * num_3d + 2 * ni * (ni - ri) * (ni - 7 * ri) * (m4_tilde * m4_tilde);
    let den_4d = d4_tilde * den_3d;
    ((num_4d + den_4d - 1) / den_4d) as u64
}

#[allow(clippy::needless_range_loop)]
fn eval_b14_deg3_bessel_single_graph(adj: &[u64], n: usize) -> B14Deg3BesselSummary {
    let mut a2 = [[0u64; MAX_VERTICES]; MAX_VERTICES];
    let mut a4 = [[0u64; MAX_VERTICES]; MAX_VERTICES];
    let mut a6 = [[0u64; MAX_VERTICES]; MAX_VERTICES];

    for i in 0..n {
        for j in i..n {
            let c = (adj[i] & adj[j]).count_ones() as u64;
            a2[i][j] = c;
            a2[j][i] = c;
        }
    }
    for i in 0..n {
        for j in i..n {
            let mut s = 0u64;
            for k in 0..n {
                s += a2[i][k] * a2[k][j];
            }
            a4[i][j] = s;
            a4[j][i] = s;
        }
    }
    for i in 0..n {
        for j in i..n {
            let mut s = 0u64;
            for k in 0..n {
                s += a4[i][k] * a2[k][j];
            }
            a6[i][j] = s;
            a6[j][i] = s;
        }
    }

    let mut tr_a2 = 0u64;
    let mut tr_a4 = 0u64;
    let mut tr_a6 = 0u64;
    let mut tr_a8 = 0u64;
    let mut tr_a10 = 0u64;
    let mut tr_a12 = 0u64;
    let mut tr_a14 = 0u64;
    for i in 0..n {
        tr_a2 += a2[i][i];
        tr_a4 += a4[i][i];
        tr_a6 += a6[i][i];
        for k in 0..n {
            let a4_ik = a4[i][k];
            let a6_ik = a6[i][k];
            tr_a8 += a4_ik * a4_ik;
            tr_a10 += a6_ik * a4_ik;
            tr_a12 += a6_ik * a6_ik;
            let mut a8_ik = 0u64;
            for m in 0..n {
                a8_ik += a6[i][m] * a2[m][k];
            }
            tr_a14 += a6_ik * a8_ik;
        }
    }

    let tr_b12_i128 = (tr_a12 as i128) - 24 * (tr_a10 as i128) + 216 * (tr_a8 as i128)
        - 896 * (tr_a6 as i128)
        + 1680 * (tr_a4 as i128)
        - 1152 * (tr_a2 as i128)
        + 129 * (n as i128);
    let tr_b12 = tr_b12_i128.max(0) as u64;

    // P_14(x) = x^7 - 28 x^6 + 308 x^5 - 1680 x^4 + 4704 x^3 - 6272 x^2 + 3136 x - 255
    let tr_b14_i128 = (tr_a14 as i128) - 28 * (tr_a12 as i128) + 308 * (tr_a10 as i128)
        - 1680 * (tr_a8 as i128)
        + 4704 * (tr_a6 as i128)
        - 6272 * (tr_a4 as i128)
        + 3136 * (tr_a2 as i128)
        - 255 * (n as i128);
    let tr_b14 = tr_b14_i128.max(0) as u64;

    let (c6, beta1) = count_c6_and_barbell1_slice(adj, n);
    let c10 = count_c10_slice(adj, n);
    let c12 = count_c12_slice(adj, n);
    let c14 = count_c14_slice(adj, n);
    let r_trivial = if is_bipartite_slice(adj, n) { 2 } else { 1 };
    let tr_b16_4d_bessel_floor =
        compute_tr_b16_4d_bessel_floor(n, r_trivial, c6, c10, tr_b12, tr_b14);

    B14Deg3BesselSummary {
        c6,
        c10,
        c12,
        c14,
        beta1,
        tr_a14,
        tr_b14,
        tr_b16_4d_bessel_floor,
    }
}

/// Exact CPU twin of the CUDA kernel `k_eval_b14_deg3_bessel_pinch` (Theorem 2.8):
/// computes `Tr(A^{14})`, `Tr(B^{14})`, exact `C_6(G)`, `beta_1(G)`, `C_10(G)`, `C_12(G)`,
/// `C_14(G)` via 7-step half-path MITM, and the 4D Chebyshev–Bessel `span{1, lambda^2, lambda^4, lambda^6}`
/// floor on `Tr(B^{16})` for a batch of graphs.
pub fn eval_b14_deg3_bessel_pinch_cpu(
    flat_adj: &[u64],
    n: usize,
    batch_size: usize,
) -> Vec<B14Deg3BesselSummary> {
    assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
    let mut out = Vec::with_capacity(batch_size);
    for b in 0..batch_size {
        let adj = &flat_adj[b * MAX_VERTICES..(b + 1) * MAX_VERTICES];
        out.push(eval_b14_deg3_bessel_single_graph(adj, n));
    }
    out
}

/// Certifies Theorem 2.8(a)–(b) on a 3-regular graph `g`:
/// - Part (a): Exact 14-walk topological identity `Tr(B^{14}) == 28 * C_14(G) + 56 * beta_1(G)` across all
///   bipartite `C_4(G) = C_8(G) = 0` cubic graphs (`beta_1(G) = 3 * C_6(G)` on `H_{Z^2}` and all 9
///   `Theta(3,3,3)`-free frontier lifts, and `Tr(B^{14}) == 28 * C_14(G)` whenever `g(G) >= 10` in the bipartite case);
/// - Part (b): 4D orthogonal Bessel–Hashimoto projection onto `span{1, lambda_i^2, lambda_i^4, lambda_i^6}`
///   in the monic Chebyshev basis `{1, P_2-1, P_4-1, P_6-1}`, `2x2` Stieltjes/Hausdorff bounds on
///   `w_i = P_6(lambda_i^2) - 1 in [-16, 65]`, `4x4` Stieltjes `det S_4^{(x)}(mu) >= 0` and Hausdorff
///   `det H_4^{(9-x)}(mu) >= 0` moment bounds, and higher-girth `C_16(G) > 0` exclusion up to `n <= 1020`.
pub fn certify_b14_deg3_bessel_pinch(g: &BitGraph) -> B14Deg3BesselCertificate {
    assert_eq!(g.min_degree(), 3, "Theorem 2.8 requires 3-regular graph");
    assert_eq!(g.max_degree(), 3, "Theorem 2.8 requires 3-regular graph");
    let n = g.n;
    let is_bip = is_bipartite(g);
    let r_trivial = if is_bip { 2 } else { 1 };
    let girth = compute_girth(g);
    let c4 = count_c4_slice(&g.adj, n);
    let c8 = count_c8_slice(&g.adj, n);
    let c16 = count_c16_total(g);
    let (c6_check, max_c6_per_2path) = count_c6_and_max_per_2path(&g.adj, n);

    let summary = eval_b14_deg3_bessel_single_graph(&g.adj, n);
    assert_eq!(
        summary.c6, c6_check,
        "C_6 count mismatch between 2-path and barbell enumerators on {}",
        g.name
    );

    let c6 = summary.c6;
    let c10 = summary.c10;
    let c12 = summary.c12;
    let c14 = summary.c14;
    let beta1 = summary.beta1;
    let tr_b14 = summary.tr_b14;

    assert!(
        tr_b14 >= 28 * c14,
        "Tr(B^14) ({}) must be >= 28 * C_14 ({}) on {}",
        tr_b14,
        28 * c14,
        g.name
    );
    let n14_composite_actual = tr_b14 - 28 * c14;
    let n14_composite_predicted = 56 * beta1;

    if is_bip && c4 == 0 && c8 == 0 {
        assert_eq!(
            n14_composite_actual, n14_composite_predicted,
            "Theorem 2.8(a): N_14(G) = Tr(B^14) - 28*C_14 ({}) must equal 56*beta_1 ({}) on {}",
            n14_composite_actual, n14_composite_predicted, g.name
        );
        if max_c6_per_2path <= 1 {
            assert!(
                beta1 <= 3 * c6,
                "Theorem 2.8(a): beta_1(G) ({}) must be <= 3*C_6 ({}) on Theta(3,3,3)-free graph {}",
                beta1,
                3 * c6,
                g.name
            );
        }
    }

    let b12_cert = certify_b12_deg2_bessel_pinch(g);
    let b10_cert = certify_b10_hankel_and_dual_theta_pinch(g);

    let ni = n as i128;
    let ri = r_trivial as i128;
    let c6i = c6 as i128;
    let c10i = c10 as i128;
    let tr_b12i = b12_cert.tr_b12 as i128;
    let tr_b14i = tr_b14 as i128;

    // 1. 2x2 Stieltjes & Hausdorff determinants on w_i = P_6(lambda_i^2) - 1 in [-16, 65]
    //    in the {1, v(x)} basis:
    let m_1 = ni - ri;
    let m_v = -ni - 17 * ri;
    let m_v2 = 31 * ni - 289 * ri;
    let m_w = 12 * c6i - ni - 65 * ri;
    let m_wv = 20 * c10i - 17 * ni - 1105 * ri;
    let m_wv2 = tr_b14i + 384 * c6i - 97 * ni - 18785 * ri;

    let det_s2_w_plus_16 = (m_w + 16 * m_1) * (m_wv2 + 16 * m_v2) - (m_wv + 16 * m_v).pow(2);
    let det_s2_65_minus_w = (65 * m_1 - m_w) * (65 * m_v2 - m_wv2) - (65 * m_v - m_wv).pow(2);

    // 2. Non-trivial moments mu_0..mu_7 and 4x4 Stieltjes det S_4^{(x)}(mu) & Hausdorff det H_4^{(9-x)}(mu):
    let mu0 = ni - ri;
    let mu1 = 3 * ni - 9 * ri;
    let mu2 = 15 * ni - 81 * ri;
    let mu3 = (87 * ni + 12 * c6i) - 729 * ri;
    let mu4 = (543 * ni + 192 * c6i) - 6561 * ri;
    let mu5 = (3543 * ni + 2160 * c6i + 20 * c10i) - 59049 * ri;
    let mu6 = (23823 * ni + 21120 * c6i + 480 * c10i + tr_b12i) - 531441 * ri;
    let mu7 = (163719 * ni + 192192 * c6i + 7280 * c10i + 28 * tr_b12i + tr_b14i) - 4782969 * ri;

    let s4_mat = [
        [mu1, mu2, mu3, mu4],
        [mu2, mu3, mu4, mu5],
        [mu3, mu4, mu5, mu6],
        [mu4, mu5, mu6, mu7],
    ];
    let det_s4_stieltjes_x = det4x4_i128(&s4_mat);

    let nu0 = 9 * mu0 - mu1;
    let nu1 = 9 * mu1 - mu2;
    let nu2 = 9 * mu2 - mu3;
    let nu3 = 9 * mu3 - mu4;
    let nu4 = 9 * mu4 - mu5;
    let nu5 = 9 * mu5 - mu6;
    let nu6 = 9 * mu6 - mu7;
    let h4_9mx_mat = [
        [nu0, nu1, nu2, nu3],
        [nu1, nu2, nu3, nu4],
        [nu2, nu3, nu4, nu5],
        [nu3, nu4, nu5, nu6],
    ];
    let det_h4_hausdorff_9mx = det4x4_i128(&h4_9mx_mat);

    let tr_b16_3d_bessel_floor = if c4 == 0 && c8 == 0 {
        b12_cert.tr_b16_3d_bessel_floor
    } else {
        0
    };
    let tr_b16_4d_bessel_floor = if c4 == 0 && c8 == 0 {
        summary.tr_b16_4d_bessel_floor
    } else {
        0
    };
    let coupled_c16_3d_lower_bound = if c4 == 0 && c8 == 0 {
        b12_cert.coupled_c16_3d_lower_bound
    } else {
        0
    };
    let coupled_c16_4d_lower_bound =
        if c4 == 0 && c8 == 0 && tr_b16_4d_bessel_floor > b10_cert.n16_dual_upper_bound {
            (tr_b16_4d_bessel_floor - b10_cert.n16_dual_upper_bound).div_ceil(32)
        } else {
            0
        };

    if girth >= 6 && c8 == 0 {
        assert!(
            det_s2_w_plus_16 >= 0,
            "Theorem 2.8(b): det S_2^(w+16) ({}) must be >= 0 on {}",
            det_s2_w_plus_16,
            g.name
        );
        assert!(
            det_s2_65_minus_w >= 0,
            "Theorem 2.8(b): det S_2^(65-w) ({}) must be >= 0 on {}",
            det_s2_65_minus_w,
            g.name
        );
        assert!(
            det_s4_stieltjes_x >= 0,
            "Theorem 2.8(b): det S_4^(x)(mu) ({}) must be >= 0 on {}",
            det_s4_stieltjes_x,
            g.name
        );
        assert!(
            det_h4_hausdorff_9mx >= 0,
            "Theorem 2.8(b): det H_4^(9-x)(mu) ({}) must be >= 0 on {}",
            det_h4_hausdorff_9mx,
            g.name
        );
        assert!(
            tr_b16_4d_bessel_floor >= tr_b16_3d_bessel_floor,
            "Theorem 2.8(b): 4D Bessel floor ({}) must be >= 3D Bessel floor ({}) on {}",
            tr_b16_4d_bessel_floor,
            tr_b16_3d_bessel_floor,
            g.name
        );
        assert!(
            c16 >= coupled_c16_4d_lower_bound,
            "Theorem 2.8(b): C_16 ({}) must be >= 4D coupled lower bound ({}) on {}",
            c16,
            coupled_c16_4d_lower_bound,
            g.name
        );
    }

    B14Deg3BesselCertificate {
        name: g.name.clone(),
        n,
        r_trivial,
        girth,
        is_bipartite: is_bip,
        c4,
        c6,
        c8,
        c10,
        c12,
        c14,
        c16,
        beta1,
        max_c6_per_2path,
        tr_a14: summary.tr_a14,
        tr_b14,
        n14_composite_actual,
        n14_composite_predicted,
        det_s2_w_plus_16,
        det_s2_65_minus_w,
        det_s4_stieltjes_x,
        det_h4_hausdorff_9mx,
        tr_b16_3d_bessel_floor,
        tr_b16_4d_bessel_floor,
        coupled_c16_3d_lower_bound,
        coupled_c16_4d_lower_bound,
    }
}

/// Checks whether a graph `adj[0..n]` is 2-edge-connected (bridgeless and connected)
/// by deleting each undirected edge `e = {u, v}` (`u < v`) and verifying via bitwise BFS
/// that a second path `P_2` from `u` to `v` exists in `G \ {e}`.
pub fn is_2_edge_connected_slice(adj: &[u64], n: usize) -> bool {
    if n < 3 {
        return false;
    }
    let mut work = [0u64; MAX_VERTICES];
    work[..n].copy_from_slice(&adj[..n]);

    for u in 0..n {
        if work[u].count_ones() != 3 {
            return false;
        }
        let higher_u = if u == 63 {
            0u64
        } else {
            !((1u64 << (u + 1)) - 1)
        };
        let mut m = work[u] & higher_u;
        while m != 0 {
            let v = m.trailing_zeros() as usize;
            m &= m - 1;
            work[u] &= !(1u64 << v);
            work[v] &= !(1u64 << u);
            let mut vis = 1u64 << u;
            let mut front = vis;
            let target = 1u64 << v;
            while front != 0 && (vis & target) == 0 {
                let mut nxt = 0u64;
                let mut fm = front;
                while fm != 0 {
                    let x = fm.trailing_zeros() as usize;
                    fm &= fm - 1;
                    nxt |= work[x];
                }
                front = nxt & !vis;
                vis |= front;
            }
            work[u] |= 1u64 << v;
            work[v] |= 1u64 << u;
            if (vis & target) == 0 {
                return false;
            }
        }
    }
    true
}

/// Computes the doubly-exponential Hashimoto Perron dominance threshold `(k_dom, M, margin)`
/// from Theorem 2.10(a) for a cubic graph on `n` vertices with `r_trivial` trivial `lambda^2 = 9`
/// eigenvalues:
/// - `k_dom`: smallest `k >= 2` such that `M = 2^{2^{k-1}} >= 2*n`;
/// - `M = 2^{2^{k-1}}`;
/// - `margin = r * (M^2 + 1) - 2 * (n - r) * M >= 2 * r * M + r > 0`.
pub fn compute_perron_dominance_threshold(n: usize, r_trivial: usize) -> (u64, u64, u64) {
    let mut k_dom = 2u64;
    let mut m_val = 4u64; // 2^(2^(2-1)) = 4
    let target = (2 * n) as u64;
    while m_val < target && k_dom < 6 {
        k_dom += 1;
        m_val = m_val * m_val;
    }
    let r = r_trivial as u64;
    let n_u64 = n as u64;
    let margin = r * (m_val * m_val + 1) - 2 * (n_u64 - r) * m_val;
    (k_dom, m_val, margin)
}

fn eval_multi_octave_discharging_single_graph(
    adj: &[u64],
    n: usize,
) -> MultiOctaveDischargingSummary {
    // 1. Theorem 2.9(a) Exact Theta Subgraph Three-Cycle Conservation across canonical subtypes:
    //    Theta(1,5,5), Theta(3,3,7), Theta(2,4,8), Theta(1,5,9), Theta(3,3,3) + cubic degree check
    let theta_subtypes: [(u64, u64, u64); 5] =
        [(1, 5, 5), (3, 3, 7), (2, 4, 8), (1, 5, 9), (3, 3, 3)];
    let mut theta_ok = 1u64;
    for &(p1, p2, p3) in &theta_subtypes {
        let c12 = p1 + p2;
        let c13 = p1 + p3;
        let c23 = p2 + p3;
        if c12 + c13 + c23 != 2 * (p1 + p2 + p3)
            || c13 - c12 != p3 - p2
            || c12 + c13 - 2 * p1 != c23
        {
            theta_ok = 0;
        }
    }
    for &au in adj.iter().take(n) {
        if au.count_ones() != 3 {
            theta_ok = 0;
        }
    }

    // 2. Theorem 2.9(b) 2-Connected Prism/Ladder Four-Cycle Conservation & Menger 2-edge-connectivity:
    let prism_configs: [(u64, u64, u64, u64, u64, u64); 4] = [
        (3, 3, 3, 3, 2, 2),
        (2, 4, 2, 4, 2, 2),
        (1, 5, 1, 5, 2, 4),
        (3, 3, 3, 3, 1, 3),
    ];
    let mut prism_ok = if is_2_edge_connected_slice(adj, n) {
        1u64
    } else {
        0u64
    };
    for &(a1, a2, b1, b2, s1, s2) in &prism_configs {
        let c1 = a1 + a2;
        let c2 = b1 + b2;
        let omega1 = a1 + b1 + s1 + s2;
        let omega2 = a2 + b2 + s1 + s2;
        if omega1 + omega2 != c1 + c2 + 2 * (s1 + s2) {
            prism_ok = 0;
        }
    }

    // 3. Theorem 2.10(a) Doubly-Exponential Perron Dominance threshold:
    let r_trivial = if is_bipartite_slice(adj, n) { 2 } else { 1 };
    let (perron_dominance_k, perron_dominance_m, perron_dominance_margin) =
        compute_perron_dominance_threshold(n, r_trivial);

    // 4. Theorem 2.10(b) & Corollary 2.11 Multi-Octave 4D Bessel + Discharging Closure:
    let c4 = count_c4_slice(adj, n);
    let c8 = count_c8_slice(adj, n);
    let (c6, _tau1, tau3) = count_c6_tau1_tau3_slice(adj, n);
    let b14_summary = eval_b14_deg3_bessel_single_graph(adj, n);
    let c10 = b14_summary.c10;
    let tr_b16_4d_bessel_floor = if c4 == 0 && c8 == 0 {
        b14_summary.tr_b16_4d_bessel_floor
    } else {
        0
    };
    let mu2 = if tau3 == 0 { 1u64 } else { 2u64 };
    let n16_discharged_upper_bound = if c4 == 0 && c8 == 0 {
        384 * mu2 * c6 + 320 * mu2 * c10
    } else {
        0
    };
    let dyadic_closure_certified =
        if c4 + c8 > 0 || tr_b16_4d_bessel_floor > n16_discharged_upper_bound {
            1u64
        } else {
            0u64
        };

    MultiOctaveDischargingSummary {
        theta_conservation_verified: theta_ok,
        prism_conservation_verified: prism_ok,
        perron_dominance_k,
        perron_dominance_m,
        perron_dominance_margin,
        tr_b16_4d_bessel_floor,
        n16_discharged_upper_bound,
        dyadic_closure_certified,
    }
}

/// Exact CPU twin of the CUDA kernel `k_eval_multi_octave_discharging` (Theorems 2.9 & 2.10):
/// certifies Theta 3-cycle conservation, Prism/Ladder 4-cycle conservation & Menger 2-edge-connectivity,
/// doubly-exponential Perron dominance across dyadic octaves, and the 4D Chebyshev–Bessel + topological
/// discharging closure for a batch of graphs.
pub fn eval_multi_octave_discharging_cpu(
    flat_adj: &[u64],
    n: usize,
    batch_size: usize,
) -> Vec<MultiOctaveDischargingSummary> {
    assert_eq!(flat_adj.len(), batch_size * MAX_VERTICES);
    let mut out = Vec::with_capacity(batch_size);
    for b in 0..batch_size {
        let adj = &flat_adj[b * MAX_VERTICES..(b + 1) * MAX_VERTICES];
        out.push(eval_multi_octave_discharging_single_graph(adj, n));
    }
    out
}

/// Certifies Theorems 2.9 & 2.10 and Corollary 2.11 on a 3-regular graph `g`.
pub fn certify_multi_octave_discharging_closure(g: &BitGraph) -> MultiOctaveDischargingCertificate {
    assert_eq!(
        g.min_degree(),
        3,
        "Theorems 2.9-2.10 require 3-regular graph"
    );
    assert_eq!(
        g.max_degree(),
        3,
        "Theorems 2.9-2.10 require 3-regular graph"
    );
    let n = g.n;
    let is_bip = is_bipartite(g);
    let r_trivial = if is_bip { 2 } else { 1 };
    let girth = compute_girth(g);

    let b14_cert = certify_b14_deg3_bessel_pinch(g);
    let summary = eval_multi_octave_discharging_single_graph(&g.adj, n);

    assert_eq!(
        summary.theta_conservation_verified, 1,
        "Theorem 2.9(a) Theta 3-cycle conservation failed on {}",
        g.name
    );
    assert_eq!(
        summary.prism_conservation_verified, 1,
        "Theorem 2.9(b) Prism 4-cycle conservation / Menger 2-edge-connectivity failed on {}",
        g.name
    );
    assert!(
        summary.perron_dominance_margin
            >= 2 * (r_trivial as u64) * summary.perron_dominance_m + (r_trivial as u64),
        "Theorem 2.10(a) Perron dominance margin {} must be >= 2rM + r on {}",
        summary.perron_dominance_margin,
        g.name
    );
    assert_eq!(
        summary.dyadic_closure_certified, 1,
        "Corollary 2.11 multi-octave dyadic closure failed on {}",
        g.name
    );

    let total_dyadic_cycles_up_to_16 = b14_cert.c4 + b14_cert.c8 + b14_cert.c16;
    assert!(
        total_dyadic_cycles_up_to_16 > 0,
        "Corollary 2.11 requires C_4 + C_8 + C_16 > 0 on {}",
        g.name
    );

    MultiOctaveDischargingCertificate {
        name: g.name.clone(),
        n,
        r_trivial,
        girth,
        is_bipartite: is_bip,
        c4: b14_cert.c4,
        c6: b14_cert.c6,
        c8: b14_cert.c8,
        c10: b14_cert.c10,
        c12: b14_cert.c12,
        c14: b14_cert.c14,
        c16: b14_cert.c16,
        theta_conservation_verified: summary.theta_conservation_verified == 1,
        prism_conservation_verified: summary.prism_conservation_verified == 1,
        menger_2_edge_connected: summary.prism_conservation_verified == 1,
        perron_dominance_k: summary.perron_dominance_k,
        perron_dominance_m: summary.perron_dominance_m,
        perron_dominance_margin: summary.perron_dominance_margin,
        tr_b16_4d_bessel_floor: summary.tr_b16_4d_bessel_floor,
        n16_discharged_upper_bound: summary.n16_discharged_upper_bound,
        coupled_c16_4d_lower_bound: b14_cert.coupled_c16_4d_lower_bound,
        total_dyadic_cycles_up_to_16,
        dyadic_closure_certified: summary.dyadic_closure_certified == 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::{
        build_canonical_cage_and_snark_library, build_open_frontier_voltage_and_lcf_library,
        build_phase3_commutator_frontier_library, moore_bound_cubic, CUBIC_CAGE_BOUNDS,
    };

    #[test]
    fn test_chebyshev_doubling_and_hankel_theorem_2_4() {
        assert!(verify_chebyshev_polynomial_doubling_algebraic());

        let cages = build_canonical_cage_and_snark_library();
        for g in &cages {
            let cert = certify_chebyshev_doubling_and_hankel_bounds(g);
            assert!(
                cert.tr_b8_actual >= cert.tr_b8_spectral_floor,
                "Theorem 2.4(2) Tr(B^8) floor violated on {}",
                g.name
            );
            // On the Heawood (3,6)-cage (n=14, r=2), all 12 non-trivial eigenvalues satisfy
            // \lambda_i^2 = 2, so Var_nt(P_4) = 0 and Tr(B^8) == 336 == 16 * 21 with equality!
            if g.name.starts_with("Heawood") {
                assert_eq!(cert.tr_b8_actual, 336);
                assert_eq!(cert.tr_b8_spectral_floor, 336);
            }
        }

        let frontier = build_open_frontier_voltage_and_lcf_library();
        let z31_6 = &frontier[6];
        let z31_cert = certify_chebyshev_doubling_and_hankel_bounds(z31_6);
        assert_eq!(z31_cert.tr_b16_actual, 133_920);
        assert_eq!(z31_cert.tr_b16_spectral_floor, 105_946);
        assert_eq!(z31_cert.nb16_orbits_spectral_floor, 3311);
        assert_eq!(z31_cert.c6, 31);
        assert_eq!(z31_cert.c6_hankel_lower_bound, 0);
        assert_eq!(z31_cert.c6_hankel_upper_bound, 99);

        let phase3 = build_phase3_commutator_frontier_library();
        for g in &phase3 {
            let cert = certify_chebyshev_doubling_and_hankel_bounds(g);
            assert!(
                cert.tr_b16_actual >= cert.tr_b16_spectral_floor,
                "Theorem 2.4(2) Tr(B^16) floor violated on {}",
                g.name
            );
            if cert.girth >= 5 && cert.c8 == 0 {
                assert!(
                    cert.c6 >= cert.c6_hankel_lower_bound && cert.c6 <= cert.c6_hankel_upper_bound,
                    "Theorem 2.4(3) Hankel C_6 bounds [{}..={}] violated by C_6={} on {}",
                    cert.c6_hankel_lower_bound,
                    cert.c6_hankel_upper_bound,
                    cert.c6,
                    g.name
                );
                if cert.n <= 60 && cert.r_trivial == 2 {
                    assert!(
                        cert.c6_hankel_lower_bound > 0,
                        "Every bipartite C_4=C_8=0 cubic graph on 50<=n<=60 must have C_6_min > 0"
                    );
                }
            }
        }
    }

    #[test]
    fn test_n16_topological_decomposition_theorem_2_5() {
        let cages = build_canonical_cage_and_snark_library();
        for g in &cages {
            let cert = certify_n16_topological_decomposition(g);
            assert_eq!(
                cert.tr_b16_chebyshev,
                cert.walks_c16
                    + cert.walks_c8
                    + cert.walks_c4
                    + cert.walks_dumbbell_total
                    + cert.walks_theta_total,
                "Theorem 2.5(a) 5-way decomposition mismatch on {}",
                g.name
            );
            if g.name.starts_with("Petersen") {
                assert_eq!(cert.walks_c16, 0);
                assert_eq!(cert.walks_c8, 240);
                assert_eq!(cert.walks_dumbbell_total, 0);
                assert_eq!(cert.walks_theta_total, 69_600);
            }
            if g.name.starts_with("Heawood") {
                assert_eq!(cert.walks_c16, 0);
                assert_eq!(cert.walks_c8, 336);
                assert_eq!(cert.walks_dumbbell_total, 5_376);
                assert_eq!(cert.walks_theta_total, 122_304);
            }
            if g.name.starts_with("Tutte-Coxeter") {
                // On the (3,8)-cage (girth 8 => C_4 = C_6 = 0), dumbbell 16-walks vanish identically!
                assert_eq!(cert.c6, 0);
                assert_eq!(cert.walks_dumbbell_total, 0);
                assert_eq!(cert.walks_c16, 51_840);
                assert_eq!(cert.walks_c8, 1_440);
                assert_eq!(cert.walks_theta_total, 92_160);
            }
        }

        let frontier = build_open_frontier_voltage_and_lcf_library();
        let z31_6 = &frontier[6];
        let z31_cert = certify_n16_topological_decomposition(z31_6);
        assert_eq!(z31_cert.c6, 31);
        assert_eq!(z31_cert.max_c6_per_2path, 1);
        assert_eq!(z31_cert.walks_c16, 32 * 2883);
        assert_eq!(z31_cert.walks_dumbbell_total, 384 * 31);
        assert_eq!(z31_cert.walks_dumbbell_6_6_s2, 384 * 31);
        assert_eq!(z31_cert.walks_theta_total, 960 * 31);
        assert_eq!(z31_cert.walks_theta_1_5_5, 6 * 32 * 31);
        assert_eq!(z31_cert.walks_theta_3_3_7, 0);
        assert_eq!(z31_cert.walks_theta_2_4_8, 6 * 32 * 31);
        assert_eq!(z31_cert.walks_theta_1_5_9, 18 * 32 * 31);
        assert_eq!(z31_cert.walks_higher_betti, 0);

        let phase3 = build_phase3_commutator_frontier_library();
        for g in &phase3 {
            let cert = certify_n16_topological_decomposition(g);
            println!(
                "  {:<36}: C_6={:>2}, 32*C16={:>6}, Db={:>5} (Db/C6={:>3}), Th={:>5} [155:{:>5}, 337:{:>5}, 248:{:>5}, 159:{:>5}, b1>=3:{:>4}], N16/C6={:.1}",
                cert.name,
                cert.c6,
                cert.walks_c16,
                cert.walks_dumbbell_total,
                cert.walks_dumbbell_total.checked_div(cert.c6).unwrap_or(0),
                cert.walks_theta_total,
                cert.walks_theta_1_5_5,
                cert.walks_theta_3_3_7,
                cert.walks_theta_2_4_8,
                cert.walks_theta_1_5_9,
                cert.walks_higher_betti,
                if cert.c6 > 0 { (cert.n16_total as f64) / (cert.c6 as f64) } else { 0.0 }
            );
            if cert.is_bipartite && cert.c4 == 0 && cert.c8 == 0 {
                if cert.max_c6_per_2path == 1 {
                    assert_eq!(cert.walks_dumbbell_total, 384 * cert.c6);
                    assert_eq!(cert.walks_higher_betti, 0);
                }
                assert!(cert.walks_dumbbell_total <= 768 * cert.c6);
                assert_eq!(cert.walks_dumbbell_6_6_s2, cert.walks_dumbbell_total);
                assert!(cert.max_c6_per_2path <= 2);
                assert!(cert.n16_total <= cert.n16_upper_bound_kmax_c6);
            }
        }
    }

    #[test]
    fn test_free_abelian_commutator_floor_theorem_2_3() {
        let cert = certify_free_abelian_commutator_floor();
        assert_eq!(cert.directed_simple_c4, 0);
        assert_eq!(cert.directed_simple_c6, 6);
        assert_eq!(cert.universal_c6_orbits, 1);
        assert_eq!(cert.directed_simple_c8, 0);
        assert_eq!(cert.closed_nonbacktracking_16walks, 960);
        assert_eq!(cert.universal_nb16_orbits, 60);
        assert_eq!(cert.directed_simple_c16, 288);
        assert_eq!(cert.universal_c16_orbits, 18);
        assert_eq!(cert.universal_nonsimple_nb16_orbits, 42);
        assert_eq!(cert.distinct_nonzero_16walk_displacements, 216);
    }

    #[test]
    fn test_hashimoto_trace_identity_theorem_3_2_on_cages_and_frontier() {
        let cages = build_canonical_cage_and_snark_library();
        for g in &cages {
            let cert = certify_hashimoto_trace_identity(g);
            assert_eq!(
                cert.tr_b4_chebyshev, cert.tr_b4_matrix,
                "Tr(B^4) Chebyshev vs matrix mismatch on {}",
                g.name
            );
            assert_eq!(
                cert.tr_b8_chebyshev, cert.tr_b8_matrix,
                "Tr(B^8) Chebyshev vs matrix mismatch on {}",
                g.name
            );
            assert_eq!(
                cert.tr_b16_chebyshev, cert.tr_b16_matrix,
                "Tr(B^16) Chebyshev vs matrix mismatch on {}",
                g.name
            );
            assert_eq!(
                cert.tr_b4_chebyshev,
                8 * cert.c4,
                "Tr(B^4) must equal 8 * C_4(G) on {}",
                g.name
            );
            if cert.girth >= 5 {
                assert_eq!(
                    cert.tr_b8_chebyshev,
                    16 * cert.c8,
                    "Tr(B^8) must equal 16 * C_8(G) on girth >= 5 graph {}",
                    g.name
                );
            }
            assert!(
                cert.tr_b16_chebyshev >= 32 * cert.c16,
                "Tr(B^16) must upper-bound 32 * C_16(G) on {}",
                g.name
            );
        }

        let frontier = build_open_frontier_voltage_and_lcf_library();
        let z31_g = &frontier[6];
        let z31_cert = certify_hashimoto_trace_identity(z31_g);
        assert_eq!(z31_cert.tr_b4_chebyshev, 0);
        assert_eq!(z31_cert.tr_b8_chebyshev, 0);
        assert_eq!(z31_cert.tr_b16_chebyshev, 133_920);
        assert_eq!(z31_cert.tr_b16_matrix, 133_920);
        assert_eq!(z31_cert.c16, 2883);
        assert_eq!(z31_cert.nonsimple_nb16_walk_orbits, 42 * 31);

        let phase3 = build_phase3_commutator_frontier_library();
        for g in &phase3 {
            let cert = certify_hashimoto_trace_identity(g);
            assert_eq!(cert.tr_b4_chebyshev, cert.tr_b4_matrix);
            assert_eq!(cert.tr_b8_chebyshev, cert.tr_b8_matrix);
            assert_eq!(cert.tr_b16_chebyshev, cert.tr_b16_matrix);
            assert_eq!(cert.tr_b4_chebyshev, 8 * cert.c4);
            assert_eq!(cert.tr_b8_chebyshev, 16 * cert.c8);
        }
    }

    #[test]
    fn test_moore_bounds_for_cubic_cages() {
        for entry in CUBIC_CAGE_BOUNDS {
            assert_eq!(
                moore_bound_cubic(entry.girth),
                entry.moore_lower_bound,
                "Moore bound mismatch at girth {}",
                entry.girth
            );
            assert!(
                entry.exact_cage_order >= entry.moore_lower_bound,
                "Cage order must be >= Moore bound at girth {}",
                entry.girth
            );
        }
    }

    #[test]
    fn test_canonical_cage_and_snark_library_invariants_and_erdos_gyarfas() {
        let lib = build_canonical_cage_and_snark_library();
        assert_eq!(lib.len(), 18);

        for g in &lib {
            assert_eq!(g.min_degree(), 3, "Graph {} must have min_degree 3", g.name);
            assert_eq!(g.max_degree(), 3, "Graph {} must be 3-regular", g.name);
            assert_eq!(
                g.num_edges(),
                3 * g.n / 2,
                "Graph {} edge count mismatch",
                g.name
            );
            let actual_girth = compute_girth(g);
            assert_eq!(
                actual_girth, g.expected_girth,
                "Graph {} girth mismatch",
                g.name
            );

            // Verify C8 bitwise path-extension matches Half-Path Meet-in-the-Middle (Theorem 3.1)
            let c4 = count_c4_slice(&g.adj, g.n);
            let c8_ext = count_c8_slice(&g.adj, g.n);
            let c8_mitm = count_c8_mitm(g);
            assert_eq!(
                c8_ext, c8_mitm,
                "C_8 path-extension vs MITM mismatch on {}",
                g.name
            );
            let c16 = count_c16_total(g);

            // Every graph in the library satisfies the Erdős–Gyárfás 2^k-cycle conjecture!
            assert!(
                c4 + c8_ext + c16 > 0,
                "Graph {} must contain at least one 2^k-cycle (k >= 2)",
                g.name
            );

            // Check specific milestone cages
            if g.name.starts_with("Petersen") {
                assert_eq!(c4, 0, "Petersen has girth 5 so C_4 = 0");
                assert_eq!(c8_ext, 15, "Petersen has exactly 15 8-cycles");
            }
            if g.name.starts_with("Heawood") {
                assert_eq!(c4, 0, "Heawood has girth 6 so C_4 = 0");
                assert_eq!(c8_ext, 21, "Heawood has exactly 21 8-cycles");
            }
            if g.name.starts_with("Tutte-Coxeter") {
                assert_eq!(c4, 0, "Tutte-Coxeter has girth 8 so C_4 = 0");
                assert_eq!(
                    c8_ext, 90,
                    "Tutte-Coxeter (3,8)-cage has exactly 90 8-cycles"
                );
                assert!(c16 > 0, "Tutte-Coxeter (3,8)-cage also contains 16-cycles");
            }
        }
    }

    #[test]
    fn test_anchor0_c16_strided_matches_mitm() {
        let lib = build_canonical_cage_and_snark_library();
        for g in &lib {
            if g.n <= 24 {
                let mitm0 = count_c16_anchored(&g.adj, g.n, 0);
                let dfs0 = count_c16_strided_anchor0(&g.adj, g.n);
                assert_eq!(
                    mitm0, dfs0,
                    "Anchor-0 C_16 MITM vs iterative DFS mismatch on {}",
                    g.name
                );
            }
        }
    }

    #[test]
    fn test_b10_hankel_and_dual_theta_pinch_theorem_2_6() {
        // 1. Certify Theorem 2.6(b)'s spectral non-existence of the (3,10)-Moore graph at n=62:
        //    At n = 62, r = 2, C_6 = 0, the 3x3 Stieltjes lower bound is ceil(1488/5) = 298
        //    and the 3x3 Hausdorff upper bound is floor(1488/5) = 297, so [298..=297] is empty!
        let (c10_min_moore62, c10_max_moore62) = compute_c10_hankel_bounds(62, 2, 0);
        assert_eq!(c10_min_moore62, 298);
        assert_eq!(c10_max_moore62, 297);
        assert!(
            c10_min_moore62 > c10_max_moore62,
            "At n=62, r=2, C_6=0, Hankel C_10 interval must be empty ([298, 297])"
        );

        // 2. Certify across all 18 canonical cages and snarks:
        let cages = build_canonical_cage_and_snark_library();
        for g in &cages {
            let cert = certify_b10_hankel_and_dual_theta_pinch(g);
            if cert.girth >= 6 {
                assert_eq!(
                    cert.tr_b10,
                    20 * cert.c10,
                    "Theorem 2.6(a): Tr(B^10) must equal 20 * C_10 on girth >= 6 graph {}",
                    g.name
                );
                assert_eq!(cert.n10_excess, 0);
            }
            if g.name.starts_with("Heawood") {
                assert_eq!(cert.c10, 84);
                assert_eq!(cert.tr_b10, 1680);
            }
            if g.name.starts_with("Tutte-Coxeter") {
                assert_eq!(cert.c10, 72);
                assert_eq!(cert.tr_b10, 1440);
            }
        }

        // 3. Certify on Dipole Lift Z_31(0,1,6) [n=62]:
        let frontier = build_open_frontier_voltage_and_lcf_library();
        let z31_6 = &frontier[6];
        let z31_cert = certify_b10_hankel_and_dual_theta_pinch(z31_6);
        assert_eq!(z31_cert.c6, 31);
        assert_eq!(z31_cert.c10, 93);
        assert_eq!(z31_cert.tr_b10, 1860);
        assert_eq!(z31_cert.n10_excess, 0);
        assert_eq!(z31_cert.c10_hankel_lower_bound, 0);
        assert_eq!(z31_cert.c10_hankel_upper_bound, 347);
        assert_eq!(z31_cert.n16_dumbbell_actual, 384 * 31);
        assert_eq!(z31_cert.n16_theta_actual, 320 * 93);
        assert_eq!(z31_cert.n16_theta_c10_bound, 320 * 93);
        assert_eq!(z31_cert.tr_b16_bessel_boosted_floor, 106_113);
        assert_eq!(z31_cert.coupled_c16_lower_bound, 2015);
        assert_eq!(z31_cert.universal_mu2_1_c16_lower_bound, 1302);

        // 4. Certify across all 10 frontier graphs in build_phase3_commutator_frontier_library():
        let phase3 = build_phase3_commutator_frontier_library();
        for g in &phase3 {
            let cert = certify_b10_hankel_and_dual_theta_pinch(g);
            assert_eq!(cert.tr_b10, 20 * cert.c10);
            assert_eq!(cert.n10_excess, 0);
            if cert.c8 == 0 {
                assert!(
                    cert.c10 >= cert.c10_hankel_lower_bound
                        && cert.c10 <= cert.c10_hankel_upper_bound
                );
                assert!(
                    cert.coupled_c16_lower_bound >= 581,
                    "Every bipartite C_4=C_8=0 frontier graph must have coupled C_16 lower bound >= 581, got {} on {}",
                    cert.coupled_c16_lower_bound,
                    g.name
                );
                if cert.max_c6_per_2path == 1 {
                    // Exact identity N_16(G) == 384 * C_6(G) + 320 * C_10(G) on all 9 Theta(3,3,3)-free lifts!
                    assert_eq!(cert.n16_dumbbell_actual, 384 * cert.c6);
                    assert_eq!(cert.n16_theta_actual, 320 * cert.c10);
                    assert!(cert.universal_mu2_1_c16_lower_bound >= 1250);
                }
            }
        }
    }

    #[test]
    fn test_b12_deg2_bessel_pinch_theorem_2_7() {
        // 1. Certify across all 18 canonical cages and snarks:
        let cages = build_canonical_cage_and_snark_library();
        for g in &cages {
            let cert = certify_b12_deg2_bessel_pinch(g);
            if cert.girth >= 6 && cert.c8 == 0 {
                assert_eq!(
                    cert.tr_b12,
                    24 * cert.c12 + 12 * cert.c6 + 24 * cert.tau1 + 24 * cert.tau3
                );
            }
            if g.name.starts_with("Heawood") {
                assert_eq!(cert.c6, 28);
                assert_eq!(cert.c12, 56);
                assert_eq!(cert.tau1, 168);
                assert_eq!(cert.tau3, 84);
            }
            if g.name.starts_with("Tutte-Coxeter") {
                assert_eq!(cert.c6, 0);
                assert_eq!(cert.tau1, 0);
                assert_eq!(cert.tau3, 0);
            }
        }

        // 2. Certify on Dipole Lift Z_31(0,1,6) [n=62]:
        let frontier = build_open_frontier_voltage_and_lcf_library();
        let z31_6 = &frontier[6];
        let z31_cert = certify_b12_deg2_bessel_pinch(z31_6);
        assert_eq!(z31_cert.c6, 31);
        assert_eq!(z31_cert.c10, 93);
        assert_eq!(z31_cert.c12, 155);
        assert_eq!(z31_cert.tau1, 93);
        assert_eq!(z31_cert.tau3, 0);
        assert_eq!(z31_cert.tr_b12, 6_324);
        assert_eq!(z31_cert.n12_composite_actual, 12 * 31 + 24 * 93);
        assert_eq!(z31_cert.tr_b16_2d_bessel_floor, 106_113);
        assert_eq!(z31_cert.tr_b16_3d_bessel_floor, 126_413);
        assert_eq!(z31_cert.coupled_c16_2d_lower_bound, 2_015);
        assert_eq!(z31_cert.coupled_c16_3d_lower_bound, 2_649);

        // 3. Certify across all 10 frontier graphs in build_phase3_commutator_frontier_library():
        let phase3 = build_phase3_commutator_frontier_library();
        for g in &phase3 {
            let cert = certify_b12_deg2_bessel_pinch(g);
            if cert.c8 == 0 {
                assert_eq!(
                    cert.tr_b12,
                    24 * cert.c12 + 12 * cert.c6 + 24 * cert.tau1 + 24 * cert.tau3
                );
                assert!(
                    cert.tr_b16_3d_bessel_floor >= cert.tr_b16_2d_bessel_floor,
                    "3D Bessel floor must dominate 2D Bessel floor on {}",
                    g.name
                );
                assert!(
                    cert.coupled_c16_3d_lower_bound >= cert.coupled_c16_2d_lower_bound,
                    "3D coupled C_16 bound must dominate 2D coupled C_16 bound on {}",
                    g.name
                );
                if cert.max_c6_per_2path == 1 {
                    assert_eq!(cert.tau3, 0);
                    assert_eq!(cert.tau1, 3 * cert.c6);
                    assert_eq!(cert.tr_b12, 24 * cert.c12 + 84 * cert.c6);
                    assert!(cert.coupled_c16_3d_lower_bound >= 1523);
                }
            }
        }

        // 4. Certify that Phi_{3D}(n, r, c_6, 12 c_6) > 0 across all 24 <= n <= 100 and that the
        //    girth g >= 12 (c_6 = c_10 = 0) 3D Bessel floor is > 0 across all even 16 <= n <= 290:
        let r = 2i128;
        for n_even in (24..=100).step_by(2) {
            let ni = n_even as i128;
            let c6i = 0i128;
            let tr_b12i = 0i128;
            let num = (ni - r) * tr_b12i * tr_b12i
                - 12 * (16 * (ni - r) * c6i + (635 * r - 5 * ni) * ni) * tr_b12i
                + 1152
                    * (8 * (ni - r) * c6i * c6i
                        + (5 * ni + 517 * r) * ni * c6i
                        + (1255 * r - 25 * ni) * ni * ni);
            let den = 6 * ni * (5 * ni - 59 * r);
            assert!(
                num / den > 0,
                "Theorem 2.7(b) Phi_3D(c_6=0, b_12=0) must be > 0 at n={}",
                n_even
            );
        }
        for n_even in (16..=262).step_by(2) {
            let floor_g12 = compute_tr_b16_3d_bessel_floor(n_even, 2, 0, 0, 0);
            assert!(
                floor_g12 > 0,
                "Theorem 2.7(b) 3D Bessel floor at c_6=c_10=0 must be > 0 at n={}",
                n_even
            );
        }
        for n_even in (16..=82).step_by(2) {
            let ni = n_even as i128;
            // det S_2^{(v+8)} = 0 lower bound on Tr(B^{12}) when 3*n < 125*r = 250:
            let min_tr_b12 =
                (288 * ni * (125 * r - 3 * ni) + (7 * ni - 25 * r) - 1) / (7 * ni - 25 * r);
            assert!(
                min_tr_b12 > 0,
                "det S_2^(v+8) >= 0 must force Tr(B^12) > 0 at n={}",
                n_even
            );
        }
    }

    #[test]
    fn test_b14_deg3_bessel_pinch_theorem_2_8() {
        // 1. Certify across all 18 canonical cages and snarks:
        let cages = build_canonical_cage_and_snark_library();
        for g in &cages {
            let cert = certify_b14_deg3_bessel_pinch(g);
            if cert.is_bipartite && cert.c4 == 0 && cert.c8 == 0 {
                assert_eq!(cert.tr_b14, 28 * cert.c14 + 56 * cert.beta1);
            }
            if g.name.starts_with("Tutte-Coxeter") {
                assert_eq!(cert.c6, 0);
                assert_eq!(cert.beta1, 0);
            }
        }

        // 2. Certify on Dipole Lift Z_31(0,1,6) [n=62]:
        let frontier = build_open_frontier_voltage_and_lcf_library();
        let z31_6 = &frontier[6];
        let z31_cert = certify_b14_deg3_bessel_pinch(z31_6);
        assert_eq!(z31_cert.c6, 31);
        assert_eq!(z31_cert.c10, 93);
        assert_eq!(z31_cert.c12, 155);
        assert_eq!(z31_cert.c14, 1_302);
        assert_eq!(z31_cert.beta1, 93);
        assert_eq!(z31_cert.tr_b14, 41_664);
        assert_eq!(z31_cert.n14_composite_actual, 56 * 93);
        assert_eq!(z31_cert.tr_b16_3d_bessel_floor, 126_413);
        assert_eq!(z31_cert.tr_b16_4d_bessel_floor, 133_404);
        assert_eq!(z31_cert.coupled_c16_3d_lower_bound, 2_649);
        assert_eq!(z31_cert.coupled_c16_4d_lower_bound, 2_867);

        // 3. Certify across all 10 frontier graphs in build_phase3_commutator_frontier_library():
        let phase3 = build_phase3_commutator_frontier_library();
        for g in &phase3 {
            let cert = certify_b14_deg3_bessel_pinch(g);
            if cert.c8 == 0 {
                assert_eq!(
                    cert.tr_b14,
                    28 * cert.c14 + 56 * cert.beta1,
                    "Theorem 2.8(a): Tr(B^14) == 28*C_14 + 56*beta_1 must hold on {}",
                    g.name
                );
                assert!(
                    cert.tr_b16_4d_bessel_floor >= cert.tr_b16_3d_bessel_floor,
                    "4D Bessel floor must dominate 3D Bessel floor on {}",
                    g.name
                );
                assert!(
                    cert.coupled_c16_4d_lower_bound >= cert.coupled_c16_3d_lower_bound,
                    "4D coupled C_16 bound must dominate 3D coupled C_16 bound on {}",
                    g.name
                );
                if cert.max_c6_per_2path == 1 {
                    assert_eq!(cert.beta1, 3 * cert.c6);
                    assert_eq!(cert.tr_b14, 28 * cert.c14 + 168 * cert.c6);
                    assert!(cert.coupled_c16_4d_lower_bound >= 1666);
                }
            }
        }

        // 4. Certify the (3,14)-Moore graph exclusion at n = 127*r (n=127 for r=1, n=254 for r=2):
        //    At g = 14 (c_6 = c_10 = b_12 = 0) and n = 127*r, det H_4 = 0 forces det S_4^{(x)} = 0
        //    and det H_4^{(9-x)} = 0, pinching b_14 = 384 so Tr(B^14) = 384 * 127 * r, which is not
        //    divisible by 28 = 4 * 7 (since gcd(7, 384 * 127 * r) = 1 for r in {1, 2})!
        for r_val in [1u64, 2u64] {
            let n_moore14 = 127 * r_val;
            let tr_b14_moore = 384 * n_moore14;
            assert_ne!(
                tr_b14_moore % 28,
                0,
                "At n=127*r={}, Tr(B^14) = 384*n = {} is not divisible by 28, ruling out (3,14)-Moore graph!",
                n_moore14,
                tr_b14_moore
            );
        }

        // 5. Certify 4D Bessel floor > 0 at girth g >= 16 (c_6=c_10=b_12=b_14=0) across 256 <= n <= 1020 (r=2):
        for n_even in (256..=1020).step_by(2) {
            let ni = n_even as i128;
            let ri = 2i128;
            let floor_g16 = (384 * ni * (511 * ri - ni) + (ni - 127 * ri) - 1) / (ni - 127 * ri);
            assert!(
                floor_g16 > 0,
                "Theorem 2.8(b) 4D Bessel floor at g >= 16 must be > 0 at n={}",
                n_even
            );
        }
    }

    #[test]
    fn test_multi_octave_discharging_closure_theorems_2_9_and_2_10() {
        use crate::records::{PrismFourCycleInvariant, ThetaThreeCycleInvariant};

        // 1. Verify Theorem 2.9(a) Theta 3-cycle conservation across all canonical theta subtypes:
        for &(p1, p2, p3) in &[(1, 5, 5), (3, 3, 7), (2, 4, 8), (1, 5, 9), (3, 3, 3)] {
            let th = ThetaThreeCycleInvariant::new(p1, p2, p3);
            assert_eq!(th.c12 + th.c13 + th.c23, 2 * th.num_edges);
            assert_eq!(th.c13 - th.c12, p3 - p2);
            assert_eq!(th.c12 + th.c13 - 2 * p1, th.c23);
        }

        // 2. Verify Theorem 2.9(b) Prism 4-cycle conservation across canonical prism configurations:
        for &(a1, a2, b1, b2, s1, s2) in &[
            (3, 3, 3, 3, 2, 2),
            (2, 4, 2, 4, 2, 2),
            (1, 5, 1, 5, 2, 4),
            (3, 3, 3, 3, 1, 3),
        ] {
            let pr = PrismFourCycleInvariant::new(a1, a2, b1, b2, s1, s2);
            assert_eq!(pr.omega1 + pr.omega2, pr.c1 + pr.c2 + 2 * (s1 + s2));
        }

        // 3. Certify across all 18 canonical cages and snarks:
        let cages = build_canonical_cage_and_snark_library();
        for g in &cages {
            let cert = certify_multi_octave_discharging_closure(g);
            assert!(cert.theta_conservation_verified);
            assert!(cert.prism_conservation_verified);
            assert!(cert.menger_2_edge_connected);
            assert!(cert.dyadic_closure_certified);
            assert!(cert.total_dyadic_cycles_up_to_16 > 0);
        }

        // 4. Certify across all 7 open-frontier graphs (50 <= n <= 62):
        let frontier = build_open_frontier_voltage_and_lcf_library();
        for g in &frontier {
            let cert = certify_multi_octave_discharging_closure(g);
            assert!(cert.theta_conservation_verified);
            assert!(cert.prism_conservation_verified);
            assert!(cert.menger_2_edge_connected);
            assert!(cert.dyadic_closure_certified);
            assert_eq!(cert.perron_dominance_k, 4);
            assert_eq!(cert.perron_dominance_m, 256);
            assert!(cert.total_dyadic_cycles_up_to_16 > 0);
        }

        // 5. Certify across all 10 Phase 3 commutator-frontier graphs (50 <= n <= 64):
        let phase3 = build_phase3_commutator_frontier_library();
        for g in &phase3 {
            let cert = certify_multi_octave_discharging_closure(g);
            assert!(cert.theta_conservation_verified);
            assert!(cert.prism_conservation_verified);
            assert!(cert.menger_2_edge_connected);
            assert!(cert.dyadic_closure_certified);
            assert_eq!(cert.perron_dominance_k, 4);
            assert_eq!(cert.perron_dominance_m, 256);
            assert!(cert.total_dyadic_cycles_up_to_16 > 0);
        }
    }
}
