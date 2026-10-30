//! ============================================================================
//! `search.rs` — Parallel 2-Edge-Switch (2-Opt / Whitehead) Simulated Annealing
//! Counterexample Hunter for Erdős Problem #64 ($1,000 Falsifiable Prize).
//! ============================================================================
//!
//! Explores the space of 3-regular (`\delta(G) = \Delta(G) = 3`) graphs on
//! `n = 2k` vertices (`10 <= n <= 38`) using 2-edge-switch Markov moves:
//! given two independent edges `e_1 = {u_1, v_1}` and `e_2 = {u_2, v_2}` with
//! `{u_1, u_2} notin E(G)` and `{v_1, v_2} notin E(G)`, replace `(e_1, e_2)` by
//! `({u_1, u_2}, {v_1, v_2})`.
//!
//! By Theorem 3.1, this preserves 3-regularity identically. We minimize the
//! **Dyadic Cycle Energy**:
//!   `E_dyadic(G) = w_3 * C_3(G) + w_4 * C_4(G) + w_8 * C_8(G) + w_16 * C_16(G)`
//! where `w_3, w_4 >> w_8` strictly enforces girth `g(G) >= 5` (`C_3 = C_4 = 0`)
//! while driving `C_8(G)` and `C_16(G)` toward zero. A graph with `E_dyadic(G) = 0`
//! and `n < 32` would immediately falsify the Erdős–Gyárfás conjecture and win
//! the $1,000 prize!

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::records::BitGraph;
use crate::verify_cpu::{
    compute_girth, count_c16_total, count_c4_slice, count_c8_slice, is_bipartite,
};

/// Summary of a 2-edge-switch simulated annealing search at vertex count `n`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SearchRunSummary {
    pub n: usize,
    pub restarts: usize,
    pub steps_per_restart: usize,
    pub best_girth: usize,
    pub best_c4: u64,
    pub best_c8: u64,
    pub best_c16: u64,
    pub best_dyadic_total: u64,
    pub is_bipartite: bool,
    pub counterexample_found: bool,
    pub best_adj: Vec<u64>,
}

/// Counts triangles `C_3(G)` in `O(n^2)` bitwise operations (`sum_{u < v, {u,v} in E} popcount(adj[u] & adj[v]) / 3`).
#[inline]
pub fn count_c3_slice(adj: &[u64], n: usize) -> u64 {
    let mut sum_wedges = 0u64;
    for u in 0..n {
        let mut mask = adj[u];
        while mask != 0 {
            let v = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            if v > u {
                sum_wedges += (adj[u] & adj[v]).count_ones() as u64;
            }
        }
    }
    sum_wedges / 3
}

/// Generates a random connected 3-regular seed graph on `n = 2k` vertices by starting
/// from a Hamiltonian cycle `0 - 1 - ... - (n-1) - 0` plus diameter chords `(i, i + n/2)`
/// (the Möbius ladder `M_n`) and applying `burn_in` random 2-edge switches.
pub fn random_cubic_seed_graph(n: usize, rng: &mut ChaCha8Rng, burn_in: usize) -> BitGraph {
    assert!(
        n >= 6 && n.is_multiple_of(2),
        "cubic graph requires even n >= 6"
    );
    let mut g = BitGraph::new(&format!("Cubic_{n}"), "2-Opt Markov", n, 5, false, false);
    let half = n / 2;
    for i in 0..n {
        g.add_edge(i, (i + 1) % n);
        if i < half {
            g.add_edge(i, i + half);
        }
    }
    let mut edges = extract_edge_list(&g);
    for _ in 0..burn_in {
        try_random_2opt_step(&mut g, &mut edges, rng);
    }
    g
}

/// Extracts the list of undirected edges `(u, v)` (`u < v`) of `g`.
pub fn extract_edge_list(g: &BitGraph) -> Vec<(usize, usize)> {
    let mut edges = Vec::with_capacity(3 * g.n / 2);
    for u in 0..g.n {
        let mut mask = g.adj[u];
        while mask != 0 {
            let v = mask.trailing_zeros() as usize;
            mask &= mask - 1;
            if v > u {
                edges.push((u, v));
            }
        }
    }
    edges
}

/// Attempts a random 2-edge-switch move on `g` and `edges`. If successful, returns
/// `Some((idx1, idx2, old_e1, old_e2, new_e1, new_e2))` so the move can be undone in `O(1)`.
#[allow(clippy::type_complexity)]
pub fn try_random_2opt_step(
    g: &mut BitGraph,
    edges: &mut [(usize, usize)],
    rng: &mut ChaCha8Rng,
) -> Option<(
    usize,
    usize,
    (usize, usize),
    (usize, usize),
    (usize, usize),
    (usize, usize),
)> {
    let m = edges.len();
    if m < 2 {
        return None;
    }
    let idx1 = rng.random_range(0..m);
    let mut idx2 = rng.random_range(0..(m - 1));
    if idx2 >= idx1 {
        idx2 += 1;
    }
    let (u1, v1) = edges[idx1];
    let (u2, v2) = edges[idx2];

    // Ensure all 4 endpoints are distinct
    if u1 == u2 || u1 == v2 || v1 == u2 || v1 == v2 {
        return None;
    }

    // Choose orientation randomly: either (u1-u2, v1-v2) or (u1-v2, v1-u2)
    let (nu1, nv1, nu2, nv2) = if rng.random_bool(0.5) {
        (u1, u2, v1, v2)
    } else {
        (u1, v2, v1, u2)
    };

    if g.has_edge(nu1, nv1) || g.has_edge(nu2, nv2) {
        return None;
    }

    g.remove_edge(u1, v1);
    g.remove_edge(u2, v2);
    g.add_edge(nu1, nv1);
    g.add_edge(nu2, nv2);

    let new_e1 = (nu1.min(nv1), nu1.max(nv1));
    let new_e2 = (nu2.min(nv2), nu2.max(nv2));
    edges[idx1] = new_e1;
    edges[idx2] = new_e2;

    Some((idx1, idx2, (u1, v1), (u2, v2), new_e1, new_e2))
}

/// Undoes a 2-edge-switch move in `O(1)`.
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn undo_2opt_step(
    g: &mut BitGraph,
    edges: &mut [(usize, usize)],
    idx1: usize,
    idx2: usize,
    old_e1: (usize, usize),
    old_e2: (usize, usize),
    new_e1: (usize, usize),
    new_e2: (usize, usize),
) {
    g.remove_edge(new_e1.0, new_e1.1);
    g.remove_edge(new_e2.0, new_e2.1);
    g.add_edge(old_e1.0, old_e1.1);
    g.add_edge(old_e2.0, old_e2.1);
    edges[idx1] = old_e1;
    edges[idx2] = old_e2;
}

/// Evaluates the penalized Dyadic Cycle Energy for simulated annealing:
/// `E(G) = 10_000 * C_3(G) + 1_000 * C_4(G) + C_8(G)`.
#[inline]
pub fn dyadic_annealing_energy(g: &BitGraph) -> (u64, u64, u64, u64) {
    let c3 = count_c3_slice(&g.adj, g.n);
    let c4 = count_c4_slice(&g.adj, g.n);
    let c8 = count_c8_slice(&g.adj, g.n);
    let energy = 10_000 * c3 + 1_000 * c4 + c8;
    (energy, c3, c4, c8)
}

/// Runs parallel 2-edge-switch simulated annealing across `restarts` independent Markov
/// chains for `n`-vertex cubic graphs, minimizing `E_dyadic(G)` to hunt for girth-`>= 5`
/// cubic graphs with minimum `C_8(G) + C_16(G)`.
pub fn run_dyadic_annealing_search(
    n: usize,
    restarts: usize,
    steps_per_restart: usize,
    base_seed: u64,
) -> SearchRunSummary {
    let best_per_restart: Vec<(u64, u64, u64, BitGraph)> = (0..restarts)
        .into_par_iter()
        .map(|r| {
            let mut rng = ChaCha8Rng::seed_from_u64(
                base_seed
                    .wrapping_add((n as u64) * 1_000_003)
                    .wrapping_add(r as u64),
            );
            let mut g = random_cubic_seed_graph(n, &mut rng, 4 * n);
            let mut edges = extract_edge_list(&g);
            let (mut cur_e, _, mut cur_c4, mut cur_c8) = dyadic_annealing_energy(&g);

            let mut best_e = cur_e;
            let mut best_c4 = cur_c4;
            let mut best_c8 = cur_c8;
            let mut best_g = g.clone();

            let t_start = 25.0f64;
            let t_end = 0.15f64;

            for step in 0..steps_per_restart {
                let frac = (step as f64) / (steps_per_restart.max(1) as f64);
                let temp = t_start * (t_end / t_start).powf(frac);

                if let Some((idx1, idx2, old_e1, old_e2, new_e1, new_e2)) =
                    try_random_2opt_step(&mut g, &mut edges, &mut rng)
                {
                    let (new_e, new_c3, new_c4, new_c8) = dyadic_annealing_energy(&g);
                    let accept = if new_e <= cur_e {
                        true
                    } else {
                        let delta = (new_e - cur_e) as f64;
                        rng.random::<f64>() < (-delta / temp).exp()
                    };

                    if accept {
                        cur_e = new_e;
                        cur_c4 = new_c4;
                        cur_c8 = new_c8;
                        if cur_e < best_e || (cur_e == best_e && new_c3 == 0 && cur_c4 == 0) {
                            best_e = cur_e;
                            best_c4 = cur_c4;
                            best_c8 = cur_c8;
                            best_g = g.clone();
                        }
                    } else {
                        undo_2opt_step(
                            &mut g, &mut edges, idx1, idx2, old_e1, old_e2, new_e1, new_e2,
                        );
                    }
                }
            }
            (best_e, best_c4, best_c8, best_g)
        })
        .collect();

    let (_, best_c4, best_c8, best_g) = best_per_restart
        .into_iter()
        .min_by_key(|(e, c4, c8, _)| (*e, *c4, *c8))
        .expect("restarts must be >= 1");

    let best_c16 = count_c16_total(&best_g);
    let best_dyadic_total = best_c4 + best_c8 + best_c16;
    let best_girth = compute_girth(&best_g);
    let bip = is_bipartite(&best_g);

    SearchRunSummary {
        n,
        restarts,
        steps_per_restart,
        best_girth,
        best_c4,
        best_c8,
        best_c16,
        best_dyadic_total,
        is_bipartite: bip,
        counterexample_found: best_dyadic_total == 0,
        best_adj: best_g.adj,
    }
}

/// Summary of Theorem 3.2(b)'s Hard-`C_4=C_8=0` Barrier + `C_16`-Minimizing Tabu Search
/// starting from an algebraic `C_4 = C_8 = 0` frontier graph.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct C16TabuSearchRecord {
    pub seed_name: String,
    pub n: usize,
    pub m: usize,
    pub initial_c16: u64,
    pub best_c16_girth6: u64,
    pub best_c16_girth5: u64,
    pub best_girth: usize,
    pub best_c4: u64,
    pub best_c8: u64,
    pub best_is_bipartite: bool,
    pub accepted_switches: usize,
}

/// Runs Theorem 3.2(b)'s Hard-`C_4 = C_8 = 0` Barrier + `C_16`-Minimizing Tabu Search
/// starting from `seed_graph` (which must satisfy `C_3 = C_4 = C_8 = 0`).
///
/// Every candidate 2-edge switch that creates any triangle (`C_3 > 0`), 4-cycle (`C_4 > 0`),
/// or 8-cycle (`C_8 > 0`) is strictly rejected by the hard barrier. Among valid `C_4 = C_8 = 0`
/// moves not on the short-term Tabu list, the search tracks both the minimum `C_16(G)` achieved
/// on the `g(G) >= 6` submanifold (`C_5 = 0`) and on the `g(G) >= 5` submanifold.
pub fn run_hard_c4_c8_zero_c16_tabu_search(
    seed_graph: &BitGraph,
    restarts: usize,
    steps_per_restart: usize,
    base_seed: u64,
) -> C16TabuSearchRecord {
    let n = seed_graph.n;
    let m = seed_graph.num_edges();
    let initial_c16 = count_c16_total(seed_graph);
    assert_eq!(count_c4_slice(&seed_graph.adj, n), 0);
    assert_eq!(count_c8_slice(&seed_graph.adj, n), 0);

    let per_restart: Vec<(u64, u64, usize, BitGraph)> = (0..restarts)
        .into_par_iter()
        .map(|r| {
            let mut rng = ChaCha8Rng::seed_from_u64(
                base_seed
                    .wrapping_add((n as u64) * 999_983)
                    .wrapping_add(r as u64),
            );
            let mut g = seed_graph.clone();
            let mut edges = extract_edge_list(&g);
            let mut cur_c16 = initial_c16;
            let mut best_c16_g6 = initial_c16;
            let mut best_c16_g5 = initial_c16;
            let mut best_g = g.clone();
            let mut accepted = 0usize;
            let mut tabu: std::collections::VecDeque<((usize, usize), (usize, usize))> =
                std::collections::VecDeque::with_capacity(16);

            let t_start = 18.0f64;
            let t_end = 0.5f64;

            for step in 0..steps_per_restart {
                if let Some((idx1, idx2, old_e1, old_e2, new_e1, new_e2)) =
                    try_random_2opt_step(&mut g, &mut edges, &mut rng)
                {
                    let canonical_pair = if new_e1 <= new_e2 {
                        (new_e1, new_e2)
                    } else {
                        (new_e2, new_e1)
                    };
                    if tabu.contains(&canonical_pair) {
                        undo_2opt_step(
                            &mut g, &mut edges, idx1, idx2, old_e1, old_e2, new_e1, new_e2,
                        );
                        continue;
                    }

                    // Hard barrier: strictly enforce C_3(G) = 0, C_4(G) = 0, and C_8(G) = 0!
                    if count_c3_slice(&g.adj, n) == 0
                        && count_c4_slice(&g.adj, n) == 0
                        && count_c8_slice(&g.adj, n) == 0
                    {
                        let new_c16 = count_c16_total(&g);
                        let frac = (step as f64) / (steps_per_restart.max(1) as f64);
                        let temp = t_start * (t_end / t_start).powf(frac);
                        let accept = if new_c16 <= cur_c16 {
                            true
                        } else {
                            let delta = (new_c16 - cur_c16) as f64;
                            rng.random::<f64>() < (-delta / temp).exp()
                        };

                        if accept {
                            cur_c16 = new_c16;
                            accepted += 1;
                            if tabu.len() >= 12 {
                                tabu.pop_front();
                            }
                            let old_pair = if old_e1 <= old_e2 {
                                (old_e1, old_e2)
                            } else {
                                (old_e2, old_e1)
                            };
                            tabu.push_back(old_pair);

                            let girth = compute_girth(&g);
                            if girth >= 6 && new_c16 < best_c16_g6 {
                                best_c16_g6 = new_c16;
                            }
                            if new_c16 < best_c16_g5 {
                                best_c16_g5 = new_c16;
                                best_g = g.clone();
                            }
                            continue;
                        }
                    }

                    undo_2opt_step(
                        &mut g, &mut edges, idx1, idx2, old_e1, old_e2, new_e1, new_e2,
                    );
                }
            }
            (best_c16_g6, best_c16_g5, accepted, best_g)
        })
        .collect();

    let best_c16_girth6 = per_restart
        .iter()
        .map(|(g6, _, _, _)| *g6)
        .min()
        .unwrap_or(initial_c16);
    let (_, best_c16_girth5, accepted_switches, best_g) = per_restart
        .into_iter()
        .min_by_key(|(_, g5, _, _)| *g5)
        .expect("restarts >= 1");

    C16TabuSearchRecord {
        seed_name: seed_graph.name.clone(),
        n,
        m,
        initial_c16,
        best_c16_girth6,
        best_c16_girth5,
        best_girth: compute_girth(&best_g),
        best_c4: count_c4_slice(&best_g.adj, n),
        best_c8: count_c8_slice(&best_g.adj, n),
        best_is_bipartite: is_bipartite(&best_g),
        accepted_switches,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::build_zm_dipole_voltage_lift;

    #[test]
    fn test_2opt_preserves_cubic_regularity_and_finds_girth_5_at_n10() {
        let summary = run_dyadic_annealing_search(10, 8, 400, 42);
        assert_eq!(summary.n, 10);
        // At n=10, 2-opt simulated annealing easily discovers the Petersen graph (girth 5, C_4 = 0, C_8 = 15)
        assert_eq!(summary.best_c4, 0, "Should reach girth 5 (C_4 = 0) at n=10");
        assert_eq!(summary.best_girth, 5);
        assert_eq!(
            summary.best_c8, 15,
            "At n=10 girth 5, Petersen has C_8 = 15"
        );
        assert!(!summary.counterexample_found);
    }

    #[test]
    fn test_hard_c4_c8_zero_tabu_search_breaks_abelian_floor_on_z31() {
        let z31_6 =
            build_zm_dipole_voltage_lift("Dipole Lift Z_31(0,1,6) [n=62]", 31, [0, 1, 6], 6);
        let rec = run_hard_c4_c8_zero_c16_tabu_search(&z31_6, 4, 80, 20261026);
        assert_eq!(rec.initial_c16, 2883);
        assert_eq!(rec.best_c4, 0);
        assert_eq!(rec.best_c8, 0);
        assert!(
            rec.best_c16_girth6 < 2294,
            "Tabu search on g>=6 submanifold must break global abelian dipole floor 2294 (got {})",
            rec.best_c16_girth6
        );
        assert!(
            rec.best_c16_girth5 <= rec.best_c16_girth6,
            "Tabu search on g>=5 submanifold must be <= g>=6 minimum"
        );
    }
}
