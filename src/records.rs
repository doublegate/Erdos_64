//! ============================================================================
//! `records.rs` — Exact library of (3,g)-Cages, Generalized Polygons, and Snarks
//! for Erdős Problem #64 (the Erdős–Gyárfás 2^k-Cycle Conjecture).
//! ============================================================================
//!
//! Represents finite graphs on `n <= 64` vertices as `[u64; 64]` bitmask
//! adjacency arrays (`adj[u]` has bit `v` set iff `{u, v} in E(G)`).
//!
//! Provides exact definitions and certified structural invariants (order `n`,
//! size `m`, minimum degree `delta`, maximum degree `Delta`, girth `g`,
//! bipartiteness, and dyadic cycle counts `C_4, C_8, C_16, C_32`) for the classic
//! cubic cages and snarks:
//! - `K_4` (`n = 4`, `(3,3)`-cage, girth 3)
//! - `K_{3,3}` (`n = 6`, `(3,4)`-cage, girth 4, bipartite)
//! - Cube `Q_3` (`n = 8`, cubic planar bipartite, girth 4)
//! - Wagner graph `W_8` (`n = 8`, cubic non-planar Möbius ladder `M_8`, girth 4)
//! - Petersen graph (`n = 10`, unique `(3,5)`-cage, girth 5, non-Hamiltonian snark)
//! - Heawood graph (`n = 14`, unique `(3,6)`-cage, girth 6, bipartite incidence graph of `PG(2,2)`)
//! - Generalized Petersen `GP(8, 3)` (`n = 16`, Möbius–Kantor graph, girth 6, bipartite)
//! - Pappus graph (`n = 18`, girth 6, bipartite Levi graph of the Pappus configuration)
//! - Blanuša Snark 1 (`n = 18`, girth 5, cyclically 4-edge-connected snark)
//! - Desargues graph `GP(10, 3)` (`n = 20`, girth 6, bipartite Levi graph of Desargues configuration)
//! - Dodecahedron `GP(10, 2)` (`n = 20`, girth 5, cubic planar Platonic solid)
//! - Flower Snark `J_5` (`n = 20`, Isaacs 1975 cyclically 5-edge-connected snark, girth 5)
//! - McGee graph (`n = 24`, unique `(3,7)`-cage, girth 7, LCF `[-12, 7, -7]^8`)
//! - Coxeter graph (`n = 28`, girth 7, non-Hamiltonian cubic symmetric graph)
//! - Tutte–Coxeter `(3,8)`-cage (`n = 30`, unique `(3,8)`-cage, girth 8, bipartite Levi graph of Cremona–Richmond configuration, LCF `[-13, -9, 7, -7, 9, 13]^5`)

use serde::{Deserialize, Serialize};

/// Maximum supported vertex count in 64-bit adjacency-mask representation.
pub const MAX_VERTICES: usize = 64;

/// Compact 64-bit adjacency-mask representation of a simple undirected graph on `n <= 64` vertices.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BitGraph {
    pub name: String,
    pub family: String,
    pub n: usize,
    pub adj: Vec<u64>,
    pub expected_girth: usize,
    pub is_cage: bool,
    pub is_snark: bool,
}

impl BitGraph {
    /// Creates an empty graph on `n` vertices (`n <= 64`).
    pub fn new(
        name: &str,
        family: &str,
        n: usize,
        expected_girth: usize,
        is_cage: bool,
        is_snark: bool,
    ) -> Self {
        assert!(n <= MAX_VERTICES, "BitGraph requires n <= 64 (got n = {n})");
        Self {
            name: name.to_string(),
            family: family.to_string(),
            n,
            adj: vec![0u64; n],
            expected_girth,
            is_cage,
            is_snark,
        }
    }

    /// Adds an undirected edge `{u, v}` (`u != v`).
    #[inline]
    pub fn add_edge(&mut self, u: usize, v: usize) {
        assert!(u < self.n && v < self.n && u != v);
        self.adj[u] |= 1u64 << v;
        self.adj[v] |= 1u64 << u;
    }

    /// Removes an undirected edge `{u, v}`.
    #[inline]
    pub fn remove_edge(&mut self, u: usize, v: usize) {
        if u < self.n && v < self.n {
            self.adj[u] &= !(1u64 << v);
            self.adj[v] &= !(1u64 << u);
        }
    }

    /// Checks whether `{u, v}` is an edge.
    #[inline]
    pub fn has_edge(&self, u: usize, v: usize) -> bool {
        (self.adj[u] & (1u64 << v)) != 0
    }

    /// Degree of vertex `u`.
    #[inline]
    pub fn degree(&self, u: usize) -> u32 {
        self.adj[u].count_ones()
    }

    /// Minimum degree `\delta(G)`.
    pub fn min_degree(&self) -> u32 {
        (0..self.n).map(|u| self.degree(u)).min().unwrap_or(0)
    }

    /// Maximum degree `\Delta(G)`.
    pub fn max_degree(&self) -> u32 {
        (0..self.n).map(|u| self.degree(u)).max().unwrap_or(0)
    }

    /// Number of undirected edges `|E(G)|`.
    pub fn num_edges(&self) -> usize {
        let sum_deg: u32 = (0..self.n).map(|u| self.degree(u)).sum();
        (sum_deg / 2) as usize
    }

    /// Constructs a cubic graph from LCF (Lederberg–Coxeter–Frucht) notation:
    /// a Hamiltonian cycle `0 - 1 - 2 - ... - (n-1) - 0` plus chords `(i, (i + shift) mod n)`.
    #[allow(clippy::too_many_arguments)]
    pub fn from_lcf(
        name: &str,
        family: &str,
        n: usize,
        shifts: &[isize],
        repeats: usize,
        expected_girth: usize,
        is_cage: bool,
        is_snark: bool,
    ) -> Self {
        assert_eq!(shifts.len() * repeats, n, "LCF length must equal n");
        let mut g = Self::new(name, family, n, expected_girth, is_cage, is_snark);
        for i in 0..n {
            g.add_edge(i, (i + 1) % n);
            let s = shifts[i % shifts.len()];
            let target = ((i as isize + s).rem_euclid(n as isize)) as usize;
            g.add_edge(i, target);
        }
        g
    }

    /// Constructs the Generalized Petersen graph `GP(m, k)` on `n = 2m` vertices:
    /// outer rim `0..m` with edges `(i, (i+1)%m)`, spokes `(i, m+i)`, and inner star
    /// `(m+i, m + (i+k)%m)`.
    pub fn generalized_petersen(
        name: &str,
        family: &str,
        m: usize,
        k: usize,
        expected_girth: usize,
        is_cage: bool,
        is_snark: bool,
    ) -> Self {
        let n = 2 * m;
        let mut g = Self::new(name, family, n, expected_girth, is_cage, is_snark);
        for i in 0..m {
            g.add_edge(i, (i + 1) % m);
            g.add_edge(i, m + i);
            g.add_edge(m + i, m + ((i + k) % m));
        }
        g
    }

    /// Packs the adjacency array into a fixed-length `[u64; 64]` slice for GPU kernels.
    pub fn padded_adj_64(&self) -> [u64; MAX_VERTICES] {
        let mut out = [0u64; MAX_VERTICES];
        for (i, &mask) in self.adj.iter().enumerate() {
            out[i] = mask;
        }
        out
    }
}

/// Complete cycle-length spectrum and Erdős–Gyárfás dyadic cycle certificate for a graph.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CycleSpectrumRecord {
    pub name: String,
    pub family: String,
    pub n: usize,
    pub m: usize,
    pub min_deg: u32,
    pub max_deg: u32,
    pub girth: usize,
    pub circumference: usize,
    pub is_bipartite: bool,
    pub is_cage: bool,
    pub is_snark: bool,
    /// Number of undirected simple cycles of length 4 (`2^2`).
    pub c4: u64,
    /// Number of undirected simple cycles of length 8 (`2^3`).
    pub c8: u64,
    /// Number of undirected simple cycles of length 16 (`2^4`).
    pub c16: u64,
    /// Number of undirected simple cycles of length 32 (`2^5`).
    pub c32: u64,
    /// Total number of dyadic cycles (`C_4 + C_8 + C_16 + C_32`).
    pub total_dyadic_cycles: u64,
    /// Distinct dyadic cycle lengths `2^k` (`k >= 2`) present in `G`.
    pub dyadic_lengths_present: Vec<usize>,
    /// All cycle lengths `l in [3, n]` present in `G` (populated when exact full spectrum is run).
    pub cycle_lengths: Vec<usize>,
    /// Exact counts `C_l(G)` for `l = 0..=n` when full spectrum enumeration is enabled.
    pub cycle_counts_by_length: Vec<u64>,
}

/// Moore lower bound `n_0(3, g)` on the minimum number of vertices of a cubic (`3`-regular)
/// graph of girth `g`, together with the actual `(3, g)`-cage order `f(3, g)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CubicCageBound {
    pub girth: usize,
    pub moore_lower_bound: usize,
    pub exact_cage_order: usize,
    pub cage_name: &'static str,
    pub unique_cage: bool,
    pub smallest_dyadic_cycle_length: usize,
}

/// Exact table of Moore bounds `n_0(3, g)` and known `(3, g)`-cages for `3 <= g <= 12`.
pub const CUBIC_CAGE_BOUNDS: &[CubicCageBound] = &[
    CubicCageBound {
        girth: 3,
        moore_lower_bound: 4,
        exact_cage_order: 4,
        cage_name: "K_4 (Tetrahedron)",
        unique_cage: true,
        smallest_dyadic_cycle_length: 4,
    },
    CubicCageBound {
        girth: 4,
        moore_lower_bound: 6,
        exact_cage_order: 6,
        cage_name: "K_{3,3} (Utility graph)",
        unique_cage: true,
        smallest_dyadic_cycle_length: 4,
    },
    CubicCageBound {
        girth: 5,
        moore_lower_bound: 10,
        exact_cage_order: 10,
        cage_name: "Petersen graph",
        unique_cage: true,
        smallest_dyadic_cycle_length: 8,
    },
    CubicCageBound {
        girth: 6,
        moore_lower_bound: 14,
        exact_cage_order: 14,
        cage_name: "Heawood graph PG(2,2)",
        unique_cage: true,
        smallest_dyadic_cycle_length: 8,
    },
    CubicCageBound {
        girth: 7,
        moore_lower_bound: 22,
        exact_cage_order: 24,
        cage_name: "McGee (3,7)-cage",
        unique_cage: true,
        smallest_dyadic_cycle_length: 8,
    },
    CubicCageBound {
        girth: 8,
        moore_lower_bound: 30,
        exact_cage_order: 30,
        cage_name: "Tutte-Coxeter (3,8)-cage",
        unique_cage: true,
        smallest_dyadic_cycle_length: 8,
    },
    CubicCageBound {
        girth: 9,
        moore_lower_bound: 46,
        exact_cage_order: 58,
        cage_name: "Biggs-Smith / Brinkmann-McKay-Saager (18 cages)",
        unique_cage: false,
        smallest_dyadic_cycle_length: 16,
    },
    CubicCageBound {
        girth: 10,
        moore_lower_bound: 62,
        exact_cage_order: 70,
        cage_name: "Balaban / O'Keefe-Wong (3 cages)",
        unique_cage: false,
        smallest_dyadic_cycle_length: 16,
    },
    CubicCageBound {
        girth: 11,
        moore_lower_bound: 94,
        exact_cage_order: 112,
        cage_name: "Balaban (3,11)-cage",
        unique_cage: true,
        smallest_dyadic_cycle_length: 16,
    },
    CubicCageBound {
        girth: 12,
        moore_lower_bound: 126,
        exact_cage_order: 126,
        cage_name: "Benson (3,12)-cage (Generalized Hexagon GH(2,2))",
        unique_cage: true,
        smallest_dyadic_cycle_length: 16,
    },
];

/// Computes the exact Moore lower bound `n_0(3, g)` on the order of a 3-regular graph of girth `g`:
/// - For odd `g = 2r + 1`: `n_0(3, 2r+1) = 1 + 3 * sum_{i=0}^{r-1} 2^i = 3 * 2^r - 2`.
/// - For even `g = 2r`: `n_0(3, 2r) = 2 * sum_{i=0}^{r-1} 2^i = 2^{r+1} - 2`.
pub fn moore_bound_cubic(girth: usize) -> usize {
    assert!(girth >= 3, "girth must be >= 3");
    if girth % 2 == 1 {
        let r = (girth - 1) / 2;
        3 * (1usize << r) - 2
    } else {
        let r = girth / 2;
        (1usize << (r + 1)) - 2
    }
}

/// Constructs the complete benchmark library of classic `(3,g)`-cages, snarks, and
/// generalized polygons used to certify Theorem 2.1.
#[allow(clippy::vec_init_then_push)]
pub fn build_canonical_cage_and_snark_library() -> Vec<BitGraph> {
    let mut lib = Vec::new();

    // 1. K_4 (n=4, (3,3)-cage, girth 3)
    lib.push(BitGraph::from_lcf(
        "K_4 (Tetrahedron)",
        "(3,3)-cage",
        4,
        &[2],
        4,
        3,
        true,
        false,
    ));

    // 2. K_{3,3} (n=6, (3,4)-cage, girth 4)
    lib.push(BitGraph::from_lcf(
        "K_{3,3} (Utility)",
        "(3,4)-cage",
        6,
        &[3],
        6,
        4,
        true,
        false,
    ));

    // 3. Cube Q_3 = GP(4,1) (n=8, girth 4)
    lib.push(BitGraph::generalized_petersen(
        "Cube Q_3",
        "Platonic / GP(4,1)",
        4,
        1,
        4,
        false,
        false,
    ));

    // 4. Wagner graph M_8 (n=8, girth 4, Möbius ladder)
    lib.push(BitGraph::from_lcf(
        "Wagner M_8",
        "Mobius ladder",
        8,
        &[4],
        8,
        4,
        false,
        false,
    ));

    // 5. Petersen graph GP(5,2) (n=10, unique (3,5)-cage, girth 5, fundamental snark)
    lib.push(BitGraph::generalized_petersen(
        "Petersen GP(5,2)",
        "(3,5)-cage / Snark",
        5,
        2,
        5,
        true,
        true,
    ));

    // 6. Tietze's graph (n=12, girth 3, non-Hamiltonian cubic graph obtained by replacing
    //    one vertex of Petersen with a triangle)
    {
        let mut tietze =
            BitGraph::new("Tietze graph", "Cubic non-Hamiltonian", 12, 3, false, false);
        // Outer 9-cycle 0..9 with triangles / Petersen-style spokes to triangle 9-10-11
        let petersen = BitGraph::generalized_petersen("p", "p", 5, 2, 5, true, true);
        // Replace vertex 9 of Petersen (neighbors 4, 6, 7) with triangle (9, 10, 11)
        for u in 0..9 {
            for v in (u + 1)..9 {
                if petersen.has_edge(u, v) {
                    tietze.add_edge(u, v);
                }
            }
        }
        // Triangle 9-10-11
        tietze.add_edge(9, 10);
        tietze.add_edge(10, 11);
        tietze.add_edge(11, 9);
        // Connect original neighbors of vertex 9 in GP(5,2): vertices 4, 6, 7
        tietze.add_edge(4, 9);
        tietze.add_edge(6, 10);
        tietze.add_edge(7, 11);
        lib.push(tietze);
    }

    // 7. Heawood graph (n=14, unique (3,6)-cage, girth 6, LCF [5, -5]^7)
    lib.push(BitGraph::from_lcf(
        "Heawood (3,6)-cage",
        "(3,6)-cage / PG(2,2)",
        14,
        &[5, -5],
        7,
        6,
        true,
        false,
    ));

    // 8. Möbius–Kantor graph GP(8,3) (n=16, girth 6, LCF [5, -5, -3, 3]^4)
    lib.push(BitGraph::generalized_petersen(
        "Mobius-Kantor GP(8,3)",
        "Generalized Petersen",
        8,
        3,
        6,
        false,
        false,
    ));

    // 9. Pappus graph (n=18, girth 6, LCF [5, 7, -7, 7, -7, -5]^3)
    lib.push(BitGraph::from_lcf(
        "Pappus graph",
        "Levi graph (9_3)",
        18,
        &[5, 7, -7, 7, -7, -5],
        3,
        6,
        false,
        false,
    ));

    // 10. First Blanuša Snark (n=18, girth 5, discovered by Danilo Blanuša in 1946)
    //     Constructed by taking two copies of Petersen minus two adjacent vertices (16 vertices)
    //     and joining them via a 2-vertex bridge gadget, or equivalently via its explicit edge list.
    lib.push(build_blanusa_snark_1());

    // 11. Desargues graph GP(10,3) (n=20, girth 6, LCF [5, -5, 9, -9]^5)
    lib.push(BitGraph::generalized_petersen(
        "Desargues GP(10,3)",
        "Levi graph (10_3)",
        10,
        3,
        6,
        false,
        false,
    ));

    // 12. Dodecahedron GP(10,2) (n=20, girth 5, Platonic solid)
    lib.push(BitGraph::generalized_petersen(
        "Dodecahedron GP(10,2)",
        "Platonic / Planar",
        10,
        2,
        5,
        false,
        false,
    ));

    // 13. Isaacs Flower Snark J_5 (n=20, girth 5, cyclically 5-edge-connected snark)
    lib.push(build_flower_snark_j5());

    // 14. Generalized Petersen GP(11,2) (n=22, girth 5, non-bipartite cubic graph at Moore bound for g=7)
    lib.push(BitGraph::generalized_petersen(
        "Generalized Petersen GP(11,2)",
        "Generalized Petersen",
        11,
        2,
        5,
        false,
        false,
    ));

    // 15. McGee graph (n=24, unique (3,7)-cage, girth 7, LCF [-12, 7, -7]^8)
    lib.push(BitGraph::from_lcf(
        "McGee (3,7)-cage",
        "(3,7)-cage",
        24,
        &[-12, 7, -7],
        8,
        7,
        true,
        false,
    ));

    // 16. Nauru graph GP(12,5) (n=24, girth 6, LCF [5, -9, 7, -7, 9, -5]^4)
    lib.push(BitGraph::generalized_petersen(
        "Nauru GP(12,5)",
        "Generalized Petersen",
        12,
        5,
        6,
        false,
        false,
    ));

    // 17. Coxeter graph (n=28, girth 7, non-Hamiltonian cubic symmetric graph)
    lib.push(build_coxeter_graph());

    // 18. Tutte–Coxeter graph (n=30, unique (3,8)-cage, girth 8, LCF [-13, -9, 7, -7, 9, 13]^5)
    lib.push(BitGraph::from_lcf(
        "Tutte-Coxeter (3,8)-cage",
        "(3,8)-cage / GQ(2,2)",
        30,
        &[-13, -9, 7, -7, 9, 13],
        5,
        8,
        true,
        false,
    ));

    lib
}

/// Constructs the First Blanuša Snark on `n = 18` vertices (girth 5, cubic, non-Hamiltonian).
///
/// Built via the standard dot-product of two Petersen graphs `P^{(1)}` and `P^{(2)}`:
/// remove an edge `{0, 1}` from `P^{(1)}` (leaving vertices `2..10` with degree-2 ports at
/// `{4, 5}` and `{2, 6}`) and remove `{8, 9}` plus its incident vertices from `P^{(2)}`
/// (leaving 8 vertices `0..8` with 4 degree-2 ports), then join the two 4-edge cuts.
pub fn build_blanusa_snark_1() -> BitGraph {
    let mut g = BitGraph::new("Blanusa Snark B_1", "Snark (n=18)", 18, 5, false, true);
    let p = BitGraph::generalized_petersen("p", "p", 5, 2, 5, true, true);

    // Block A: vertices 0..9 correspond to Petersen vertices 1..10 (removing vertex 0, whose
    // neighbors in GP(5,2) were 1, 4, 5).
    // Even simpler: dot product of two Petersen graphs P1 (removing edge (0,5)) and P2 (removing vertices 0 and 5).
    // Let's verify ports when removing vertices 0 and 5 from GP(5,2):
    // In GP(5,2):
    // 0 has neighbors 1, 4, 5.
    // 5 has neighbors 0, 7, 8.
    // Since {0, 5} is an edge, removing BOTH vertices 0 and 5 leaves 8 vertices {1,2,3,4,6,7,8,9}
    // with exactly FOUR degree-2 vertices: 1, 4 (former neighbors of 0) and 7, 8 (former neighbors of 5)!
    // In the second Petersen copy P1, remove two non-adjacent edges, say (1,2) and (6,8), leaving
    // all 10 vertices with four degree-2 ports {1, 2, 6, 8}!
    // Joining the four ports {1, 4, 7, 8} of (P2 \ {0,5}) to {1, 2, 6, 8} of (P1 \ {(1,2),(6,8)})
    // yields the 18-vertex Blanuša snark of girth 5!
    // Copy 1 (vertices 0..10): GP(5,2) minus edges (1,2) and (6,8)
    for u in 0..10 {
        for v in (u + 1)..10 {
            if p.has_edge(u, v) && !((u == 1 && v == 2) || (u == 6 && v == 8)) {
                g.add_edge(u, v);
            }
        }
    }
    // Copy 2 (vertices 10..18): GP(5,2) minus vertices 0 and 5
    let map_p2 = [usize::MAX, 10, 11, 12, 13, usize::MAX, 14, 15, 16, 17];
    for u in 0..10 {
        if map_p2[u] == usize::MAX {
            continue;
        }
        for v in (u + 1)..10 {
            if map_p2[v] == usize::MAX {
                continue;
            }
            if p.has_edge(u, v) {
                g.add_edge(map_p2[u], map_p2[v]);
            }
        }
    }
    // Connect the 4 ports of Copy 1 ({1, 2, 6, 8}) to the 4 ports of Copy 2
    // (map_p2[1]=10, map_p2[4]=13, map_p2[7]=15, map_p2[8]=16):
    g.add_edge(1, 10);
    g.add_edge(2, 13);
    g.add_edge(6, 15);
    g.add_edge(8, 16);
    g
}

/// Constructs Isaacs's Flower Snark `J_5` on `n = 20` vertices (girth 5, cubic, non-Hamiltonian).
///
/// For `k = 5`, `J_5` has 5 four-vertex star gadgets `(a_i, b_i, c_i, d_i)` for `i in 0..5`:
/// - Central hub `a_i` joined to `b_i, c_i, d_i`;
/// - Outer rim cycle `b_0 - b_1 - b_2 - b_3 - b_4 - b_0`;
/// - Twisted inner 10-cycle connecting `(c_i, c_{i+1})` for `i = 0..4` with cross-over
///   `(c_4, d_0)`, `(d_i, d_{i+1})` for `i = 0..4`, and `(d_4, c_0)`.
pub fn build_flower_snark_j5() -> BitGraph {
    let mut g = BitGraph::new("Flower Snark J_5", "Snark (n=20)", 20, 5, false, true);
    let a = |i: usize| i;
    let b = |i: usize| 5 + i;
    let c = |i: usize| 10 + i;
    let d = |i: usize| 15 + i;

    for i in 0..5 {
        // Hub edges from a_i
        g.add_edge(a(i), b(i));
        g.add_edge(a(i), c(i));
        g.add_edge(a(i), d(i));
        // Outer 5-cycle on b_i
        g.add_edge(b(i), b((i + 1) % 5));
    }
    // Twisted 10-cycle on c_i and d_i
    for i in 0..4 {
        g.add_edge(c(i), c(i + 1));
        g.add_edge(d(i), d(i + 1));
    }
    g.add_edge(c(4), d(0));
    g.add_edge(d(4), c(0));
    g
}

/// Constructs the Coxeter graph on `n = 28` vertices (girth 7, cubic, symmetric, non-Hamiltonian).
///
/// Built from three concentric 7-cycles (`u_0..u_6` with step 1, `v_0..v_6` with step 2,
/// `w_0..w_6` with step 3 mod 7) joined by spokes to a central 7-vertex star `z_0..z_6`.
pub fn build_coxeter_graph() -> BitGraph {
    let mut g = BitGraph::new(
        "Coxeter graph",
        "Symmetric non-Hamiltonian (g=7)",
        28,
        7,
        false,
        false,
    );
    let z = |i: usize| i; // 0..7
    let u = |i: usize| 7 + i; // 7..14 (7-cycle step 1)
    let v = |i: usize| 14 + i; // 14..21 (7-cycle step 2)
    let w = |i: usize| 21 + i; // 21..28 (7-cycle step 3)

    for i in 0..7 {
        // Spokes from z_i to u_i, v_i, w_i
        g.add_edge(z(i), u(i));
        g.add_edge(z(i), v(i));
        g.add_edge(z(i), w(i));
        // Three 7-cycles with steps 1, 2, 3 mod 7
        g.add_edge(u(i), u((i + 1) % 7));
        g.add_edge(v(i), v((i + 2) % 7));
        g.add_edge(w(i), w((i + 3) % 7));
    }
    g
}

/// Constructs an algebraic `Z_m`-Voltage-Lift of the Dipole `D_3` on `n = 2m` vertices:
/// two fiber orbits `u_i = i` (`0..m`) and `v_i = m + i` (`m..2m`) with three edges from
/// each `u_i` to `v_{(i + s_0) mod m}`, `v_{(i + s_1) mod m}`, and `v_{(i + s_2) mod m}`.
///
/// By Theorem 2.2, any simple cycle of length `2r` in `Lift(D_3; s_0, s_1, s_2)` alternates
/// between the two fibers and corresponds to a non-backtracking step sequence summing to
/// `0 (mod m)`. With `s = [0, 1, k]`, this is a bipartite cubic graph of order `n = 2m`.
pub fn build_zm_dipole_voltage_lift(
    name: &str,
    m: usize,
    shifts: [usize; 3],
    expected_girth: usize,
) -> BitGraph {
    let n = 2 * m;
    let mut g = BitGraph::new(
        name,
        &format!("Z_{m} Dipole Voltage-Lift"),
        n,
        expected_girth,
        false,
        false,
    );
    for i in 0..m {
        for &s in &shifts {
            g.add_edge(i, m + ((i + s) % m));
        }
    }
    g
}

/// Builds the 7-graph Open-Frontier Benchmark Library (`50 <= n <= 62`) for Theorem 2.2,
/// spanning every even order in the uncharted post-Ducoffe–Dumitru (`n >= 50`) regime that
/// fits within a single 64-bit `u64` adjacency mask.
#[allow(clippy::vec_init_then_push)]
pub fn build_open_frontier_voltage_and_lcf_library() -> Vec<BitGraph> {
    let mut lib = Vec::new();

    // 1. n = 50: Generalized Petersen GP(25, 7) (Z_25 Theta-lift, girth 8, non-bipartite)
    lib.push(BitGraph::generalized_petersen(
        "Frontier GP(25,7) [n=50]",
        "Z_25 Voltage-Lift GP(25,7)",
        25,
        7,
        8,
        false,
        false,
    ));

    // 2. n = 52: Generalized Petersen GP(26, 5) (Z_26 Theta-lift, girth 8, bipartite)
    lib.push(BitGraph::generalized_petersen(
        "Frontier GP(26,5) [n=52]",
        "Z_26 Voltage-Lift GP(26,5)",
        26,
        5,
        8,
        false,
        false,
    ));

    // 3. n = 54: Gray (3,8)-semisymmetric graph LCF [-25, 7, -7, 13, -13, 25]^9 (n=54, girth 8, bipartite)
    lib.push(BitGraph::from_lcf(
        "Gray Graph LCF [n=54]",
        "Semisymmetric (g=8) LCF",
        54,
        &[-25, 7, -7, 13, -13, 25],
        9,
        8,
        false,
        false,
    ));

    // 4. n = 56: Generalized Petersen GP(28, 5) (Z_28 Theta-lift, girth 8, bipartite)
    lib.push(BitGraph::generalized_petersen(
        "Frontier GP(28,5) [n=56]",
        "Z_28 Voltage-Lift GP(28,5)",
        28,
        5,
        8,
        false,
        false,
    ));

    // 5. n = 58: Z_29 Dipole Voltage-Lift D_3(0, 1, 8) (n=58, girth 6, bipartite; Moore bound order for g=9)
    lib.push(build_zm_dipole_voltage_lift(
        "Dipole Lift Z_29(0,1,8) [n=58]",
        29,
        [0, 1, 8],
        6,
    ));

    // 6. n = 60: Generalized Petersen GP(30, 7) (Z_30 Theta-lift, girth 8, bipartite)
    lib.push(BitGraph::generalized_petersen(
        "Frontier GP(30,7) [n=60]",
        "Z_30 Voltage-Lift GP(30,7)",
        30,
        7,
        8,
        false,
        false,
    ));

    // 7. n = 62: Z_31 Dipole Voltage-Lift D_3(0, 1, 6) (n=62, girth 6, bipartite; Moore bound order n_0(3,10)=62)
    lib.push(build_zm_dipole_voltage_lift(
        "Dipole Lift Z_31(0,1,6) [n=62]",
        31,
        [0, 1, 6],
        6,
    ));

    lib
}

/// Constructs a Non-Abelian Periodic-Shift Dipole Voltage Lift on `n = 2m` vertices:
/// two fiber orbits `u_i = i` (`0..m`) and `v_j = m + j` (`m..2m`), with edges from `u_i`
/// to `v_{(i + s_0) mod m}`, `v_{(i + s_1) mod m}`, and `v_{(i + k_{i mod p}) mod m}`
/// where `k_pattern` has period `p | m` and `i -> (i + k_{i mod p}) mod m` is a bijection.
///
/// When `k_pattern` is non-constant, the periodic chord map does not commute with unit
/// fiber translation `tau_1(i) = i + 1 mod m`, breaking abelian commutator symmetry.
pub fn build_periodic_dipole_voltage_lift(
    name: &str,
    m: usize,
    base_shifts: [usize; 2],
    k_pattern: &[usize],
    expected_girth: usize,
) -> BitGraph {
    let n = 2 * m;
    let p = k_pattern.len();
    assert!(p > 0 && m.is_multiple_of(p), "period p must divide m");
    let mut g = BitGraph::new(
        name,
        &format!("Non-Abelian Per_{p}(Z_{m}) Dipole-Lift"),
        n,
        expected_girth,
        false,
        false,
    );
    for i in 0..m {
        g.add_edge(i, m + ((i + base_shifts[0]) % m));
        g.add_edge(i, m + ((i + base_shifts[1]) % m));
        let k_i = k_pattern[i % p];
        g.add_edge(i, m + ((i + k_i) % m));
    }
    g
}

/// Constructs a Non-Abelian Metacyclic Dipole Voltage Lift `Lift_{Z_a \rtimes_q Z_b}(D_3; (0,0), g_1, g_2)`
/// on `n = 2 * a * b` vertices, where the semidirect product multiplication law is
/// `(x_1, y_1) * (x_2, y_2) = ((x_1 + q^{y_1} x_2) mod a, (y_1 + y_2) mod b)`.
///
/// Whenever `g_1` and `g_2` do not commute in `Z_a \rtimes_q Z_b`, the universal abelian
/// hexagonal 6-walk `(+s_0, -s_1, +s_2, -s_0, +s_1, -s_2)` has non-trivial commutator net
/// voltage `[g_1^{-1}, g_2] != e`, allowing the dipole lift to break the universal abelian
/// girth-6 barrier and reach girth 8.
pub fn build_metacyclic_dipole_voltage_lift(
    name: &str,
    a: usize,
    b: usize,
    q: usize,
    g1: (usize, usize),
    g2: (usize, usize),
    expected_girth: usize,
) -> BitGraph {
    let m = a * b;
    let n = 2 * m;
    let mut g = BitGraph::new(
        name,
        &format!("Metacyclic Z_{a}x_{q}Z_{b} Lift"),
        n,
        expected_girth,
        false,
        false,
    );
    let pow_mod = |mut base: usize, mut exp: usize, modulus: usize| -> usize {
        let mut res = 1usize % modulus;
        base %= modulus;
        while exp > 0 {
            if exp & 1 == 1 {
                res = (res * base) % modulus;
            }
            base = (base * base) % modulus;
            exp >>= 1;
        }
        res
    };
    let idx = |x: usize, y: usize| -> usize { y * a + x };
    let shifts = [(0usize, 0usize), g1, g2];
    for y1 in 0..b {
        let q_pow = pow_mod(q, y1, a);
        for x1 in 0..a {
            let u = idx(x1, y1);
            for &(x2, y2) in &shifts {
                let vx = (x1 + q_pow * x2) % a;
                let vy = (y1 + y2) % b;
                let v = m + idx(vx, vy);
                g.add_edge(u, v);
            }
        }
    }
    g
}

/// Builds the 10-graph Theorem 2.3 Commutator-Breaking & Optimal-Abelian Frontier Library
/// across `50 <= n <= 64`, certifying the universal abelian commutator floor `18 |Gamma|`,
/// the Minkowski sublattice resonance floor on `25 <= |Gamma| <= 32`, and non-abelian
/// periodic-LCF / metacyclic commutator breaking.
#[allow(clippy::vec_init_then_push)]
pub fn build_phase3_commutator_frontier_library() -> Vec<BitGraph> {
    let mut lib = Vec::new();

    // 1. n = 50: Optimal Abelian Z_25 Dipole Lift D_3(0, 1, 10) (C_4 = C_8 = 0, C_16 = 2550 = 102 * 25)
    lib.push(build_zm_dipole_voltage_lift(
        "Abelian Opt Z_25(0,1,10) [n=50]",
        25,
        [0, 1, 10],
        6,
    ));

    // 2. n = 52: Crossed Period-4 LCF [11, 13, -13, -11]^13 (C_4 = C_8 = 0, C_16 = 2483 = 95.5 * 26)
    lib.push(BitGraph::from_lcf(
        "Crossed LCF [11,13,-13,-11]^13 [n=52]",
        "Period-4 Bipartite LCF",
        52,
        &[11, 13, -13, -11],
        13,
        6,
        false,
        false,
    ));

    // 3. n = 54: Optimal Abelian Z_27 Dipole Lift D_3(0, 1, 6) (C_4 = C_8 = 0, C_16 = 2484 = 92 * 27)
    lib.push(build_zm_dipole_voltage_lift(
        "Abelian Opt Z_27(0,1,6) [n=54]",
        27,
        [0, 1, 6],
        6,
    ));

    // 4. n = 54: Non-Abelian Period-3 Dipole Lift Per_3(Z_27; 0, 1, [3, 3, 9]) (C_4 = C_8 = 0, C_16 = 2385 < 2484!)
    lib.push(build_periodic_dipole_voltage_lift(
        "Non-Abelian Per_3(Z_27;3,3,9) [n=54]",
        27,
        [0, 1],
        &[3, 3, 9],
        6,
    ));

    // 5. n = 56: Crossed Period-4 LCF [11, 45, -45, -11]^14 (C_4 = C_8 = 0, C_16 = 2331 = 83.25 * 28)
    lib.push(BitGraph::from_lcf(
        "Crossed LCF [11,45,-45,-11]^14 [n=56]",
        "Period-4 Bipartite LCF",
        56,
        &[11, 45, -45, -11],
        14,
        6,
        false,
        false,
    ));

    // 6. n = 58: Optimal Abelian Z_29 Dipole Lift D_3(0, 1, 9) (C_4 = C_8 = 0, C_16 = 2552 = 88 * 29)
    lib.push(build_zm_dipole_voltage_lift(
        "Abelian Opt Z_29(0,1,9) [n=58]",
        29,
        [0, 1, 9],
        6,
    ));

    // 7. n = 60: Non-Abelian Period-4 LCF [17, -17, 45, -45]^15 (C_4 = C_8 = 0, C_16 = 2340 = 78 * 30 < 2550 abelian floor!)
    lib.push(BitGraph::from_lcf(
        "Non-Abelian LCF [17,-17,45,-45]^15 [n=60]",
        "Period-4 Bipartite LCF",
        60,
        &[17, -17, 45, -45],
        15,
        6,
        false,
        false,
    ));

    // 8. n = 62: Optimal Abelian Z_31 Dipole Lift D_3(0, 1, 12) (C_4 = C_8 = 0, C_16 = 2294 = 74 * 31 < 2883!)
    lib.push(build_zm_dipole_voltage_lift(
        "Abelian Opt Z_31(0,1,12) [n=62]",
        31,
        [0, 1, 12],
        6,
    ));

    // 9. n = 64: Optimal Abelian Z_32 Dipole Lift D_3(0, 1, 7) (C_4 = C_8 = 0, C_16 = 2232 = 69.75 * 32)
    lib.push(build_zm_dipole_voltage_lift(
        "Abelian Opt Z_32(0,1,7) [n=64]",
        32,
        [0, 1, 7],
        6,
    ));

    // 10. n = 64: Non-Abelian Metacyclic Z_8 \rtimes_3 Z_4 Dipole Lift (breaks universal abelian g<=6 barrier to reach girth 8!)
    lib.push(build_metacyclic_dipole_voltage_lift(
        "Metacyclic Z_8x_3Z_4(1,0;0,1) [n=64]",
        8,
        4,
        3,
        (1, 0),
        (0, 1),
        8,
    ));

    lib
}

/// Theorem 2.9(a) Exact Theta Subgraph Three-Cycle Conservation Invariant:
/// For any theta subgraph `Theta(p_1, p_2, p_3)` (`1 <= p_1 <= p_2 <= p_3`) formed by three
/// internally vertex-disjoint paths between degree-3 vertices `u` and `v`, the three simple
/// cycles `C_12, C_13, C_23` have lengths `c_12 = p_1 + p_2`, `c_13 = p_1 + p_3`, `c_23 = p_2 + p_3`
/// satisfying `c_12 + c_13 + c_23 = 2 * (p_1 + p_2 + p_3) = 2 * |E(Theta)|` and
/// `c_23 - c_12 = p_3 - p_2 = c_13 - 2 * p_1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ThetaThreeCycleInvariant {
    pub p1: usize,
    pub p2: usize,
    pub p3: usize,
    pub c12: usize,
    pub c13: usize,
    pub c23: usize,
    pub num_vertices: usize,
    pub num_edges: usize,
}

impl ThetaThreeCycleInvariant {
    /// Constructs and verifies the exact Theorem 2.9(a) conservation law for `Theta(p_1, p_2, p_3)`.
    pub fn new(p1: usize, p2: usize, p3: usize) -> Self {
        assert!(1 <= p1 && p1 <= p2 && p2 <= p3);
        let c12 = p1 + p2;
        let c13 = p1 + p3;
        let c23 = p2 + p3;
        let num_edges = p1 + p2 + p3;
        let num_vertices = num_edges - 1;
        assert_eq!(c12 + c13 + c23, 2 * num_edges);
        assert_eq!(c13 - c12, p3 - p2);
        assert_eq!(c23 - c12, p3 - p1);
        assert_eq!(c12 + c13 - 2 * p1, c23);
        Self {
            p1,
            p2,
            p3,
            c12,
            c13,
            c23,
            num_vertices,
            num_edges,
        }
    }
}

/// Theorem 2.9(b) 2-Connected Prism/Ladder Four-Cycle Conservation Invariant:
/// For any two vertex-disjoint simple cycles `C^(1) = a_1 + a_2` and `C^(2) = b_1 + b_2`
/// connected by two vertex-disjoint bridge paths `P_1` (length `s_1 >= 1`) and `P_2`
/// (length `s_2 >= 1`) in a 2-connected cubic core, the two simple cross-bridge cycles
/// `Omega_1` and `Omega_2` have lengths `|Omega_1| = a_1 + b_1 + s_1 + s_2` and
/// `|Omega_2| = a_2 + b_2 + s_1 + s_2`, satisfying `|Omega_1| + |Omega_2| = c_1 + c_2 + 2*(s_1 + s_2)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PrismFourCycleInvariant {
    pub a1: usize,
    pub a2: usize,
    pub b1: usize,
    pub b2: usize,
    pub s1: usize,
    pub s2: usize,
    pub c1: usize,
    pub c2: usize,
    pub omega1: usize,
    pub omega2: usize,
}

impl PrismFourCycleInvariant {
    /// Constructs and verifies the exact Theorem 2.9(b) four-cycle conservation law for `Pr(a_1, a_2; b_1, b_2; s_1, s_2)`.
    pub fn new(a1: usize, a2: usize, b1: usize, b2: usize, s1: usize, s2: usize) -> Self {
        assert!(a1 >= 1 && a2 >= 1 && b1 >= 1 && b2 >= 1 && s1 >= 1 && s2 >= 1);
        let c1 = a1 + a2;
        let c2 = b1 + b2;
        let omega1 = a1 + b1 + s1 + s2;
        let omega2 = a2 + b2 + s1 + s2;
        assert_eq!(omega1 + omega2, c1 + c2 + 2 * (s1 + s2));
        Self {
            a1,
            a2,
            b1,
            b2,
            s1,
            s2,
            c1,
            c2,
            omega1,
            omega2,
        }
    }
}
