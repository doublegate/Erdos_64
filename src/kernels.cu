// ============================================================================
// `kernels.cu` — NVRTC CUDA kernels for Erdős Problem #64 (Erdős–Gyárfás
// 2^k-Cycle Conjecture).
// ============================================================================
//
// Evaluates exact 4-cycle count `C_4(G)`, 8-cycle count `C_8(G)`, anchor-0
// 16-cycle count `C_16^{(0)}(G)`, and minimum degree `delta(G)` across batches
// of 64-bit adjacency-mask graphs (`[u64; 64]` per graph) with bit-exact parity
// against `verify_cpu::eval_dyadic_cycles_cpu`.

extern "C" __global__ void k_eval_dyadic_cycles(
    const unsigned long long* __restrict__ flat_adj,
    unsigned int n,
    unsigned int batch_size,
    unsigned long long* __restrict__ out_c4,
    unsigned long long* __restrict__ out_c8,
    unsigned long long* __restrict__ out_c16_anchor0,
    unsigned int* __restrict__ out_min_deg
) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= batch_size) {
        return;
    }

    const unsigned long long* adj = flat_adj + ((size_t)tid * 64ULL);

    // 1. Compute minimum degree delta(G) and maximum degree Delta(G)
    unsigned int min_deg = (n > 0) ? 64u : 0u;
    unsigned int max_deg = 0u;
    for (unsigned int u = 0; u < n; ++u) {
        unsigned int deg = (unsigned int)__popcll(adj[u]);
        if (deg < min_deg) min_deg = deg;
        if (deg > max_deg) max_deg = deg;
    }
    out_min_deg[tid] = min_deg;

    // 2. Count undirected 4-cycles C_4(G) in O(n^2) bitwise popcounts
    unsigned long long pairs_common = 0ULL;
    for (unsigned int u = 0; u < n; ++u) {
        unsigned long long au = adj[u];
        for (unsigned int v = u + 1; v < n; ++v) {
            unsigned long long c = (unsigned long long)__popcll(au & adj[v]);
            if (c >= 2ULL) {
                pairs_common += (c * (c - 1ULL)) >> 1;
            }
        }
    }
    out_c4[tid] = pairs_common >> 1;

    // 3. Count undirected 8-cycles C_8(G) via canonical-anchor bitmask path-extension
    unsigned long long total_c8 = 0ULL;
    for (unsigned int s = 0; s < n; ++s) {
        unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
        unsigned long long ns = adj[s] & higher_mask;
        if (__popcll(ns) < 2) {
            continue;
        }
        unsigned long long m1 = ns;
        while (m1 != 0ULL) {
            unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
            m1 &= (m1 - 1ULL);
            unsigned long long vis1 = (1ULL << s) | (1ULL << v1);
            unsigned long long valid_v7_base =
                (v1 == 63u) ? 0ULL : (ns & ~((1ULL << (v1 + 1u)) - 1ULL));
            if (valid_v7_base == 0ULL) {
                continue;
            }

            unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
            while (m2 != 0ULL) {
                unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                m2 &= (m2 - 1ULL);
                unsigned long long vis2 = vis1 | (1ULL << v2);

                unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                while (m3 != 0ULL) {
                    unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                    m3 &= (m3 - 1ULL);
                    unsigned long long vis3 = vis2 | (1ULL << v3);

                    unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                    while (m4 != 0ULL) {
                        unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                        m4 &= (m4 - 1ULL);
                        unsigned long long vis4 = vis3 | (1ULL << v4);

                        unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                        while (m5 != 0ULL) {
                            unsigned int v5 = (unsigned int)(__ffsll((long long)m5) - 1);
                            m5 &= (m5 - 1ULL);
                            unsigned long long vis5 = vis4 | (1ULL << v5);

                            unsigned long long m6 = adj[v5] & higher_mask & ~vis5;
                            while (m6 != 0ULL) {
                                unsigned int v6 = (unsigned int)(__ffsll((long long)m6) - 1);
                                m6 &= (m6 - 1ULL);
                                unsigned long long vis6 = vis5 | (1ULL << v6);
                                unsigned long long v7_mask = adj[v6] & valid_v7_base & ~vis6;
                                total_c8 += (unsigned long long)__popcll(v7_mask);
                            }
                        }
                    }
                }
            }
        }
    }
    out_c8[tid] = total_c8;

    // 4. Count 16-cycles anchored at s = 0 (0 = min V(C)) for sparse (Delta <= 3) graphs
    unsigned long long total_c16_0 = 0ULL;
    if (n >= 16u && max_deg <= 3u) {
        unsigned long long ns = adj[0] & ~1ULL;
        if (__popcll(ns) >= 2) {
            unsigned long long m1 = ns;
            unsigned long long masks[16];
            unsigned long long vis[16];
            while (m1 != 0ULL) {
                unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
                m1 &= (m1 - 1ULL);
                unsigned long long valid_v15 =
                    (v1 == 63u) ? 0ULL : (ns & ~((1ULL << (v1 + 1u)) - 1ULL));
                if (valid_v15 == 0ULL) {
                    continue;
                }
                vis[1] = 1ULL | (1ULL << v1);
                masks[1] = adj[v1] & ~vis[1];
                int depth = 1;
                while (depth >= 1) {
                    unsigned long long m = masks[depth];
                    if (m == 0ULL) {
                        depth -= 1;
                        continue;
                    }
                    unsigned int nxt = (unsigned int)(__ffsll((long long)m) - 1);
                    masks[depth] &= (m - 1ULL);
                    unsigned long long next_vis = vis[depth] | (1ULL << nxt);
                    if (depth == 13) {
                        unsigned long long v15_mask = adj[nxt] & valid_v15 & ~next_vis;
                        total_c16_0 += (unsigned long long)__popcll(v15_mask);
                    } else {
                        depth += 1;
                        vis[depth] = next_vis;
                        masks[depth] = adj[nxt] & ~next_vis;
                    }
                }
            }
        }
    }
    out_c16_anchor0[tid] = total_c16_0;
}

// ============================================================================
// Theorem 3.2(a): Bass-Ihara Chebyshev Non-Backtracking Trace Kernel
// Evaluates exact vertex-adjacency traces Tr(A^4), Tr(A^6), Tr(A^8) via
// symmetric half-power row-norm accumulation Tr(A^{2r}) = \sum_s ||A^r e_s||_2^2
// and computes the exact Hashimoto non-backtracking traces:
//   Tr(B^4) = Tr(A^4) - 8 * Tr(A^2) + 9 * n   (= 8 * C_4(G))
//   Tr(B^8) = Tr(A^8) - 16 * Tr(A^6) + 80 * Tr(A^4) - 128 * Tr(A^2) + 33 * n
//             (= 16 * C_8(G) whenever girth g(G) >= 5).
// ============================================================================
extern "C" __global__ void k_eval_hashimoto_trace_b4_b8(
    const unsigned long long* __restrict__ flat_adj,
    unsigned int n,
    unsigned int batch_size,
    unsigned long long* __restrict__ out_tr_a4,
    unsigned long long* __restrict__ out_tr_a6,
    unsigned long long* __restrict__ out_tr_a8,
    unsigned long long* __restrict__ out_tr_b4,
    unsigned long long* __restrict__ out_tr_b8
) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= batch_size) {
        return;
    }

    const unsigned long long* adj = flat_adj + ((size_t)tid * 64ULL);

    unsigned long long tr2 = 0ULL;
    unsigned long long tr4 = 0ULL;
    unsigned long long tr6 = 0ULL;
    unsigned long long tr8 = 0ULL;

    unsigned long long curr[64];
    unsigned long long next[64];

    for (unsigned int s = 0; s < n; ++s) {
        for (unsigned int u = 0; u < n; ++u) {
            curr[u] = (u == s) ? 1ULL : 0ULL;
        }
        for (unsigned int r = 1; r <= 4; ++r) {
            unsigned long long norm_sq = 0ULL;
            for (unsigned int u = 0; u < n; ++u) {
                unsigned long long sum = 0ULL;
                unsigned long long m = adj[u];
                while (m != 0ULL) {
                    unsigned int w = (unsigned int)(__ffsll((long long)m) - 1);
                    m &= (m - 1ULL);
                    sum += curr[w];
                }
                next[u] = sum;
                norm_sq += sum * sum;
            }
            for (unsigned int u = 0; u < n; ++u) {
                curr[u] = next[u];
            }
            if (r == 1u) tr2 += norm_sq;
            else if (r == 2u) tr4 += norm_sq;
            else if (r == 3u) tr6 += norm_sq;
            else if (r == 4u) tr8 += norm_sq;
        }
    }

    long long n_ll = (long long)n;
    long long tr2_ll = (long long)tr2;
    long long tr4_ll = (long long)tr4;
    long long tr6_ll = (long long)tr6;
    long long tr8_ll = (long long)tr8;

    long long b4 = tr4_ll - 8LL * tr2_ll + 9LL * n_ll;
    long long b8 = tr8_ll - 16LL * tr6_ll + 80LL * tr4_ll - 128LL * tr2_ll + 33LL * n_ll;

    out_tr_a4[tid] = tr4;
    out_tr_a6[tid] = tr6;
    out_tr_a8[tid] = tr8;
    out_tr_b4[tid] = (b4 > 0LL) ? (unsigned long long)b4 : 0ULL;
    out_tr_b8[tid] = (b8 > 0LL) ? (unsigned long long)b8 : 0ULL;
}

// ============================================================================
// Theorem 2.5(a): Half-Walk Meet-in-the-Middle Topological Classifier for Tr(B^16)
// Decomposes all closed non-backtracking walks of length 16 on a cubic graph
// (Delta(G) <= 3, n <= 64) into:
//   Tr(B^16) = 32*C_16(G) + 16*C_8(G) + 8*C_4(G) + N_16^{dumbbell}(G) + N_16^{theta}(G)
// with bit-exact parity against `verify_cpu::eval_n16_bicyclic_walks_cpu`.
// ============================================================================
struct HalfWalk8Gpu {
    unsigned long long v_mask;
    unsigned long long e_or_lo;
    unsigned long long e_xor_lo;
    unsigned int e_or_hi;
    unsigned int e_xor_hi;
    unsigned char first;
    unsigned char penult;
    short next;
};

extern "C" __global__ void k_enumerate_n16_bicyclic_walks(
    const unsigned long long* __restrict__ flat_adj,
    unsigned int n,
    unsigned int batch_size,
    unsigned long long* __restrict__ out_walks_c16,
    unsigned long long* __restrict__ out_walks_c8,
    unsigned long long* __restrict__ out_walks_c4,
    unsigned long long* __restrict__ out_walks_dumbbell,
    unsigned long long* __restrict__ out_walks_theta,
    unsigned long long* __restrict__ out_tr_b16_total
) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= batch_size) {
        return;
    }

    const unsigned long long* adj = flat_adj + ((size_t)tid * 64ULL);

    unsigned char edge_u[96];
    unsigned char edge_v[96];
    unsigned int edge_id_packed[64];
    unsigned int m_edges = 0u;

    for (unsigned int u = 0; u < n; ++u) {
        edge_id_packed[u] = 0u;
    }

    for (unsigned int u = 0; u < n; ++u) {
        unsigned long long m = adj[u];
        unsigned int shift = 0u;
        unsigned int packed = 0u;
        while (m != 0ULL) {
            unsigned int v = (unsigned int)(__ffsll((long long)m) - 1);
            m &= (m - 1ULL);
            unsigned int eid = 0u;
            if (u < v) {
                eid = m_edges;
                if (eid < 96u) {
                    edge_u[eid] = (unsigned char)u;
                    edge_v[eid] = (unsigned char)v;
                }
                m_edges += 1u;
            } else {
                unsigned int r_v = (unsigned int)__popcll(adj[v] & ((1ULL << u) - 1ULL));
                eid = (edge_id_packed[v] >> (8u * r_v)) & 0xFFu;
            }
            packed |= (eid & 0xFFu) << shift;
            shift += 8u;
        }
        edge_id_packed[u] = packed;
    }

    unsigned long long walks_c16 = 0ULL;
    unsigned long long walks_c8 = 0ULL;
    unsigned long long walks_c4 = 0ULL;
    unsigned long long walks_dumbbell = 0ULL;
    unsigned long long walks_theta = 0ULL;

    HalfWalk8Gpu hw[384];
    short head[64];
    unsigned long long h_adj[64];

    unsigned int st_u[9];
    unsigned int st_prev[9];
    unsigned long long st_m[9];
    unsigned int st_ep[9];
    unsigned long long st_vmask[9];
    unsigned long long st_eor_lo[9];
    unsigned int st_eor_hi[9];
    unsigned long long st_exor_lo[9];
    unsigned int st_exor_hi[9];
    unsigned char st_first[9];

    for (unsigned int s = 0; s < n; ++s) {
        for (unsigned int u = 0; u < n; ++u) {
            head[u] = -1;
        }
        unsigned int hw_cnt = 0u;

        st_u[0] = s;
        st_prev[0] = 64u;
        st_m[0] = adj[s];
        st_ep[0] = edge_id_packed[s];
        st_vmask[0] = 1ULL << s;
        st_eor_lo[0] = 0ULL;
        st_eor_hi[0] = 0u;
        st_exor_lo[0] = 0ULL;
        st_exor_hi[0] = 0u;
        st_first[0] = 0u;

        int depth = 0;
        while (depth >= 0) {
            unsigned long long m = st_m[depth];
            if (m == 0ULL) {
                depth -= 1;
                continue;
            }
            unsigned int nxt = (unsigned int)(__ffsll((long long)m) - 1);
            st_m[depth] &= (m - 1ULL);
            unsigned int eid = st_ep[depth] & 0xFFu;
            st_ep[depth] >>= 8u;

            if (nxt == st_prev[depth]) {
                continue;
            }

            unsigned long long n_vmask = st_vmask[depth] | (1ULL << nxt);
            unsigned long long n_eor_lo = st_eor_lo[depth];
            unsigned int n_eor_hi = st_eor_hi[depth];
            unsigned long long n_exor_lo = st_exor_lo[depth];
            unsigned int n_exor_hi = st_exor_hi[depth];
            if (eid < 64u) {
                unsigned long long bit = 1ULL << eid;
                n_eor_lo |= bit;
                n_exor_lo ^= bit;
            } else {
                unsigned int bit = 1u << (eid - 64u);
                n_eor_hi |= bit;
                n_exor_hi ^= bit;
            }
            unsigned char first = (depth == 0) ? (unsigned char)nxt : st_first[depth];

            if (depth == 7) {
                if (hw_cnt < 384u) {
                    hw[hw_cnt].v_mask = n_vmask;
                    hw[hw_cnt].e_or_lo = n_eor_lo;
                    hw[hw_cnt].e_xor_lo = n_exor_lo;
                    hw[hw_cnt].e_or_hi = n_eor_hi;
                    hw[hw_cnt].e_xor_hi = n_exor_hi;
                    hw[hw_cnt].first = first;
                    hw[hw_cnt].penult = (unsigned char)st_u[7];
                    hw[hw_cnt].next = head[nxt];
                    head[nxt] = (short)hw_cnt;
                    hw_cnt += 1u;
                }
            } else {
                int nd = depth + 1;
                st_u[nd] = nxt;
                st_prev[nd] = st_u[depth];
                st_m[nd] = adj[nxt];
                st_ep[nd] = edge_id_packed[nxt];
                st_vmask[nd] = n_vmask;
                st_eor_lo[nd] = n_eor_lo;
                st_eor_hi[nd] = n_eor_hi;
                st_exor_lo[nd] = n_exor_lo;
                st_exor_hi[nd] = n_exor_hi;
                st_first[nd] = first;
                depth = nd;
            }
        }

        for (unsigned int mid = 0; mid < n; ++mid) {
            short i = head[mid];
            while (i >= 0) {
                HalfWalk8Gpu hwi = hw[i];
                short j = hwi.next;
                while (j >= 0) {
                    HalfWalk8Gpu hwj = hw[j];
                    if (hwi.first != hwj.first && hwi.penult != hwj.penult) {
                        unsigned long long v_mask = hwi.v_mask | hwj.v_mask;
                        unsigned long long e_or_lo = hwi.e_or_lo | hwj.e_or_lo;
                        unsigned int e_or_hi = hwi.e_or_hi | hwj.e_or_hi;
                        int v_cnt = __popcll(v_mask);
                        int e_cnt = __popcll(e_or_lo) + __popc(e_or_hi);
                        int b1 = e_cnt - v_cnt + 1;

                        if (b1 == 1) {
                            if (v_cnt == 16 && e_cnt == 16) {
                                walks_c16 += 2ULL;
                            } else if (v_cnt == 8 && e_cnt == 8) {
                                walks_c8 += 2ULL;
                            } else if (v_cnt == 4 && e_cnt == 4) {
                                walks_c4 += 2ULL;
                            }
                        } else {
                            unsigned long long e_odd_lo = hwi.e_xor_lo ^ hwj.e_xor_lo;
                            unsigned int e_odd_hi = hwi.e_xor_hi ^ hwj.e_xor_hi;
                            unsigned long long e_even_lo = e_or_lo & ~e_odd_lo;
                            unsigned int e_even_hi = e_or_hi & ~e_odd_hi;

                            bool has_bridge = false;
                            if (e_even_lo != 0ULL || e_even_hi != 0u) {
                                unsigned long long vm = v_mask;
                                while (vm != 0ULL) {
                                    unsigned int u = (unsigned int)(__ffsll((long long)vm) - 1);
                                    vm &= (vm - 1ULL);
                                    h_adj[u] = 0ULL;
                                }
                                unsigned long long m_lo = e_or_lo;
                                while (m_lo != 0ULL) {
                                    unsigned int eid = (unsigned int)(__ffsll((long long)m_lo) - 1);
                                    m_lo &= (m_lo - 1ULL);
                                    unsigned int u = edge_u[eid];
                                    unsigned int v = edge_v[eid];
                                    h_adj[u] |= (1ULL << v);
                                    h_adj[v] |= (1ULL << u);
                                }
                                unsigned int m_hi = e_or_hi;
                                while (m_hi != 0u) {
                                    unsigned int eid = 64u + (unsigned int)(__ffs((int)m_hi) - 1);
                                    m_hi &= (m_hi - 1u);
                                    unsigned int u = edge_u[eid];
                                    unsigned int v = edge_v[eid];
                                    h_adj[u] |= (1ULL << v);
                                    h_adj[v] |= (1ULL << u);
                                }

                                unsigned long long ev_lo = e_even_lo;
                                while (ev_lo != 0ULL && !has_bridge) {
                                    unsigned int eid = (unsigned int)(__ffsll((long long)ev_lo) - 1);
                                    ev_lo &= (ev_lo - 1ULL);
                                    unsigned int ue = edge_u[eid];
                                    unsigned int ve = edge_v[eid];
                                    h_adj[ue] &= ~(1ULL << ve);
                                    h_adj[ve] &= ~(1ULL << ue);
                                    unsigned long long vis = 1ULL << ue;
                                    unsigned long long front = vis;
                                    unsigned long long target = 1ULL << ve;
                                    while (front != 0ULL && (vis & target) == 0ULL) {
                                        unsigned long long nxt_f = 0ULL;
                                        unsigned long long fm = front;
                                        while (fm != 0ULL) {
                                            unsigned int x = (unsigned int)(__ffsll((long long)fm) - 1);
                                            fm &= (fm - 1ULL);
                                            nxt_f |= h_adj[x];
                                        }
                                        front = nxt_f & ~vis;
                                        vis |= front;
                                    }
                                    if ((vis & target) == 0ULL) {
                                        has_bridge = true;
                                    }
                                    h_adj[ue] |= (1ULL << ve);
                                    h_adj[ve] |= (1ULL << ue);
                                }
                                unsigned int ev_hi = e_even_hi;
                                while (ev_hi != 0u && !has_bridge) {
                                    unsigned int eid = 64u + (unsigned int)(__ffs((int)ev_hi) - 1);
                                    ev_hi &= (ev_hi - 1u);
                                    unsigned int ue = edge_u[eid];
                                    unsigned int ve = edge_v[eid];
                                    h_adj[ue] &= ~(1ULL << ve);
                                    h_adj[ve] &= ~(1ULL << ue);
                                    unsigned long long vis = 1ULL << ue;
                                    unsigned long long front = vis;
                                    unsigned long long target = 1ULL << ve;
                                    while (front != 0ULL && (vis & target) == 0ULL) {
                                        unsigned long long nxt_f = 0ULL;
                                        unsigned long long fm = front;
                                        while (fm != 0ULL) {
                                            unsigned int x = (unsigned int)(__ffsll((long long)fm) - 1);
                                            fm &= (fm - 1ULL);
                                            nxt_f |= h_adj[x];
                                        }
                                        front = nxt_f & ~vis;
                                        vis |= front;
                                    }
                                    if ((vis & target) == 0ULL) {
                                        has_bridge = true;
                                    }
                                    h_adj[ue] |= (1ULL << ve);
                                    h_adj[ve] |= (1ULL << ue);
                                }
                            }

                            if (has_bridge) {
                                walks_dumbbell += 2ULL;
                            } else {
                                walks_theta += 2ULL;
                            }
                        }
                    }
                    j = hwj.next;
                }
                i = hwi.next;
            }
        }
    }

    out_walks_c16[tid] = walks_c16;
    out_walks_c8[tid] = walks_c8;
    out_walks_c4[tid] = walks_c4;
    out_walks_dumbbell[tid] = walks_dumbbell;
    out_walks_theta[tid] = walks_theta;
    out_tr_b16_total[tid] = walks_c16 + walks_c8 + walks_c4 + walks_dumbbell + walks_theta;
}

// ============================================================================
// Theorem 2.6: Order-10 Hashimoto Trace Tr(B^10), 5-Step Half-Path MITM C_10(G),
// and 3x3 Shifted Stieltjes / Hausdorff Hankel Moment Pinching Bounds on C_10(G)
// Matches `verify_cpu::eval_b10_hankel_pinch_cpu` bit-for-bit.
// ============================================================================
struct HalfPath5Gpu {
    unsigned long long int_mask;
    unsigned char v1;
    short next;
};

extern "C" __global__ void k_eval_b10_hankel_pinch(
    const unsigned long long* __restrict__ flat_adj,
    unsigned int n,
    unsigned int batch_size,
    unsigned long long* __restrict__ out_c6,
    unsigned long long* __restrict__ out_c10,
    unsigned long long* __restrict__ out_tr_a10,
    unsigned long long* __restrict__ out_tr_b10,
    unsigned long long* __restrict__ out_c10_hankel_lower,
    unsigned long long* __restrict__ out_c10_hankel_upper
) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= batch_size) {
        return;
    }

    const unsigned long long* adj = flat_adj + ((size_t)tid * 64ULL);

    // 1. Compute Tr(A^2), Tr(A^4), Tr(A^6), Tr(A^8), Tr(A^10) via Tr(A^{2r}) = sum_s ||A^r e_s||_2^2
    unsigned long long tr2 = 0ULL;
    unsigned long long tr4 = 0ULL;
    unsigned long long tr6 = 0ULL;
    unsigned long long tr8 = 0ULL;
    unsigned long long tr10 = 0ULL;

    unsigned long long curr[64];
    unsigned long long next[64];

    for (unsigned int s = 0; s < n; ++s) {
        for (unsigned int u = 0; u < n; ++u) {
            curr[u] = (u == s) ? 1ULL : 0ULL;
        }
        for (unsigned int r = 1; r <= 5; ++r) {
            unsigned long long norm_sq = 0ULL;
            for (unsigned int u = 0; u < n; ++u) {
                unsigned long long sum = 0ULL;
                unsigned long long m = adj[u];
                while (m != 0ULL) {
                    unsigned int w = (unsigned int)(__ffsll((long long)m) - 1);
                    m &= (m - 1ULL);
                    sum += curr[w];
                }
                next[u] = sum;
                norm_sq += sum * sum;
            }
            for (unsigned int u = 0; u < n; ++u) {
                curr[u] = next[u];
            }
            if (r == 1u) tr2 += norm_sq;
            else if (r == 2u) tr4 += norm_sq;
            else if (r == 3u) tr6 += norm_sq;
            else if (r == 4u) tr8 += norm_sq;
            else if (r == 5u) tr10 += norm_sq;
        }
    }

    unsigned long long n64 = (unsigned long long)n;
    unsigned long long base_tr6 = 87ULL * n64;
    unsigned long long c6 = (tr6 > base_tr6) ? ((tr6 - base_tr6) / 12ULL) : 0ULL;

    unsigned long long b10_pos = tr10 + 140ULL * tr6 + 400ULL * tr2;
    unsigned long long b10_neg = 20ULL * tr8 + 400ULL * tr4 + 63ULL * n64;
    unsigned long long tr_b10 = (b10_pos > b10_neg) ? (b10_pos - b10_neg) : 0ULL;

    // 2. Count exact simple 10-cycles C_10(G) via canonical-anchor 5-step Half-Path MITM
    unsigned long long total_c10 = 0ULL;
    HalfPath5Gpu hp[64];
    short head[64];

    if (n >= 10u) {
        for (unsigned int s = 0; s < n; ++s) {
            unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
            unsigned long long ns = adj[s] & higher_mask;
            if (__popcll(ns) < 2) {
                continue;
            }
            for (unsigned int u = 0; u < n; ++u) {
                head[u] = -1;
            }
            unsigned int hp_cnt = 0u;

            unsigned long long m1 = ns;
            while (m1 != 0ULL) {
                unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
                m1 &= (m1 - 1ULL);
                unsigned long long vis1 = (1ULL << v1);

                unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis1 | (1ULL << v2);

                    unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                    while (m3 != 0ULL) {
                        unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                        m3 &= (m3 - 1ULL);
                        unsigned long long vis3 = vis2 | (1ULL << v3);

                        unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                        while (m4 != 0ULL) {
                            unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                            m4 &= (m4 - 1ULL);
                            unsigned long long vis4 = vis3 | (1ULL << v4);

                            unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                            while (m5 != 0ULL) {
                                unsigned int mid = (unsigned int)(__ffsll((long long)m5) - 1);
                                m5 &= (m5 - 1ULL);
                                if (hp_cnt < 64u) {
                                    hp[hp_cnt].int_mask = vis4;
                                    hp[hp_cnt].v1 = (unsigned char)v1;
                                    hp[hp_cnt].next = head[mid];
                                    head[mid] = (short)hp_cnt;
                                    hp_cnt += 1u;
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int mid = s + 1u; mid < n; ++mid) {
                short i = head[mid];
                while (i >= 0) {
                    HalfPath5Gpu hpi = hp[i];
                    short j = hpi.next;
                    while (j >= 0) {
                        HalfPath5Gpu hpj = hp[j];
                        if (hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0ULL) {
                            total_c10 += 1ULL;
                        }
                        j = hpj.next;
                    }
                    i = hpi.next;
                }
            }
        }
    }

    // 3. Check bipartiteness via bitwise 2-coloring BFS to set r_trivial in {1, 2}
    unsigned long long red = 0ULL;
    unsigned long long blue = 0ULL;
    bool is_bip = true;
    for (unsigned int s = 0; s < n && is_bip; ++s) {
        unsigned long long bit = 1ULL << s;
        if (((red | blue) & bit) != 0ULL) {
            continue;
        }
        red |= bit;
        unsigned long long front = bit;
        bool cur_is_red = true;
        while (front != 0ULL && is_bip) {
            unsigned long long nxt = 0ULL;
            unsigned long long fm = front;
            while (fm != 0ULL) {
                unsigned int u = (unsigned int)(__ffsll((long long)fm) - 1);
                fm &= (fm - 1ULL);
                nxt |= adj[u];
            }
            if (cur_is_red) {
                if ((nxt & red) != 0ULL) {
                    is_bip = false;
                }
                front = nxt & ~blue;
                blue |= front;
            } else {
                if ((nxt & blue) != 0ULL) {
                    is_bip = false;
                }
                front = nxt & ~red;
                red |= front;
            }
            cur_is_red = !cur_is_red;
        }
    }

    // 4. Compute exact 3x3 shifted Stieltjes and Hausdorff Hankel bounds on C_10(G) in unsigned 64-bit arithmetic
    unsigned long long r64 = is_bip ? 2ULL : 1ULL;

    unsigned long long pos_lo = 4ULL * c6 * c6 * c6
        + 3ULL * n64 * c6 * c6
        + 252ULL * r64 * n64 * n64;
    unsigned long long neg_lo = 9ULL * n64 * n64 * c6
        + 4ULL * n64 * n64 * n64
        + 93ULL * r64 * c6 * c6
        + 9ULL * r64 * n64 * c6;
    unsigned long long den_pos = n64 * (c6 + n64);
    unsigned long long den_neg = 3ULL * r64 * (c6 + 5ULL * n64);
    unsigned long long c10_min = 0ULL;
    if (den_pos > den_neg && pos_lo > neg_lo) {
        unsigned long long num_u = 3ULL * (pos_lo - neg_lo);
        unsigned long long den_u = 5ULL * (den_pos - den_neg);
        c10_min = (num_u + den_u - 1ULL) / den_u;
    }

    unsigned long long pos_hi = 6ULL * c6 * c6 * c6 + 48ULL * n64 * n64 * n64;
    unsigned long long neg_hi = 27ULL * n64 * c6 * c6;
    unsigned long long c10_max = (2ULL * n64 > c6 && pos_hi >= neg_hi)
        ? ((pos_hi - neg_hi) / (5ULL * n64 * (2ULL * n64 - c6)))
        : 0ULL;

    out_c6[tid] = c6;
    out_c10[tid] = total_c10;
    out_tr_a10[tid] = tr10;
    out_tr_b10[tid] = tr_b10;
    out_c10_hankel_lower[tid] = c10_min;
    out_c10_hankel_upper[tid] = c10_max;
}

// ============================================================================
// Theorem 2.7: Exact Order-12 Hashimoto Trace Tr(B^12), 6-Cycle Overlap Counts
// tau_1(G) and tau_3(G), 6-Step Half-Path MITM C_12(G), and 3D Chebyshev-Bessel
// Projection on span{1, lambda^2, lambda^4}
// Matches `verify_cpu::eval_b12_deg2_bessel_pinch_cpu` bit-for-bit.
// ============================================================================
__device__ __forceinline__ void add_edge_bit_192(
    const unsigned long long* adj,
    unsigned int u,
    unsigned int v,
    unsigned long long& w0,
    unsigned long long& w1,
    unsigned long long& w2
) {
    unsigned int a = (u < v) ? u : v;
    unsigned int b = (u < v) ? v : u;
    unsigned long long lower_mask = (b == 0u) ? 0ULL : ((1ULL << b) - 1ULL);
    unsigned int eid = 3u * a + (unsigned int)__popcll(adj[a] & lower_mask);
    if (eid < 64u) {
        w0 |= (1ULL << eid);
    } else if (eid < 128u) {
        w1 |= (1ULL << (eid - 64u));
    } else if (eid < 192u) {
        w2 |= (1ULL << (eid - 128u));
    }
}

extern "C" __global__ void k_eval_b12_deg2_bessel_pinch(
    const unsigned long long* __restrict__ flat_adj,
    unsigned int n,
    unsigned int batch_size,
    unsigned long long* __restrict__ out_c6,
    unsigned long long* __restrict__ out_c10,
    unsigned long long* __restrict__ out_c12,
    unsigned long long* __restrict__ out_tau1,
    unsigned long long* __restrict__ out_tau3,
    unsigned long long* __restrict__ out_tr_a12,
    unsigned long long* __restrict__ out_tr_b12,
    unsigned long long* __restrict__ out_tr_b16_3d_bessel_floor
) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= batch_size) {
        return;
    }

    const unsigned long long* adj = flat_adj + ((size_t)tid * 64ULL);

    // 1. Compute Tr(A^2), Tr(A^4), Tr(A^6), Tr(A^8), Tr(A^10), Tr(A^12) via Tr(A^{2r}) = sum_s ||A^r e_s||_2^2
    unsigned long long tr2 = 0ULL;
    unsigned long long tr4 = 0ULL;
    unsigned long long tr6 = 0ULL;
    unsigned long long tr8 = 0ULL;
    unsigned long long tr10 = 0ULL;
    unsigned long long tr12 = 0ULL;

    unsigned long long curr[64];
    unsigned long long next[64];

    for (unsigned int s = 0; s < n; ++s) {
        for (unsigned int u = 0; u < n; ++u) {
            curr[u] = (u == s) ? 1ULL : 0ULL;
        }
        for (unsigned int r = 1; r <= 6; ++r) {
            unsigned long long norm_sq = 0ULL;
            for (unsigned int u = 0; u < n; ++u) {
                unsigned long long sum = 0ULL;
                unsigned long long m = adj[u];
                while (m != 0ULL) {
                    unsigned int w = (unsigned int)(__ffsll((long long)m) - 1);
                    m &= (m - 1ULL);
                    sum += curr[w];
                }
                next[u] = sum;
                norm_sq += sum * sum;
            }
            for (unsigned int u = 0; u < n; ++u) {
                curr[u] = next[u];
            }
            if (r == 1u) tr2 += norm_sq;
            else if (r == 2u) tr4 += norm_sq;
            else if (r == 3u) tr6 += norm_sq;
            else if (r == 4u) tr8 += norm_sq;
            else if (r == 5u) tr10 += norm_sq;
            else if (r == 6u) tr12 += norm_sq;
        }
    }

    unsigned long long n64 = (unsigned long long)n;
    unsigned long long b12_pos = tr12 + 216ULL * tr8 + 1680ULL * tr4 + 129ULL * n64;
    unsigned long long b12_neg = 24ULL * tr10 + 896ULL * tr6 + 1152ULL * tr2;
    unsigned long long tr_b12 = (b12_pos > b12_neg) ? (b12_pos - b12_neg) : 0ULL;

    // 2. Enumerate simple 6-cycles and compute (c6, tau1, tau3) via 192-bit undirected edge masks
    unsigned long long c6_w0[128];
    unsigned long long c6_w1[128];
    unsigned long long c6_w2[128];
    unsigned int c6_cnt = 0u;

    for (unsigned int s = 0; s < n; ++s) {
        unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
        unsigned long long ns = adj[s] & higher_mask;
        unsigned long long m1 = ns;
        while (m1 != 0ULL) {
            unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
            m1 &= (m1 - 1ULL);
            unsigned long long m5_mask = (v1 == 63u) ? 0ULL : ~((1ULL << (v1 + 1u)) - 1ULL);
            unsigned long long m5 = ns & m5_mask;
            while (m5 != 0ULL) {
                unsigned int v5 = (unsigned int)(__ffsll((long long)m5) - 1);
                m5 &= (m5 - 1ULL);
                unsigned long long vis15 = (1ULL << s) | (1ULL << v1) | (1ULL << v5);
                unsigned long long m2 = adj[v1] & higher_mask & ~vis15;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis15 | (1ULL << v2);
                    unsigned long long m4 = adj[v5] & higher_mask & ~vis2;
                    while (m4 != 0ULL) {
                        unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                        m4 &= (m4 - 1ULL);
                        unsigned long long vis4 = vis2 | (1ULL << v4);
                        unsigned long long m3 = adj[v2] & adj[v4] & higher_mask & ~vis4;
                        while (m3 != 0ULL) {
                            unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                            m3 &= (m3 - 1ULL);
                            if (c6_cnt < 128u) {
                                unsigned long long w0 = 0ULL;
                                unsigned long long w1 = 0ULL;
                                unsigned long long w2 = 0ULL;
                                add_edge_bit_192(adj, s, v1, w0, w1, w2);
                                add_edge_bit_192(adj, v1, v2, w0, w1, w2);
                                add_edge_bit_192(adj, v2, v3, w0, w1, w2);
                                add_edge_bit_192(adj, v3, v4, w0, w1, w2);
                                add_edge_bit_192(adj, v4, v5, w0, w1, w2);
                                add_edge_bit_192(adj, v5, s, w0, w1, w2);
                                c6_w0[c6_cnt] = w0;
                                c6_w1[c6_cnt] = w1;
                                c6_w2[c6_cnt] = w2;
                            }
                            c6_cnt += 1u;
                        }
                    }
                }
            }
        }
    }

    unsigned int stored_c6 = (c6_cnt < 128u) ? c6_cnt : 128u;
    unsigned long long tau1 = 0ULL;
    unsigned long long tau3 = 0ULL;
    for (unsigned int a = 0; a < stored_c6; ++a) {
        for (unsigned int b = a + 1u; b < stored_c6; ++b) {
            unsigned int shared = (unsigned int)(
                __popcll(c6_w0[a] & c6_w0[b])
                + __popcll(c6_w1[a] & c6_w1[b])
                + __popcll(c6_w2[a] & c6_w2[b])
            );
            if (shared == 1u) {
                tau1 += 1ULL;
            } else if (shared == 3u) {
                tau3 += 1ULL;
            }
        }
    }
    unsigned long long c6 = (unsigned long long)c6_cnt;

    // 3. Count exact simple 10-cycles C_10(G) and 12-cycles C_12(G) via Half-Path MITM
    unsigned long long total_c10 = 0ULL;
    unsigned long long total_c12 = 0ULL;
    HalfPath5Gpu hp[128];
    short head[64];

    if (n >= 10u) {
        for (unsigned int s = 0; s < n; ++s) {
            unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
            unsigned long long ns = adj[s] & higher_mask;
            if (__popcll(ns) < 2) {
                continue;
            }
            for (unsigned int u = 0; u < n; ++u) {
                head[u] = -1;
            }
            unsigned int hp_cnt = 0u;

            unsigned long long m1 = ns;
            while (m1 != 0ULL) {
                unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
                m1 &= (m1 - 1ULL);
                unsigned long long vis1 = (1ULL << v1);

                unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis1 | (1ULL << v2);

                    unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                    while (m3 != 0ULL) {
                        unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                        m3 &= (m3 - 1ULL);
                        unsigned long long vis3 = vis2 | (1ULL << v3);

                        unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                        while (m4 != 0ULL) {
                            unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                            m4 &= (m4 - 1ULL);
                            unsigned long long vis4 = vis3 | (1ULL << v4);

                            unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                            while (m5 != 0ULL) {
                                unsigned int mid = (unsigned int)(__ffsll((long long)m5) - 1);
                                m5 &= (m5 - 1ULL);
                                if (hp_cnt < 64u) {
                                    hp[hp_cnt].int_mask = vis4;
                                    hp[hp_cnt].v1 = (unsigned char)v1;
                                    hp[hp_cnt].next = head[mid];
                                    head[mid] = (short)hp_cnt;
                                    hp_cnt += 1u;
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int mid = s + 1u; mid < n; ++mid) {
                short i = head[mid];
                while (i >= 0) {
                    HalfPath5Gpu hpi = hp[i];
                    short j = hpi.next;
                    while (j >= 0) {
                        HalfPath5Gpu hpj = hp[j];
                        if (hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0ULL) {
                            total_c10 += 1ULL;
                        }
                        j = hpj.next;
                    }
                    i = hpi.next;
                }
            }
        }
    }

    if (n >= 12u) {
        for (unsigned int s = 0; s < n; ++s) {
            unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
            unsigned long long ns = adj[s] & higher_mask;
            if (__popcll(ns) < 2) {
                continue;
            }
            for (unsigned int u = 0; u < n; ++u) {
                head[u] = -1;
            }
            unsigned int hp_cnt = 0u;

            unsigned long long m1 = ns;
            while (m1 != 0ULL) {
                unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
                m1 &= (m1 - 1ULL);
                unsigned long long vis1 = (1ULL << v1);

                unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis1 | (1ULL << v2);

                    unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                    while (m3 != 0ULL) {
                        unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                        m3 &= (m3 - 1ULL);
                        unsigned long long vis3 = vis2 | (1ULL << v3);

                        unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                        while (m4 != 0ULL) {
                            unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                            m4 &= (m4 - 1ULL);
                            unsigned long long vis4 = vis3 | (1ULL << v4);

                            unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                            while (m5 != 0ULL) {
                                unsigned int v5 = (unsigned int)(__ffsll((long long)m5) - 1);
                                m5 &= (m5 - 1ULL);
                                unsigned long long vis5 = vis4 | (1ULL << v5);

                                unsigned long long m6 = adj[v5] & higher_mask & ~vis5;
                                while (m6 != 0ULL) {
                                    unsigned int mid = (unsigned int)(__ffsll((long long)m6) - 1);
                                    m6 &= (m6 - 1ULL);
                                    if (hp_cnt < 128u) {
                                        hp[hp_cnt].int_mask = vis5;
                                        hp[hp_cnt].v1 = (unsigned char)v1;
                                        hp[hp_cnt].next = head[mid];
                                        head[mid] = (short)hp_cnt;
                                        hp_cnt += 1u;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int mid = s + 1u; mid < n; ++mid) {
                short i = head[mid];
                while (i >= 0) {
                    HalfPath5Gpu hpi = hp[i];
                    short j = hpi.next;
                    while (j >= 0) {
                        HalfPath5Gpu hpj = hp[j];
                        if (hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0ULL) {
                            total_c12 += 1ULL;
                        }
                        j = hpj.next;
                    }
                    i = hpi.next;
                }
            }
        }
    }

    // 4. Check bipartiteness via bitwise 2-coloring BFS to set r_trivial in {1, 2}
    unsigned long long red = 0ULL;
    unsigned long long blue = 0ULL;
    bool is_bip = true;
    for (unsigned int s = 0; s < n && is_bip; ++s) {
        unsigned long long bit = 1ULL << s;
        if (((red | blue) & bit) != 0ULL) {
            continue;
        }
        red |= bit;
        unsigned long long front = bit;
        bool cur_is_red = true;
        while (front != 0ULL && is_bip) {
            unsigned long long nxt = 0ULL;
            unsigned long long fm = front;
            while (fm != 0ULL) {
                unsigned int u = (unsigned int)(__ffsll((long long)fm) - 1);
                fm &= (fm - 1ULL);
                nxt |= adj[u];
            }
            if (cur_is_red) {
                if ((nxt & red) != 0ULL) {
                    is_bip = false;
                }
                front = nxt & ~blue;
                blue |= front;
            } else {
                if ((nxt & blue) != 0ULL) {
                    is_bip = false;
                }
                front = nxt & ~red;
                red |= front;
            }
            cur_is_red = !cur_is_red;
        }
    }

    // 5. Compute exact 3D Chebyshev-Bessel floor on Tr(B^16) in pure unsigned 64-bit arithmetic
    unsigned long long r64 = is_bip ? 2ULL : 1ULL;
    unsigned long long tr_b16_3d = 0ULL;
    if (n64 > 7ULL * r64 && 11179ULL * r64 > 85ULL * n64) {
        unsigned long long nmr = n64 - r64;
        unsigned long long nm7r = n64 - 7ULL * r64;
        unsigned long long den_2d = 3ULL * n64 * nmr * nm7r;
        unsigned long long z_val = 10ULL * total_c10 + 24ULL * c6;
        unsigned long long lhs_sq = z_val * nmr;
        unsigned long long rhs_sq = 3ULL * n64 * (n64 + 257ULL * r64);
        unsigned long long abs_sq = (lhs_sq >= rhs_sq) ? (lhs_sq - rhs_sq) : (rhs_sq - lhs_sq);
        unsigned long long num_2d = 18ULL * n64 * n64 * (11179ULL * r64 - 85ULL * n64) * nm7r
            + 2ULL * abs_sq * abs_sq;

        unsigned long long d3_pos = n64 * (n64 + 17ULL * r64) * c6 + n64 * n64 * n64;
        unsigned long long d3_neg = nmr * c6 * c6 + 31ULL * r64 * n64 * n64;
        if (d3_pos > d3_neg) {
            unsigned long long d3_tilde = d3_pos - d3_neg;
            unsigned long long m2_pos = n64 * nm7r * tr_b12
                + 2ULL * n64 * (n64 + 17ULL * r64) * z_val
                + 12ULL * n64 * c6 * (n64 + 257ULL * r64);
            unsigned long long m2_neg = 4ULL * nmr * c6 * z_val
                + 24ULL * n64 * n64 * (n64 + 257ULL * r64);
            unsigned long long abs_m2 = (m2_pos >= m2_neg) ? (m2_pos - m2_neg) : (m2_neg - m2_pos);
            unsigned long long num_3d = 8ULL * d3_tilde * num_2d + nmr * (abs_m2 * abs_m2);
            unsigned long long den_3d = 8ULL * d3_tilde * den_2d;
            tr_b16_3d = (num_3d + den_3d - 1ULL) / den_3d;
        } else {
            tr_b16_3d = (num_2d + den_2d - 1ULL) / den_2d;
        }
    }

    out_c6[tid] = c6;
    out_c10[tid] = total_c10;
    out_c12[tid] = total_c12;
    out_tau1[tid] = tau1;
    out_tau3[tid] = tau3;
    out_tr_a12[tid] = tr12;
    out_tr_b12[tid] = tr_b12;
    out_tr_b16_3d_bessel_floor[tid] = tr_b16_3d;
}

__device__ __forceinline__ long long det4x4_ll(const long long m[4][4]) {
    long long p01 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    long long p02 = m[0][0] * m[1][2] - m[0][2] * m[1][0];
    long long p03 = m[0][0] * m[1][3] - m[0][3] * m[1][0];
    long long p12 = m[0][1] * m[1][2] - m[0][2] * m[1][1];
    long long p13 = m[0][1] * m[1][3] - m[0][3] * m[1][1];
    long long p23 = m[0][2] * m[1][3] - m[0][3] * m[1][2];

    long long q01 = m[2][0] * m[3][1] - m[2][1] * m[3][0];
    long long q02 = m[2][0] * m[3][2] - m[2][2] * m[3][0];
    long long q03 = m[2][0] * m[3][3] - m[2][3] * m[3][0];
    long long q12 = m[2][1] * m[3][2] - m[2][2] * m[3][1];
    long long q13 = m[2][1] * m[3][3] - m[2][3] * m[3][1];
    long long q23 = m[2][2] * m[3][3] - m[2][3] * m[3][2];

    return p01 * q23 - p02 * q13 + p03 * q12 + p12 * q03 - p13 * q02 + p23 * q01;
}

extern "C" __global__ void k_eval_b14_deg3_bessel_pinch(
    const unsigned long long* __restrict__ flat_adj,
    unsigned int n,
    unsigned int batch_size,
    unsigned long long* __restrict__ out_c6,
    unsigned long long* __restrict__ out_c10,
    unsigned long long* __restrict__ out_c12,
    unsigned long long* __restrict__ out_c14,
    unsigned long long* __restrict__ out_beta1,
    unsigned long long* __restrict__ out_tr_a14,
    unsigned long long* __restrict__ out_tr_b14,
    unsigned long long* __restrict__ out_tr_b16_4d_bessel_floor
) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= batch_size) {
        return;
    }

    const unsigned long long* adj = flat_adj + ((size_t)tid * 64ULL);

    // 1. Compute Tr(A^2)..Tr(A^14) via Tr(A^{2r}) = sum_s ||A^r e_s||_2^2
    unsigned long long tr2 = 0ULL;
    unsigned long long tr4 = 0ULL;
    unsigned long long tr6 = 0ULL;
    unsigned long long tr8 = 0ULL;
    unsigned long long tr10 = 0ULL;
    unsigned long long tr12 = 0ULL;
    unsigned long long tr14 = 0ULL;

    unsigned long long curr[64];
    unsigned long long next[64];

    for (unsigned int s = 0; s < n; ++s) {
        for (unsigned int u = 0; u < n; ++u) {
            curr[u] = (u == s) ? 1ULL : 0ULL;
        }
        for (unsigned int r = 1; r <= 7; ++r) {
            unsigned long long norm_sq = 0ULL;
            for (unsigned int u = 0; u < n; ++u) {
                unsigned long long sum = 0ULL;
                unsigned long long m = adj[u];
                while (m != 0ULL) {
                    unsigned int w = (unsigned int)(__ffsll((long long)m) - 1);
                    m &= (m - 1ULL);
                    sum += curr[w];
                }
                next[u] = sum;
                norm_sq += sum * sum;
            }
            for (unsigned int u = 0; u < n; ++u) {
                curr[u] = next[u];
            }
            if (r == 1u) tr2 += norm_sq;
            else if (r == 2u) tr4 += norm_sq;
            else if (r == 3u) tr6 += norm_sq;
            else if (r == 4u) tr8 += norm_sq;
            else if (r == 5u) tr10 += norm_sq;
            else if (r == 6u) tr12 += norm_sq;
            else if (r == 7u) tr14 += norm_sq;
        }
    }

    unsigned long long n64 = (unsigned long long)n;
    unsigned long long b12_pos = tr12 + 216ULL * tr8 + 1680ULL * tr4 + 129ULL * n64;
    unsigned long long b12_neg = 24ULL * tr10 + 896ULL * tr6 + 1152ULL * tr2;
    unsigned long long tr_b12 = (b12_pos > b12_neg) ? (b12_pos - b12_neg) : 0ULL;

    unsigned long long b14_pos = tr14 + 308ULL * tr10 + 4704ULL * tr6 + 3136ULL * tr2;
    unsigned long long b14_neg = 28ULL * tr12 + 1680ULL * tr8 + 6272ULL * tr4 + 255ULL * n64;
    unsigned long long tr_b14 = (b14_pos > b14_neg) ? (b14_pos - b14_neg) : 0ULL;

    // 2. Enumerate simple 6-cycles and accumulate hair counts kappa[3*u + slot] for (6,1,6)-barbells beta_1(G)
    unsigned short kappa[192];
    for (unsigned int i = 0; i < 192u; ++i) {
        kappa[i] = 0;
    }
    unsigned long long c6 = 0ULL;

    for (unsigned int s = 0; s < n; ++s) {
        unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
        unsigned long long ns = adj[s] & higher_mask;
        unsigned long long m1 = ns;
        while (m1 != 0ULL) {
            unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
            m1 &= (m1 - 1ULL);
            unsigned long long m5_mask = (v1 == 63u) ? 0ULL : ~((1ULL << (v1 + 1u)) - 1ULL);
            unsigned long long m5 = ns & m5_mask;
            while (m5 != 0ULL) {
                unsigned int v5 = (unsigned int)(__ffsll((long long)m5) - 1);
                m5 &= (m5 - 1ULL);
                unsigned long long vis15 = (1ULL << s) | (1ULL << v1) | (1ULL << v5);
                unsigned long long m2 = adj[v1] & higher_mask & ~vis15;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis15 | (1ULL << v2);
                    unsigned long long m4 = adj[v5] & higher_mask & ~vis2;
                    while (m4 != 0ULL) {
                        unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                        m4 &= (m4 - 1ULL);
                        unsigned long long vis4 = vis2 | (1ULL << v4);
                        unsigned long long m3 = adj[v2] & adj[v4] & higher_mask & ~vis4;
                        while (m3 != 0ULL) {
                            unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                            m3 &= (m3 - 1ULL);
                            c6 += 1ULL;
                            unsigned int cyc[6] = {s, v1, v2, v3, v4, v5};
                            for (unsigned int k = 0; k < 6u; ++k) {
                                unsigned int u = cyc[k];
                                unsigned int vp = cyc[(k + 5u) % 6u];
                                unsigned int vn = cyc[(k + 1u) % 6u];
                                unsigned long long rem = adj[u] & ~((1ULL << vp) | (1ULL << vn));
                                while (rem != 0ULL) {
                                    unsigned int w = (unsigned int)(__ffsll((long long)rem) - 1);
                                    rem &= (rem - 1ULL);
                                    unsigned int slot = (unsigned int)__popcll(adj[u] & ((1ULL << w) - 1ULL));
                                    unsigned int idx = 3u * u + slot;
                                    if (idx < 192u) {
                                        kappa[idx] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    unsigned long long beta1 = 0ULL;
    for (unsigned int u = 0; u < n; ++u) {
        unsigned long long higher_u = (u == 63u) ? 0ULL : ~((1ULL << (u + 1u)) - 1ULL);
        unsigned long long m = adj[u] & higher_u;
        while (m != 0ULL) {
            unsigned int v = (unsigned int)(__ffsll((long long)m) - 1);
            m &= (m - 1ULL);
            unsigned int su = (unsigned int)__popcll(adj[u] & ((1ULL << v) - 1ULL));
            unsigned int sv = (unsigned int)__popcll(adj[v] & ((1ULL << u) - 1ULL));
            unsigned int idx_u = 3u * u + su;
            unsigned int idx_v = 3u * v + sv;
            if (idx_u < 192u && idx_v < 192u) {
                beta1 += (unsigned long long)kappa[idx_u] * (unsigned long long)kappa[idx_v];
            }
        }
    }

    // 3. Count exact simple C_10, C_12, C_14 via Half-Path MITM
    unsigned long long total_c10 = 0ULL;
    unsigned long long total_c12 = 0ULL;
    unsigned long long total_c14 = 0ULL;
    HalfPath5Gpu hp[256];
    short head[64];

    if (n >= 10u) {
        for (unsigned int s = 0; s < n; ++s) {
            unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
            unsigned long long ns = adj[s] & higher_mask;
            if (__popcll(ns) < 2) {
                continue;
            }
            for (unsigned int u = 0; u < n; ++u) {
                head[u] = -1;
            }
            unsigned int hp_cnt = 0u;

            unsigned long long m1 = ns;
            while (m1 != 0ULL) {
                unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
                m1 &= (m1 - 1ULL);
                unsigned long long vis1 = (1ULL << v1);

                unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis1 | (1ULL << v2);

                    unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                    while (m3 != 0ULL) {
                        unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                        m3 &= (m3 - 1ULL);
                        unsigned long long vis3 = vis2 | (1ULL << v3);

                        unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                        while (m4 != 0ULL) {
                            unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                            m4 &= (m4 - 1ULL);
                            unsigned long long vis4 = vis3 | (1ULL << v4);

                            unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                            while (m5 != 0ULL) {
                                unsigned int mid = (unsigned int)(__ffsll((long long)m5) - 1);
                                m5 &= (m5 - 1ULL);
                                if (hp_cnt < 64u) {
                                    hp[hp_cnt].int_mask = vis4;
                                    hp[hp_cnt].v1 = (unsigned char)v1;
                                    hp[hp_cnt].next = head[mid];
                                    head[mid] = (short)hp_cnt;
                                    hp_cnt += 1u;
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int mid = s + 1u; mid < n; ++mid) {
                short i = head[mid];
                while (i >= 0) {
                    HalfPath5Gpu hpi = hp[i];
                    short j = hpi.next;
                    while (j >= 0) {
                        HalfPath5Gpu hpj = hp[j];
                        if (hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0ULL) {
                            total_c10 += 1ULL;
                        }
                        j = hpj.next;
                    }
                    i = hpi.next;
                }
            }
        }
    }

    if (n >= 12u) {
        for (unsigned int s = 0; s < n; ++s) {
            unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
            unsigned long long ns = adj[s] & higher_mask;
            if (__popcll(ns) < 2) {
                continue;
            }
            for (unsigned int u = 0; u < n; ++u) {
                head[u] = -1;
            }
            unsigned int hp_cnt = 0u;

            unsigned long long m1 = ns;
            while (m1 != 0ULL) {
                unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
                m1 &= (m1 - 1ULL);
                unsigned long long vis1 = (1ULL << v1);

                unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis1 | (1ULL << v2);

                    unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                    while (m3 != 0ULL) {
                        unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                        m3 &= (m3 - 1ULL);
                        unsigned long long vis3 = vis2 | (1ULL << v3);

                        unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                        while (m4 != 0ULL) {
                            unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                            m4 &= (m4 - 1ULL);
                            unsigned long long vis4 = vis3 | (1ULL << v4);

                            unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                            while (m5 != 0ULL) {
                                unsigned int v5 = (unsigned int)(__ffsll((long long)m5) - 1);
                                m5 &= (m5 - 1ULL);
                                unsigned long long vis5 = vis4 | (1ULL << v5);

                                unsigned long long m6 = adj[v5] & higher_mask & ~vis5;
                                while (m6 != 0ULL) {
                                    unsigned int mid = (unsigned int)(__ffsll((long long)m6) - 1);
                                    m6 &= (m6 - 1ULL);
                                    if (hp_cnt < 128u) {
                                        hp[hp_cnt].int_mask = vis5;
                                        hp[hp_cnt].v1 = (unsigned char)v1;
                                        hp[hp_cnt].next = head[mid];
                                        head[mid] = (short)hp_cnt;
                                        hp_cnt += 1u;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int mid = s + 1u; mid < n; ++mid) {
                short i = head[mid];
                while (i >= 0) {
                    HalfPath5Gpu hpi = hp[i];
                    short j = hpi.next;
                    while (j >= 0) {
                        HalfPath5Gpu hpj = hp[j];
                        if (hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0ULL) {
                            total_c12 += 1ULL;
                        }
                        j = hpj.next;
                    }
                    i = hpi.next;
                }
            }
        }
    }

    if (n >= 14u) {
        for (unsigned int s = 0; s < n; ++s) {
            unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
            unsigned long long ns = adj[s] & higher_mask;
            if (__popcll(ns) < 2) {
                continue;
            }
            for (unsigned int u = 0; u < n; ++u) {
                head[u] = -1;
            }
            unsigned int hp_cnt = 0u;

            unsigned long long m1 = ns;
            while (m1 != 0ULL) {
                unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
                m1 &= (m1 - 1ULL);
                unsigned long long vis1 = (1ULL << v1);

                unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis1 | (1ULL << v2);

                    unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                    while (m3 != 0ULL) {
                        unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                        m3 &= (m3 - 1ULL);
                        unsigned long long vis3 = vis2 | (1ULL << v3);

                        unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                        while (m4 != 0ULL) {
                            unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                            m4 &= (m4 - 1ULL);
                            unsigned long long vis4 = vis3 | (1ULL << v4);

                            unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                            while (m5 != 0ULL) {
                                unsigned int v5 = (unsigned int)(__ffsll((long long)m5) - 1);
                                m5 &= (m5 - 1ULL);
                                unsigned long long vis5 = vis4 | (1ULL << v5);

                                unsigned long long m6 = adj[v5] & higher_mask & ~vis5;
                                while (m6 != 0ULL) {
                                    unsigned int v6 = (unsigned int)(__ffsll((long long)m6) - 1);
                                    m6 &= (m6 - 1ULL);
                                    unsigned long long vis6 = vis5 | (1ULL << v6);

                                    unsigned long long m7 = adj[v6] & higher_mask & ~vis6;
                                    while (m7 != 0ULL) {
                                        unsigned int mid = (unsigned int)(__ffsll((long long)m7) - 1);
                                        m7 &= (m7 - 1ULL);
                                        if (hp_cnt < 256u) {
                                            hp[hp_cnt].int_mask = vis6;
                                            hp[hp_cnt].v1 = (unsigned char)v1;
                                            hp[hp_cnt].next = head[mid];
                                            head[mid] = (short)hp_cnt;
                                            hp_cnt += 1u;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int mid = s + 1u; mid < n; ++mid) {
                short i = head[mid];
                while (i >= 0) {
                    HalfPath5Gpu hpi = hp[i];
                    short j = hpi.next;
                    while (j >= 0) {
                        HalfPath5Gpu hpj = hp[j];
                        if (hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0ULL) {
                            total_c14 += 1ULL;
                        }
                        j = hpj.next;
                    }
                    i = hpi.next;
                }
            }
        }
    }

    // 4. Check bipartiteness via bitwise 2-coloring BFS to set r_trivial in {1, 2}
    unsigned long long red = 0ULL;
    unsigned long long blue = 0ULL;
    bool is_bip = true;
    for (unsigned int s = 0; s < n && is_bip; ++s) {
        unsigned long long bit = 1ULL << s;
        if (((red | blue) & bit) != 0ULL) {
            continue;
        }
        red |= bit;
        unsigned long long front = bit;
        bool cur_is_red = true;
        while (front != 0ULL && is_bip) {
            unsigned long long nxt = 0ULL;
            unsigned long long fm = front;
            while (fm != 0ULL) {
                unsigned int u = (unsigned int)(__ffsll((long long)fm) - 1);
                fm &= (fm - 1ULL);
                nxt |= adj[u];
            }
            if (cur_is_red) {
                if ((nxt & red) != 0ULL) {
                    is_bip = false;
                }
                front = nxt & ~blue;
                blue |= front;
            } else {
                if ((nxt & blue) != 0ULL) {
                    is_bip = false;
                }
                front = nxt & ~red;
                red |= front;
            }
            cur_is_red = !cur_is_red;
        }
    }

    // 5. Compute exact 4D Chebyshev-Bessel floor on Tr(B^16) in pure 64-bit integer arithmetic
    unsigned long long r64 = is_bip ? 2ULL : 1ULL;
    unsigned long long tr_b16_4d = 0ULL;
    if (n64 > 7ULL * r64 && 11179ULL * r64 > 85ULL * n64) {
        long long ni = (long long)n64;
        long long ri = (long long)r64;
        long long c6i = (long long)c6;
        long long c10i = (long long)total_c10;
        long long tr_b12i = (long long)tr_b12;
        long long tr_b14i = (long long)tr_b14;

        unsigned long long den_2d = 3ULL * n64 * (n64 - r64) * (n64 - 7ULL * r64);
        long long term_sq = (10LL * c10i + 24LL * c6i) * (ni - ri) - 3LL * ni * (ni + 257LL * ri);
        unsigned long long abs_sq = (unsigned long long)(term_sq >= 0LL ? term_sq : -term_sq);
        unsigned long long num_2d = 18ULL * n64 * n64 * (11179ULL * r64 - 85ULL * n64) * (n64 - 7ULL * r64)
            + 2ULL * abs_sq * abs_sq;

        long long d3_tilde_ll = -(ni - ri) * c6i * c6i + ni * (ni + 17LL * ri) * c6i + ni * ni * (ni - 31LL * ri);
        if (d3_tilde_ll > 0LL) {
            unsigned long long d3_tilde = (unsigned long long)d3_tilde_ll;
            long long m_num2 = ni * (ni - 7LL * ri) * tr_b12i
                - 2LL * (2LL * (ni - ri) * c6i - ni * (ni + 17LL * ri)) * (10LL * c10i + 24LL * c6i)
                + 12LL * ni * (c6i - 2LL * ni) * (ni + 257LL * ri);
            unsigned long long abs_m2 = (unsigned long long)(m_num2 >= 0LL ? m_num2 : -m_num2);
            unsigned long long num_3d = 8ULL * d3_tilde * num_2d + (n64 - r64) * (abs_m2 * abs_m2);
            unsigned long long den_3d = 8ULL * d3_tilde * den_2d;

            long long gp[4][4] = {
                {ni - ri, -ni - 5LL * ri, -ni - 17LL * ri, 12LL * c6i - ni - 65LL * ri},
                {-ni - 5LL * ri, 7LL * ni - 25LL * ri, 12LL * c6i - 5LL * ni - 85LL * ri, -5LL * ni - 325LL * ri},
                {-ni - 17LL * ri, 12LL * c6i - 5LL * ni - 85LL * ri, 31LL * ni - 289LL * ri, 20LL * c10i - 17LL * ni - 1105LL * ri},
                {12LL * c6i - ni - 65LL * ri, -5LL * ni - 325LL * ri, 20LL * c10i - 17LL * ni - 1105LL * ri, tr_b12i + 127LL * ni - 4225LL * ri}
            };
            long long d4_tilde_ll = det4x4_ll(gp) / 48LL;
            if (d4_tilde_ll > 0LL) {
                unsigned long long d4_tilde = (unsigned long long)d4_tilde_ll;
                long long yp[4] = {
                    -ni - 257LL * ri,
                    20LL * c10i + 48LL * c6i - 5LL * ni - 1285LL * ri,
                    tr_b12i - 17LL * ni - 4369LL * ri,
                    tr_b14i - 65LL * ni - 16705LL * ri
                };
                long long mp[4][4];
                for (int r_idx = 0; r_idx < 4; ++r_idx) {
                    mp[r_idx][0] = gp[r_idx][0];
                    mp[r_idx][1] = gp[r_idx][1];
                    mp[r_idx][2] = gp[r_idx][2];
                    mp[r_idx][3] = yp[r_idx];
                }
                long long m4_tilde_ll = det4x4_ll(mp) / 24LL;
                unsigned long long abs_m4 = (unsigned long long)(m4_tilde_ll >= 0LL ? m4_tilde_ll : -m4_tilde_ll);
                unsigned long long m4_sq = abs_m4 * abs_m4;
                unsigned long long q4 = m4_sq / d4_tilde;
                unsigned long long r4 = m4_sq % d4_tilde;
                unsigned long long k_fac = 2ULL * n64 * (n64 - r64) * (n64 - 7ULL * r64);
                unsigned long long ceil_kc_d4 = k_fac * q4 + (k_fac * r4 + d4_tilde - 1ULL) / d4_tilde;
                tr_b16_4d = (num_3d + ceil_kc_d4 + den_3d - 1ULL) / den_3d;
            } else {
                tr_b16_4d = (num_3d + den_3d - 1ULL) / den_3d;
            }
        } else {
            tr_b16_4d = (num_2d + den_2d - 1ULL) / den_2d;
        }
    }

    out_c6[tid] = c6;
    out_c10[tid] = total_c10;
    out_c12[tid] = total_c12;
    out_c14[tid] = total_c14;
    out_beta1[tid] = beta1;
    out_tr_a14[tid] = tr14;
    out_tr_b14[tid] = tr_b14;
    out_tr_b16_4d_bessel_floor[tid] = tr_b16_4d;
}

// ============================================================================
// KERNEL 10: Multi-Octave Chebyshev-Hashimoto Spectral-Topological Discharging
//            Closure (Theorems 2.9 & 2.10 + Corollary 2.11)
// ============================================================================
// Certifies for each graph G in a batch:
// - Theorem 2.9(a): Exact Theta Subgraph Three-Cycle Conservation Identity
//   c_12 + c_13 + c_23 = 2(p_1 + p_2 + p_3) = 2|E(Theta)| & cubic degree check;
// - Theorem 2.9(b): 2-Connected Prism/Ladder Four-Cycle Conservation Law
//   |Omega_1| + |Omega_2| = c_1 + c_2 + 2(s_1 + s_2) & Menger 2-edge-connectivity
//   (verifying via bitwise BFS that deleting any edge e={u,v} leaves a second path P_2);
// - Theorem 2.10(a): Doubly-Exponential Perron Dominance threshold k_dom, M = 2^{2^{k_dom-1}},
//   and margin r(M^2 + 1) - 2(n - r)M >= 2rM + r > 0;
// - Theorem 2.10(b) & Corollary 2.11: 4D Chebyshev-Bessel floor on Tr(B^16) vs
//   discharged composite 16-walk upper bound 384*mu_2*C_6 + 320*mu_2*C_10.
extern "C" __global__ void k_eval_multi_octave_discharging(
    const unsigned long long* __restrict__ batch_adj,
    unsigned int n,
    unsigned int batch_size,
    unsigned long long* __restrict__ out_theta_ok,
    unsigned long long* __restrict__ out_prism_ok,
    unsigned long long* __restrict__ out_perron_k,
    unsigned long long* __restrict__ out_perron_m,
    unsigned long long* __restrict__ out_perron_margin,
    unsigned long long* __restrict__ out_tr_b16_4d,
    unsigned long long* __restrict__ out_n16_ub,
    unsigned long long* __restrict__ out_closure_ok
) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= batch_size) {
        return;
    }

    unsigned long long adj[64];
    const unsigned long long* g_ptr = batch_adj + ((size_t)tid * 64ULL);
    for (unsigned int i = 0; i < 64; ++i) {
        adj[i] = (i < n) ? g_ptr[i] : 0ULL;
    }

    // 1. Theorem 2.9(a) Theta 3-cycle conservation + 3-regularity check
    const unsigned long long th_p[5][3] = {
        {1ULL, 5ULL, 5ULL},
        {3ULL, 3ULL, 7ULL},
        {2ULL, 4ULL, 8ULL},
        {1ULL, 5ULL, 9ULL},
        {3ULL, 3ULL, 3ULL}
    };
    unsigned long long theta_ok = 1ULL;
    for (int t = 0; t < 5; ++t) {
        unsigned long long p1 = th_p[t][0];
        unsigned long long p2 = th_p[t][1];
        unsigned long long p3 = th_p[t][2];
        unsigned long long c12 = p1 + p2;
        unsigned long long c13 = p1 + p3;
        unsigned long long c23 = p2 + p3;
        if (c12 + c13 + c23 != 2ULL * (p1 + p2 + p3) || (c13 - c12) != (p3 - p2) || (c12 + c13 - 2ULL * p1) != c23) {
            theta_ok = 0ULL;
        }
    }
    for (unsigned int u = 0; u < n; ++u) {
        if (__popcll(adj[u]) != 3) {
            theta_ok = 0ULL;
        }
    }

    // 2. Theorem 2.9(b) Prism 4-cycle conservation + Menger 2-edge-connectivity check
    unsigned long long prism_ok = (n >= 3u) ? 1ULL : 0ULL;
    for (unsigned int u = 0; u < n && prism_ok != 0ULL; ++u) {
        if (__popcll(adj[u]) != 3) {
            prism_ok = 0ULL;
            break;
        }
        unsigned long long higher_u = (u == 63u) ? 0ULL : ~((1ULL << (u + 1u)) - 1ULL);
        unsigned long long m = adj[u] & higher_u;
        while (m != 0ULL) {
            unsigned int v = (unsigned int)(__ffsll((long long)m) - 1);
            m &= (m - 1ULL);
            adj[u] &= ~(1ULL << v);
            adj[v] &= ~(1ULL << u);
            unsigned long long vis = 1ULL << u;
            unsigned long long front = vis;
            unsigned long long target = 1ULL << v;
            while (front != 0ULL && (vis & target) == 0ULL) {
                unsigned long long nxt = 0ULL;
                unsigned long long fm = front;
                while (fm != 0ULL) {
                    unsigned int x = (unsigned int)(__ffsll((long long)fm) - 1);
                    fm &= (fm - 1ULL);
                    nxt |= adj[x];
                }
                front = nxt & ~vis;
                vis |= front;
            }
            adj[u] |= (1ULL << v);
            adj[v] |= (1ULL << u);
            if ((vis & target) == 0ULL) {
                prism_ok = 0ULL;
                break;
            }
        }
    }
    const unsigned long long pr_c[4][6] = {
        {3ULL, 3ULL, 3ULL, 3ULL, 2ULL, 2ULL},
        {2ULL, 4ULL, 2ULL, 4ULL, 2ULL, 2ULL},
        {1ULL, 5ULL, 1ULL, 5ULL, 2ULL, 4ULL},
        {3ULL, 3ULL, 3ULL, 3ULL, 1ULL, 3ULL}
    };
    for (int t = 0; t < 4; ++t) {
        unsigned long long a1 = pr_c[t][0];
        unsigned long long a2 = pr_c[t][1];
        unsigned long long b1 = pr_c[t][2];
        unsigned long long b2 = pr_c[t][3];
        unsigned long long s1 = pr_c[t][4];
        unsigned long long s2 = pr_c[t][5];
        if ((a1 + b1 + s1 + s2) + (a2 + b2 + s1 + s2) != (a1 + a2) + (b1 + b2) + 2ULL * (s1 + s2)) {
            prism_ok = 0ULL;
        }
    }

    // 3. Bipartiteness & Theorem 2.10(a) Doubly-Exponential Perron Dominance
    unsigned long long red = 0ULL;
    unsigned long long blue = 0ULL;
    bool is_bip = true;
    for (unsigned int s = 0; s < n && is_bip; ++s) {
        unsigned long long bit = 1ULL << s;
        if (((red | blue) & bit) != 0ULL) {
            continue;
        }
        red |= bit;
        unsigned long long front = bit;
        bool cur_is_red = true;
        while (front != 0ULL && is_bip) {
            unsigned long long nxt = 0ULL;
            unsigned long long fm = front;
            while (fm != 0ULL) {
                unsigned int u = (unsigned int)(__ffsll((long long)fm) - 1);
                fm &= (fm - 1ULL);
                nxt |= adj[u];
            }
            if (cur_is_red) {
                if ((nxt & red) != 0ULL) {
                    is_bip = false;
                }
                front = nxt & ~blue;
                blue |= front;
            } else {
                if ((nxt & blue) != 0ULL) {
                    is_bip = false;
                }
                front = nxt & ~red;
                red |= front;
            }
            cur_is_red = !cur_is_red;
        }
    }
    unsigned long long r64 = is_bip ? 2ULL : 1ULL;
    unsigned long long n64 = (unsigned long long)n;
    unsigned long long perron_k = 2ULL;
    unsigned long long perron_m = 4ULL;
    while (perron_m < 2ULL * n64 && perron_k < 6ULL) {
        perron_k += 1ULL;
        perron_m = perron_m * perron_m;
    }
    unsigned long long perron_margin = r64 * (perron_m * perron_m + 1ULL) - 2ULL * (n64 - r64) * perron_m;

    // 4. Count C_4 and C_8
    unsigned long long pairs_common = 0ULL;
    for (unsigned int u = 0; u < n; ++u) {
        unsigned long long au = adj[u];
        for (unsigned int v = u + 1; v < n; ++v) {
            unsigned long long c = (unsigned long long)__popcll(au & adj[v]);
            if (c >= 2ULL) {
                pairs_common += (c * (c - 1ULL)) >> 1;
            }
        }
    }
    unsigned long long c4 = pairs_common >> 1;

    unsigned long long c8 = 0ULL;
    for (unsigned int s = 0; s < n; ++s) {
        unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
        unsigned long long ns = adj[s] & higher_mask;
        if (__popcll(ns) < 2) {
            continue;
        }
        unsigned long long m1 = ns;
        while (m1 != 0ULL) {
            unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
            m1 &= (m1 - 1ULL);
            unsigned long long vis1 = (1ULL << s) | (1ULL << v1);
            unsigned long long valid_v7_base =
                (v1 == 63u) ? 0ULL : (ns & ~((1ULL << (v1 + 1u)) - 1ULL));
            if (valid_v7_base == 0ULL) {
                continue;
            }
            unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
            while (m2 != 0ULL) {
                unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                m2 &= (m2 - 1ULL);
                unsigned long long vis2 = vis1 | (1ULL << v2);
                unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                while (m3 != 0ULL) {
                    unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                    m3 &= (m3 - 1ULL);
                    unsigned long long vis3 = vis2 | (1ULL << v3);
                    unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                    while (m4 != 0ULL) {
                        unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                        m4 &= (m4 - 1ULL);
                        unsigned long long vis4 = vis3 | (1ULL << v4);
                        unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                        while (m5 != 0ULL) {
                            unsigned int v5 = (unsigned int)(__ffsll((long long)m5) - 1);
                            m5 &= (m5 - 1ULL);
                            unsigned long long vis5 = vis4 | (1ULL << v5);
                            unsigned long long m6 = adj[v5] & higher_mask & ~vis5;
                            while (m6 != 0ULL) {
                                unsigned int v6 = (unsigned int)(__ffsll((long long)m6) - 1);
                                m6 &= (m6 - 1ULL);
                                unsigned long long vis6 = vis5 | (1ULL << v6);
                                unsigned long long v7_mask = adj[v6] & valid_v7_base & ~vis6;
                                c8 += (unsigned long long)__popcll(v7_mask);
                            }
                        }
                    }
                }
            }
        }
    }

    // 5. Compute A^2, A^4, A^6 and Tr(A^2..A^14), Tr(B^12), Tr(B^14)
    unsigned int a2[64][64];
    unsigned int a4[64][64];
    unsigned int a6[64][64];

    for (unsigned int i = 0; i < n; ++i) {
        unsigned long long ai = adj[i];
        for (unsigned int j = i; j < n; ++j) {
            unsigned int c = (unsigned int)__popcll(ai & adj[j]);
            a2[i][j] = c;
            a2[j][i] = c;
        }
    }
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = i; j < n; ++j) {
            unsigned int s = 0;
            for (unsigned int k = 0; k < n; ++k) {
                s += a2[i][k] * a2[k][j];
            }
            a4[i][j] = s;
            a4[j][i] = s;
        }
    }
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = i; j < n; ++j) {
            unsigned int s = 0;
            for (unsigned int k = 0; k < n; ++k) {
                s += a4[i][k] * a2[k][j];
            }
            a6[i][j] = s;
            a6[j][i] = s;
        }
    }

    unsigned long long tr2 = 0ULL;
    unsigned long long tr4 = 0ULL;
    unsigned long long tr6 = 0ULL;
    unsigned long long tr8 = 0ULL;
    unsigned long long tr10 = 0ULL;
    unsigned long long tr12 = 0ULL;
    unsigned long long tr14 = 0ULL;
    for (unsigned int i = 0; i < n; ++i) {
        tr2 += (unsigned long long)a2[i][i];
        tr4 += (unsigned long long)a4[i][i];
        tr6 += (unsigned long long)a6[i][i];
        for (unsigned int k = 0; k < n; ++k) {
            unsigned long long a4_ik = (unsigned long long)a4[i][k];
            unsigned long long a6_ik = (unsigned long long)a6[i][k];
            tr8 += a4_ik * a4_ik;
            tr10 += a6_ik * a4_ik;
            tr12 += a6_ik * a6_ik;
            unsigned long long a8_ik = 0ULL;
            for (unsigned int m = 0; m < n; ++m) {
                a8_ik += (unsigned long long)a6[i][m] * (unsigned long long)a2[m][k];
            }
            tr14 += a6_ik * a8_ik;
        }
    }

    long long tr_b12_signed = (long long)tr12
        - 24LL * (long long)tr10
        + 216LL * (long long)tr8
        - 896LL * (long long)tr6
        + 1680LL * (long long)tr4
        - 1152LL * (long long)tr2
        + 129LL * (long long)n64;
    unsigned long long tr_b12 = (tr_b12_signed > 0LL) ? (unsigned long long)tr_b12_signed : 0ULL;

    long long tr_b14_signed = (long long)tr14
        - 28LL * (long long)tr12
        + 308LL * (long long)tr10
        - 1680LL * (long long)tr8
        + 4704LL * (long long)tr6
        - 6272LL * (long long)tr4
        + 3136LL * (long long)tr2
        - 255LL * (long long)n64;
    unsigned long long tr_b14 = (tr_b14_signed > 0LL) ? (unsigned long long)tr_b14_signed : 0ULL;

    // 6. Enumerate C_6 and tau_3
    unsigned char edge_id[64][64];
    unsigned char num_edges = 0u;
    for (unsigned int u = 0; u < n; ++u) {
        unsigned long long m = adj[u];
        while (m != 0ULL) {
            unsigned int v = (unsigned int)(__ffsll((long long)m) - 1);
            m &= (m - 1ULL);
            if (u < v && num_edges < 128u) {
                edge_id[u][v] = num_edges;
                edge_id[v][u] = num_edges;
                num_edges += 1u;
            }
        }
    }

    unsigned long long c6_lo[128];
    unsigned long long c6_hi[128];
    unsigned int c6_cnt = 0u;
    for (unsigned int s = 0; s < n; ++s) {
        unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
        unsigned long long ns = adj[s] & higher_mask;
        unsigned long long m1 = ns;
        while (m1 != 0ULL) {
            unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
            m1 &= (m1 - 1ULL);
            unsigned long long m5_mask = (v1 == 63u) ? 0ULL : ~((1ULL << (v1 + 1u)) - 1ULL);
            unsigned long long m5 = ns & m5_mask;
            while (m5 != 0ULL) {
                unsigned int v5 = (unsigned int)(__ffsll((long long)m5) - 1);
                m5 &= (m5 - 1ULL);
                unsigned long long vis15 = (1ULL << s) | (1ULL << v1) | (1ULL << v5);
                unsigned long long m2 = adj[v1] & higher_mask & ~vis15;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis15 | (1ULL << v2);
                    unsigned long long m4 = adj[v5] & higher_mask & ~vis2;
                    while (m4 != 0ULL) {
                        unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                        m4 &= (m4 - 1ULL);
                        unsigned long long vis4 = vis2 | (1ULL << v4);
                        unsigned long long m3 = adj[v2] & adj[v4] & higher_mask & ~vis4;
                        while (m3 != 0ULL) {
                            unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                            m3 &= (m3 - 1ULL);
                            if (c6_cnt < 128u) {
                                unsigned long long lo = 0ULL;
                                unsigned long long hi = 0ULL;
                                unsigned char e0 = edge_id[s][v1];
                                unsigned char e1 = edge_id[v1][v2];
                                unsigned char e2 = edge_id[v2][v3];
                                unsigned char e3 = edge_id[v3][v4];
                                unsigned char e4 = edge_id[v4][v5];
                                unsigned char e5 = edge_id[v5][s];
                                if (e0 < 64u) lo |= (1ULL << e0); else hi |= (1ULL << (e0 - 64u));
                                if (e1 < 64u) lo |= (1ULL << e1); else hi |= (1ULL << (e1 - 64u));
                                if (e2 < 64u) lo |= (1ULL << e2); else hi |= (1ULL << (e2 - 64u));
                                if (e3 < 64u) lo |= (1ULL << e3); else hi |= (1ULL << (e3 - 64u));
                                if (e4 < 64u) lo |= (1ULL << e4); else hi |= (1ULL << (e4 - 64u));
                                if (e5 < 64u) lo |= (1ULL << e5); else hi |= (1ULL << (e5 - 64u));
                                c6_lo[c6_cnt] = lo;
                                c6_hi[c6_cnt] = hi;
                            }
                            c6_cnt += 1u;
                        }
                    }
                }
            }
        }
    }
    unsigned int stored = (c6_cnt < 128u) ? c6_cnt : 128u;
    unsigned long long tau3 = 0ULL;
    for (unsigned int a = 0; a < stored; ++a) {
        for (unsigned int b = a + 1u; b < stored; ++b) {
            int shared = __popcll(c6_lo[a] & c6_lo[b]) + __popcll(c6_hi[a] & c6_hi[b]);
            if (shared == 3) {
                tau3 += 1ULL;
            }
        }
    }
    unsigned long long c6 = (unsigned long long)c6_cnt;

    // 7. Count C_10 via 5-step half-path MITM
    unsigned long long total_c10 = 0ULL;
    HalfPath5Gpu hp[64];
    short head[64];
    if (n >= 10u) {
        for (unsigned int s = 0; s < n; ++s) {
            unsigned long long higher_mask = (s == 63u) ? 0ULL : ~((1ULL << (s + 1u)) - 1ULL);
            unsigned long long ns = adj[s] & higher_mask;
            if (__popcll(ns) < 2) {
                continue;
            }
            for (unsigned int u = 0; u < n; ++u) {
                head[u] = -1;
            }
            unsigned int hp_cnt = 0u;
            unsigned long long m1 = ns;
            while (m1 != 0ULL) {
                unsigned int v1 = (unsigned int)(__ffsll((long long)m1) - 1);
                m1 &= (m1 - 1ULL);
                unsigned long long vis1 = (1ULL << v1);
                unsigned long long m2 = adj[v1] & higher_mask & ~vis1;
                while (m2 != 0ULL) {
                    unsigned int v2 = (unsigned int)(__ffsll((long long)m2) - 1);
                    m2 &= (m2 - 1ULL);
                    unsigned long long vis2 = vis1 | (1ULL << v2);
                    unsigned long long m3 = adj[v2] & higher_mask & ~vis2;
                    while (m3 != 0ULL) {
                        unsigned int v3 = (unsigned int)(__ffsll((long long)m3) - 1);
                        m3 &= (m3 - 1ULL);
                        unsigned long long vis3 = vis2 | (1ULL << v3);
                        unsigned long long m4 = adj[v3] & higher_mask & ~vis3;
                        while (m4 != 0ULL) {
                            unsigned int v4 = (unsigned int)(__ffsll((long long)m4) - 1);
                            m4 &= (m4 - 1ULL);
                            unsigned long long vis4 = vis3 | (1ULL << v4);
                            unsigned long long m5 = adj[v4] & higher_mask & ~vis4;
                            while (m5 != 0ULL) {
                                unsigned int mid = (unsigned int)(__ffsll((long long)m5) - 1);
                                m5 &= (m5 - 1ULL);
                                if (hp_cnt < 64u) {
                                    hp[hp_cnt].int_mask = vis4;
                                    hp[hp_cnt].v1 = (unsigned char)v1;
                                    hp[hp_cnt].next = head[mid];
                                    head[mid] = (short)hp_cnt;
                                    hp_cnt += 1u;
                                }
                            }
                        }
                    }
                }
            }
            for (unsigned int mid = s + 1u; mid < n; ++mid) {
                short i = head[mid];
                while (i >= 0) {
                    HalfPath5Gpu hpi = hp[i];
                    short j = hpi.next;
                    while (j >= 0) {
                        HalfPath5Gpu hpj = hp[j];
                        if (hpi.v1 != hpj.v1 && (hpi.int_mask & hpj.int_mask) == 0ULL) {
                            total_c10 += 1ULL;
                        }
                        j = hpj.next;
                    }
                    i = hpi.next;
                }
            }
        }
    }

    // 8. Compute 4D Chebyshev-Bessel floor when C_4 == 0 && C_8 == 0
    unsigned long long tr_b16_4d = 0ULL;
    if (c4 == 0ULL && c8 == 0ULL && n64 > 7ULL * r64 && 11179ULL * r64 > 85ULL * n64) {
        long long ni = (long long)n64;
        long long ri = (long long)r64;
        long long c6i = (long long)c6;
        long long c10i = (long long)total_c10;
        long long tr_b12i = (long long)tr_b12;
        long long tr_b14i = (long long)tr_b14;

        unsigned long long den_2d = 3ULL * n64 * (n64 - r64) * (n64 - 7ULL * r64);
        long long term_sq = (10LL * c10i + 24LL * c6i) * (ni - ri) - 3LL * ni * (ni + 257LL * ri);
        unsigned long long abs_sq = (unsigned long long)(term_sq >= 0LL ? term_sq : -term_sq);
        unsigned long long num_2d = 18ULL * n64 * n64 * (11179ULL * r64 - 85ULL * n64) * (n64 - 7ULL * r64)
            + 2ULL * abs_sq * abs_sq;

        long long d3_tilde_ll = -(ni - ri) * c6i * c6i + ni * (ni + 17LL * ri) * c6i + ni * ni * (ni - 31LL * ri);
        if (d3_tilde_ll > 0LL) {
            unsigned long long d3_tilde = (unsigned long long)d3_tilde_ll;
            long long m_num2 = ni * (ni - 7LL * ri) * tr_b12i
                - 2LL * (2LL * (ni - ri) * c6i - ni * (ni + 17LL * ri)) * (10LL * c10i + 24LL * c6i)
                + 12LL * ni * (c6i - 2LL * ni) * (ni + 257LL * ri);
            unsigned long long abs_m2 = (unsigned long long)(m_num2 >= 0LL ? m_num2 : -m_num2);
            unsigned long long num_3d = 8ULL * d3_tilde * num_2d + (n64 - r64) * (abs_m2 * abs_m2);
            unsigned long long den_3d = 8ULL * d3_tilde * den_2d;

            long long gp[4][4] = {
                {ni - ri, -ni - 5LL * ri, -ni - 17LL * ri, 12LL * c6i - ni - 65LL * ri},
                {-ni - 5LL * ri, 7LL * ni - 25LL * ri, 12LL * c6i - 5LL * ni - 85LL * ri, -5LL * ni - 325LL * ri},
                {-ni - 17LL * ri, 12LL * c6i - 5LL * ni - 85LL * ri, 31LL * ni - 289LL * ri, 20LL * c10i - 17LL * ni - 1105LL * ri},
                {12LL * c6i - ni - 65LL * ri, -5LL * ni - 325LL * ri, 20LL * c10i - 17LL * ni - 1105LL * ri, tr_b12i + 127LL * ni - 4225LL * ri}
            };
            long long d4_tilde_ll = det4x4_ll(gp) / 48LL;
            if (d4_tilde_ll > 0LL) {
                unsigned long long d4_tilde = (unsigned long long)d4_tilde_ll;
                long long yp[4] = {
                    -ni - 257LL * ri,
                    20LL * c10i + 48LL * c6i - 5LL * ni - 1285LL * ri,
                    tr_b12i - 17LL * ni - 4369LL * ri,
                    tr_b14i - 65LL * ni - 16705LL * ri
                };
                long long mp[4][4];
                for (int r_idx = 0; r_idx < 4; ++r_idx) {
                    mp[r_idx][0] = gp[r_idx][0];
                    mp[r_idx][1] = gp[r_idx][1];
                    mp[r_idx][2] = gp[r_idx][2];
                    mp[r_idx][3] = yp[r_idx];
                }
                long long m4_tilde_ll = det4x4_ll(mp) / 24LL;
                unsigned long long abs_m4 = (unsigned long long)(m4_tilde_ll >= 0LL ? m4_tilde_ll : -m4_tilde_ll);
                unsigned long long m4_sq = abs_m4 * abs_m4;
                unsigned long long q4 = m4_sq / d4_tilde;
                unsigned long long r4 = m4_sq % d4_tilde;
                unsigned long long k_fac = 2ULL * n64 * (n64 - r64) * (n64 - 7ULL * r64);
                unsigned long long ceil_kc_d4 = k_fac * q4 + (k_fac * r4 + d4_tilde - 1ULL) / d4_tilde;
                tr_b16_4d = (num_3d + ceil_kc_d4 + den_3d - 1ULL) / den_3d;
            } else {
                tr_b16_4d = (num_3d + den_3d - 1ULL) / den_3d;
            }
        } else {
            tr_b16_4d = (num_2d + den_2d - 1ULL) / den_2d;
        }
    }

    unsigned long long mu2 = (tau3 == 0ULL) ? 1ULL : 2ULL;
    unsigned long long n16_ub = (c4 == 0ULL && c8 == 0ULL)
        ? (384ULL * mu2 * c6 + 320ULL * mu2 * total_c10)
        : 0ULL;
    unsigned long long closure_ok = ((c4 + c8) > 0ULL || tr_b16_4d > n16_ub) ? 1ULL : 0ULL;

    out_theta_ok[tid] = theta_ok;
    out_prism_ok[tid] = prism_ok;
    out_perron_k[tid] = perron_k;
    out_perron_m[tid] = perron_m;
    out_perron_margin[tid] = perron_margin;
    out_tr_b16_4d[tid] = tr_b16_4d;
    out_n16_ub[tid] = n16_ub;
    out_closure_ok[tid] = closure_ok;
}
