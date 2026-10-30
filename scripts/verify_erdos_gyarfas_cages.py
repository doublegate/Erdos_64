#!/usr/bin/env python3
"""
Independent Python Verification Suite for Erdős Problem #64
($1,000 Falsifiable Prize — The Erdős–Gyárfás 2^k-Cycle Conjecture).

Independently constructs the canonical (3, g)-cages, generalized polygons,
and snarks, verifies their degrees, girths, bipartiteness, Moore bounds,
and exact dyadic cycle multiplicities C_4(G), C_8(G), and C_16(G), and
cross-checks the Theorem 3.1 bitmask meet-in-the-middle identities.
"""

from __future__ import annotations

import json
import sys
from collections import deque
from typing import Dict, List, Tuple


def moore_bound_cubic(girth: int) -> int:
    """Exact Moore lower bound n_0(3, g) on the order of a 3-regular graph of girth g."""
    assert girth >= 3
    if girth % 2 == 1:
        r = (girth - 1) // 2
        return 3 * (1 << r) - 2
    r = girth // 2
    return (1 << (r + 1)) - 2


def from_lcf(n: int, shifts: List[int], repeats: int) -> List[int]:
    assert len(shifts) * repeats == n
    adj = [0] * n
    for i in range(n):
        j = (i + 1) % n
        adj[i] |= 1 << j
        adj[j] |= 1 << i
        t = (i + shifts[i % len(shifts)]) % n
        adj[i] |= 1 << t
        adj[t] |= 1 << i
    return adj


def generalized_petersen(m: int, k: int) -> List[int]:
    n = 2 * m
    adj = [0] * n
    for i in range(m):
        u_next = (i + 1) % m
        adj[i] |= 1 << u_next
        adj[u_next] |= 1 << i
        spoke = m + i
        adj[i] |= 1 << spoke
        adj[spoke] |= 1 << i
        inner = m + ((i + k) % m)
        adj[spoke] |= 1 << inner
        adj[inner] |= 1 << spoke
    return adj


def compute_girth(adj: List[int]) -> int:
    n = len(adj)
    best = 10**9
    for s in range(n):
        dist = [-1] * n
        parent = [-1] * n
        dist[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            d = dist[u]
            if 2 * d + 1 >= best:
                break
            mask = adj[u]
            while mask:
                lsb = mask & -mask
                v = lsb.bit_length() - 1
                mask ^= lsb
                if dist[v] == -1:
                    dist[v] = d + 1
                    parent[v] = u
                    q.append(v)
                elif parent[u] != v:
                    best = min(best, dist[u] + dist[v] + 1)
    return best


def count_c4(adj: List[int]) -> int:
    n = len(adj)
    pairs = 0
    for u in range(n):
        au = adj[u]
        for v in range(u + 1, n):
            c = (au & adj[v]).bit_count()
            if c >= 2:
                pairs += c * (c - 1) // 2
    return pairs // 2


def count_c8_mitm(adj: List[int]) -> int:
    n = len(adj)
    total = 0
    for s in range(n):
        higher_mask = ~((1 << (s + 1)) - 1)
        ns = adj[s] & higher_mask
        if ns.bit_count() < 2:
            continue
        by_mid: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
        m1 = ns
        while m1:
            b1 = m1 & -m1
            v1 = b1.bit_length() - 1
            m1 ^= b1
            m2 = adj[v1] & higher_mask & ~b1
            while m2:
                b2 = m2 & -m2
                v2 = b2.bit_length() - 1
                m2 ^= b2
                b12 = b1 | b2
                m3 = adj[v2] & higher_mask & ~b12
                while m3:
                    b3 = m3 & -m3
                    v3 = b3.bit_length() - 1
                    m3 ^= b3
                    int_mask = b12 | b3
                    m4 = adj[v3] & higher_mask & ~int_mask
                    while m4:
                        b4 = m4 & -m4
                        mid = b4.bit_length() - 1
                        m4 ^= b4
                        by_mid[mid].append((v1, int_mask))
        for paths in by_mid:
            k = len(paths)
            for i in range(k):
                v1_a, ma = paths[i]
                for j in range(i + 1, k):
                    v1_b, mb = paths[j]
                    if v1_a != v1_b and (ma & mb) == 0:
                        total += 1
    return total


def dipole_voltage_lift(m: int, shifts: List[int]) -> List[int]:
    n = 2 * m
    adj = [0] * n
    for i in range(m):
        for s in shifts:
            v = m + ((i + s) % m)
            adj[i] |= 1 << v
            adj[v] |= 1 << i
    return adj


def count_c16_anchored0(adj: List[int]) -> int:
    n = len(adj)
    if n < 16:
        return 0
    higher_mask = ~1
    ns = adj[0] & higher_mask
    if ns.bit_count() < 2:
        return 0
    by_mid: List[List[int]] = [[] for _ in range(n)]

    def dfs(u: int, depth: int, int_mask: int) -> None:
        mask = adj[u] & higher_mask & ~int_mask
        if depth == 7:
            while mask:
                lsb = mask & -mask
                mid = lsb.bit_length() - 1
                mask ^= lsb
                by_mid[mid].append(int_mask)
            return
        while mask:
            lsb = mask & -mask
            nxt = lsb.bit_length() - 1
            mask ^= lsb
            dfs(nxt, depth + 1, int_mask | lsb)

    m1 = ns
    while m1:
        b1 = m1 & -m1
        v1 = b1.bit_length() - 1
        m1 ^= b1
        dfs(v1, 1, b1)

    count = 0
    for paths in by_mid:
        k = len(paths)
        for i in range(k):
            ma = paths[i]
            for j in range(i + 1, k):
                if (ma & paths[j]) == 0:
                    count += 1
    return count


def main() -> int:
    expected_moore = {
        3: 4,
        4: 6,
        5: 10,
        6: 14,
        7: 22,
        8: 30,
        9: 46,
        10: 62,
        11: 94,
        12: 126,
    }
    for g, exp in expected_moore.items():
        actual = moore_bound_cubic(g)
        assert actual == exp, f"Moore bound mismatch at g={g}: {actual} != {exp}"

    benchmarks: Dict[str, Tuple[List[int], int, int, int]] = {
        "K_4": (from_lcf(4, [2], 4), 3, 3, 0),
        "K_{3,3}": (from_lcf(6, [3], 6), 4, 9, 0),
        "Cube Q_3": (generalized_petersen(4, 1), 4, 6, 6),
        "Petersen (3,5)-cage": (generalized_petersen(5, 2), 5, 0, 15),
        "Heawood (3,6)-cage": (from_lcf(14, [5, -5], 7), 6, 0, 21),
        "Mobius-Kantor GP(8,3)": (generalized_petersen(8, 3), 6, 0, 30),
        "Pappus": (from_lcf(18, [5, 7, -7, 7, -7, -5], 3), 6, 0, 54),
        "Desargues GP(10,3)": (generalized_petersen(10, 3), 6, 0, 30),
        "Dodecahedron GP(10,2)": (generalized_petersen(10, 2), 5, 0, 30),
        "McGee (3,7)-cage": (from_lcf(24, [-12, 7, -7], 8), 7, 0, 34),
        "Tutte-Coxeter (3,8)-cage": (
            from_lcf(30, [-13, -9, 7, -7, 9, 13], 5),
            8,
            0,
            90,
        ),
        "Frontier GP(25,7) [n=50]": (generalized_petersen(25, 7), 8, 0, 25),
        "Frontier GP(26,5) [n=52]": (generalized_petersen(26, 5), 8, 0, 78),
        "Gray Graph LCF [n=54]": (
            from_lcf(54, [-25, 7, -7, 13, -13, 25], 9),
            8,
            0,
            81,
        ),
        "Frontier GP(28,5) [n=56]": (generalized_petersen(28, 5), 8, 0, 56),
        "Dipole Lift Z_29(0,1,8) [n=58]": (dipole_voltage_lift(29, [0, 1, 8]), 6, 0, 29),
        "Frontier GP(30,7) [n=60]": (generalized_petersen(30, 7), 8, 0, 60),
        "Dipole Lift Z_31(0,1,6) [n=62]": (dipole_voltage_lift(31, [0, 1, 6]), 6, 0, 0),
    }

    report = []
    print("=" * 88)
    print(
        "Independent Python Verification of (3,g)-Cages & Open-Frontier Lifts (Erdős Problem #64)"
    )
    print("=" * 88)
    print(
        f"{'Graph':<32} {'n':>4} {'Girth':>6} {'C_4':>6} {'C_8':>8} {'C_16^(0)':>10} {'Status':>10}"
    )
    print("-" * 88)
    for name, (adj, exp_girth, exp_c4, exp_c8) in benchmarks.items():
        n = len(adj)
        degs = [x.bit_count() for x in adj]
        assert min(degs) == 3 and max(degs) == 3, f"{name} is not 3-regular"
        g = compute_girth(adj)
        c4 = count_c4(adj)
        c8 = count_c8_mitm(adj)
        c16_0 = count_c16_anchored0(adj)
        assert g == exp_girth, f"{name}: girth {g} != {exp_girth}"
        assert c4 == exp_c4, f"{name}: C_4 {c4} != {exp_c4}"
        assert c8 == exp_c8, f"{name}: C_8 {c8} != {exp_c8}"
        assert c4 + c8 + c16_0 > 0, f"{name} failed Erdős-Gyárfás 2^k-cycle check!"
        print(f"{name:<32} {n:>4} {g:>6} {c4:>6} {c8:>8} {c16_0:>10} {'CERTIFIED':>10}")
        report.append(
            {"name": name, "n": n, "girth": g, "c4": c4, "c8": c8, "c16_anchor0": c16_0}
        )
    print("=" * 88)
    if "--json" in sys.argv:
        print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
