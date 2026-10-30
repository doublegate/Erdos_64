#!/usr/bin/env python3
"""
Reproducible Publication Figure Generator for Erdős Problem #64 (Phases 1–3).
Generates PDF and 300-DPI PNG figures in `paper/figures/` following Okabe-Ito
colorblind-safe scientific visualization standards:
  1. fig1_cubic_cage_moore_bounds.{pdf,png}
  2. fig2_cage_dyadic_multiplicities.{pdf,png}
  3. fig3_2opt_simulated_annealing_landscape.{pdf,png}
  4. fig4_petersen_heawood_cycle_spectra.{pdf,png}
  5. fig5_phase3_commutator_and_tabu_reduction.{pdf,png}
"""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parent.parent / "paper" / "figures"

# Okabe-Ito colorblind-safe palette
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
RED = "#D55E00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
GRAY = "#555555"


def save_both(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{stem}.pdf", dpi=300)
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=300)
    plt.close(fig)


def plot_fig1_moore_bounds() -> None:
    girths = np.arange(3, 13)
    moore = np.array([4, 6, 10, 14, 22, 30, 46, 62, 94, 126])
    cage = np.array([4, 6, 10, 14, 24, 30, 58, 70, 112, 126])
    min_dyadic = [4, 4, 8, 8, 8, 8, 16, 16, 16, 16]

    fig, ax = plt.subplots(figsize=(8.2, 4.6), layout="constrained")
    ax.plot(girths, moore, marker="o", lw=2.2, color=BLUE, label=r"Cubic Moore Bound $n_0(3, g)$")
    ax.plot(girths, cage, marker="s", lw=2.2, ls="--", color=RED, label=r"Exact $(3, g)$-Cage Order $f(3, g)$")
    for g, c_ord, d_len in zip(girths, cage, min_dyadic):
        ax.annotate(
            f"min $2^k={d_len}$",
            (g, c_ord),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            fontsize=8,
            color="#222222",
        )
    ax.axhspan(50, 64, color=GREEN, alpha=0.14, label=r"Open 64-Bit Cubic Frontier ($50 \leq n \leq 64$)")
    ax.set_yscale("log", base=2)
    ax.set_ylim(3.2, 185)
    ax.set_xticks(girths)
    ax.set_xlabel(r"Girth $g(G)$", fontsize=10.5)
    ax.set_ylabel(r"Vertex Count $n$ ($\log_2$ scale)", fontsize=10.5)
    ax.set_title(r"Proposition 2.1 & Theorem 2.1: Cubic Moore Bounds $n_0(3,g)$ vs. Exact $(3,g)$-Cage Orders", fontsize=11)
    ax.grid(True, which="both", ls=":", alpha=0.45)
    ax.legend(loc="upper left", fontsize=9)
    save_both(fig, "fig1_cubic_cage_moore_bounds")


def plot_fig2_cage_multiplicities() -> None:
    names = [
        "K4 (n=4)", "K3,3 (n=6)", "Q3 (n=8)", "M8 (n=8)", "Petersen (10)",
        "Tietze (12)", "Heawood (14)", "GP(8,3) (16)", "Pappus (18)", "Blanusa B1 (18)",
        "Desargues (20)", "Dodec (20)", "Flower J5 (20)", "GP(11,2) (22)",
        "McGee (24)", "Nauru (24)", "Coxeter (28)", "Tutte-Cox (30)"
    ]
    c4 = np.array([3, 9, 6, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    c8 = np.array([0, 0, 6, 5, 15, 9, 21, 30, 54, 12, 30, 30, 15, 33, 34, 54, 21, 90])
    c16 = np.array([0, 0, 0, 0, 0, 0, 0, 6, 216, 101, 360, 90, 180, 220, 492, 1062, 609, 1620])

    x = np.arange(len(names))
    w = 0.27
    fig, ax = plt.subplots(figsize=(10.5, 4.8), layout="constrained")
    ax.bar(x - w, c4 + 0.5, width=w, color=BLUE, label=r"$C_4(G)$ ($+0.5$ offset)")
    ax.bar(x, c8 + 0.5, width=w, color=ORANGE, label=r"$C_8(G)$ ($+0.5$ offset)")
    ax.bar(x + w, c16 + 0.5, width=w, color=RED, label=r"$C_{16}(G)$ ($+0.5$ offset)")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=38, ha="right", fontsize=8.5)
    ax.set_ylabel(r"Undirected Cycle Multiplicity ($\log_{10}$ scale)", fontsize=10)
    ax.set_title(r"Theorem 2.1: Certified Dyadic Cycle Multiplicities Across 18 Canonical $(3,g)$-Cages & Snarks", fontsize=11)
    ax.grid(True, axis="y", ls=":", alpha=0.45)
    ax.legend(loc="upper left", fontsize=9)
    save_both(fig, "fig2_cage_dyadic_multiplicities")


def plot_fig3_2opt_seesaw() -> None:
    ns_small = [10, 12, 14, 16, 18, 20, 22, 24]
    c8_small = [15, 13, 9, 9, 8, 6, 3, 4]
    c16_small = [0, 0, 0, 8, 100, 144, 303, 438]

    ns_front = [50, 52, 54, 56, 58, 60, 62]
    c8_front = [0, 0, 1, 0, 1, 0, 0]
    c16_front = [1317, 1384, 1323, 1400, 1427, 1416, 1519]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.4), layout="constrained")
    ax1.plot(ns_small, c8_small, marker="o", lw=2.2, color=BLUE, label=r"Min $C_8(G)$ ($g \geq 5$)")
    ax1_r = ax1.twinx()
    ax1_r.plot(ns_small, c16_small, marker="s", lw=2.2, ls="--", color=RED, label=r"Induced $C_{16}(G)$")
    ax1.set_xlabel(r"Order $n$ ($10 \leq n \leq 24$)", fontsize=10)
    ax1.set_ylabel(r"Minimum $C_8(G)$", color=BLUE, fontsize=10)
    ax1_r.set_ylabel(r"Induced $C_{16}(G)$", color=RED, fontsize=10)
    ax1.set_title(r"(a) Classical Regime ($10 \leq n \leq 24$): $C_8 \leftrightarrow C_{16}$ Seesaw", fontsize=10.5)
    ax1.grid(True, ls=":", alpha=0.45)

    ax2.plot(ns_front, c8_front, marker="o", lw=2.2, color=BLUE, label=r"Min $C_8(G)$ ($0$ on 5 orders!)")
    ax2.plot(ns_front, c16_front, marker="s", lw=2.2, ls="--", color=RED, label=r"Induced $C_{16}(G)$ ($1317$–$1519$)")
    ax2.set_xlabel(r"Frontier Order $n$ ($50 \leq n \leq 62$)", fontsize=10)
    ax2.set_ylabel(r"Cycle Multiplicity", fontsize=10)
    ax2.set_title(r"(b) Open Frontier ($50 \leq n \leq 62$): $C_4=C_8=0$ vs. $C_{16}$ Surge", fontsize=10.5)
    ax2.grid(True, ls=":", alpha=0.45)
    ax2.legend(loc="center right", fontsize=9)
    save_both(fig, "fig3_2opt_simulated_annealing_landscape")


def plot_fig4_cycle_spectra() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.2), layout="constrained")
    specs = [
        ("Petersen (3,5)-Cage (n=10, g=5)", [5, 6, 8, 9], [12, 10, 15, 20]),
        ("Heawood (3,6)-Cage (n=14, g=6)", [6, 8, 10, 12, 14], [28, 21, 84, 168, 24]),
        ("Blanusa Snark B1 (n=18, g=5)", list(range(5, 19)), [8, 14, 28, 12, 46, 127, 148, 207, 292, 239, 178, 101, 44, 0]),
        ("Isaacs Flower Snark J5 (n=20, g=5)", list(range(5, 21)), [5, 15, 20, 15, 35, 102, 165, 265, 365, 470, 519, 180, 295, 195, 60, 0]),
    ]
    for ax, (title, lengths, counts) in zip(axes.flat, specs):
        colors = [RED if (l >= 4 and (l & (l - 1)) == 0) else BLUE for l in lengths]
        ax.bar(lengths, counts, color=colors, edgecolor="#222222", lw=0.6)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r"Cycle Length $\ell$", fontsize=9.5)
        ax.set_ylabel(r"$C_\ell(G)$", fontsize=9.5)
        ax.set_xticks(lengths)
        ax.grid(True, axis="y", ls=":", alpha=0.45)
    save_both(fig, "fig4_petersen_heawood_cycle_spectra")


def plot_fig5_phase3_commutator_and_tabu() -> None:
    ns = np.array([50, 52, 54, 56, 58, 60, 62, 64])
    m_vals = ns // 2
    univ_floor = 18 * m_vals
    abelian_min = np.array([2550, 2483, 2484, 2331, 2552, 2550, 2294, 2232])
    nonabelian_lcf = np.array([2550, 2483, 2385, 2331, 2552, 2340, 2294, 2232])
    tabu_g6 = np.array([1417, 1463, 1418, 1569, 1307, 1342, 1202, 1271])
    tabu_g5 = np.array([1143, 1204, 1036, 1099, 1098, 1149, 1055, 1071])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.7), layout="constrained")

    # Panel A: Commutator Hierarchy & Tabu Reduction across 50 <= n <= 64
    ax1.plot(ns, abelian_min, marker="o", lw=2.2, color=RED, label=r"Global Abelian Dipole Min $C_{16}$ ($C_4=C_8=0$)")
    ax1.plot(ns, nonabelian_lcf, marker="D", lw=2.0, ls="--", color=ORANGE, label=r"Non-Abelian / Period-4 LCF ($n=54: 2385,\, n=60: 2340$)")
    ax1.plot(ns, tabu_g6, marker="s", lw=2.2, color=BLUE, label=r"Tabu Search ($g \geq 6,\, C_4=C_8=0$: min $1202$ at $n=62$)")
    ax1.plot(ns, tabu_g5, marker="^", lw=2.2, color=GREEN, label=r"Tabu Search ($g \geq 5,\, C_4=C_8=0$: min $1036$ at $n=54$)")
    ax1.plot(ns, univ_floor, marker="v", lw=1.8, ls=":", color=PURPLE, label=r"Free Abelian $\mathbb{Z}^2$ Commutator Floor $18 |\Gamma| = 9n$")
    ax1.scatter([62], [2883], color="#111111", marker="*", s=110, zorder=5, label=r"Phase 2 $\mathbb{Z}_{31}(0,1,6)$ ($C_{16}=2883=93\times 31$)")

    ax1.set_xticks(ns)
    ax1.set_xlabel(r"Frontier Order $n = 2|\Gamma|$ ($50 \leq n \leq 64$)", fontsize=10)
    ax1.set_ylabel(r"16-Cycle Multiplicity $C_{16}(G)$ (with $C_4=C_8=0$)", fontsize=10)
    ax1.set_title(r"(a) Theorems 2.3 & 3.2(b): Breaking the Abelian $C_{16}$ Floor", fontsize=10.5)
    ax1.grid(True, ls=":", alpha=0.45)
    ax1.legend(loc="upper right", fontsize=7.8)

    # Panel B: Theorem 3.2(a) Hashimoto Spectral Trace Decomposition Tr(B^16)/32 = C_16 + N_16/32
    labels = [
        "Z25(1,10)\nn=50, g=6",
        "Per3(Z27)\nn=54, g=6",
        "LCF4_60\nn=60, g=6",
        "Z31(1,6)\nn=62, g=6",
        "Z31(1,12)\nn=62, g=6",
        "Z32(1,7)\nn=64, g=6",
        "Z8x3Z4\nn=64, g=8",
    ]
    c16_bars = np.array([2550, 2385, 2340, 2883, 2294, 2232, 3652])
    db_bars = np.array([300, 396, 360, 372, 372, 384, 0])
    th_bars = np.array([1300, 1215, 1500, 930, 1550, 1600, 40])
    spectral_floor_orbits = np.array([3537, 3460, 3348, 3311, 3311, 3275])
    x = np.arange(len(labels))

    ax2.bar(x, c16_bars, width=0.56, color=BLUE, label=r"Simple 16-Cycles $C_{16}(G)$")
    ax2.bar(
        x,
        db_bars,
        width=0.56,
        bottom=c16_bars,
        color=GREEN,
        edgecolor="#222222",
        lw=0.5,
        label=r"Dumbbell Walks $\mathcal{N}_{16}^{\mathrm{dumbbell}}(G)/32 = 12\,C_6$",
    )
    ax2.bar(
        x,
        th_bars,
        width=0.56,
        bottom=c16_bars + db_bars,
        color=ORANGE,
        hatch="//",
        edgecolor="#222222",
        lw=0.5,
        label=r"Bridgeless Theta Walks $\mathcal{N}_{16}^{\mathrm{theta}}(G)/32$",
    )
    ax2.plot(
        x[:6],
        spectral_floor_orbits,
        marker="o",
        lw=2.0,
        ls="--",
        color=RED,
        zorder=6,
        label=r"Thm 2.4(2) Spectral Floor $\lceil \mathrm{Tr}(B^{16})_{\min}/32 \rceil$",
    )
    ax2.set_ylim(0, 5250)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=7.6)
    ax2.set_ylabel(r"Hashimoto Normalized Trace $\mathrm{Tr}(B^{16}) / 32$", fontsize=10)
    ax2.set_title(r"(b) Theorems 2.4–2.5 & 3.2(a): $\mathrm{Tr}(B^{16})/32$ Topological Split", fontsize=10.5)
    ax2.grid(True, axis="y", ls=":", alpha=0.45)
    ax2.legend(loc="upper left", fontsize=7.3)

    save_both(fig, "fig5_phase3_commutator_and_tabu_reduction")


def plot_fig6_chebyshev_bessel_cascade() -> None:
    labels = [
        "Z25(1,10)\nn=50",
        "LCF4_52\nn=52",
        "Z27(1,6)\nn=54",
        "Per3_54\nn=54",
        "LCF4_56\nn=56",
        "Z29(1,9)\nn=58",
        "LCF4_60\nn=60",
        "Z31(1,6)\nn=62",
        "Z31(1,12)\nn=62",
        "Z32(1,7)\nn=64",
    ]
    lb_1d = np.array([1937, 1919, 1859, 639, 1799, 1739, 1680, 2010, 1450, 1485])
    lb_2d = np.array([1939, 1922, 1860, 640, 1802, 1743, 1686, 2015, 1465, 1506])
    lb_3d = np.array([2184, 2124, 2088, 864, 1953, 1868, 1759, 2649, 1581, 1523])
    lb_4d = np.array([2419, 2278, 2213, 931, 2095, 2059, 1896, 2867, 1778, 1666])
    actual_c16 = np.array([2550, 2483, 2484, 2385, 2331, 2552, 2340, 2883, 2294, 2232])

    x = np.arange(len(labels))
    w = 0.16
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.6), layout="constrained")

    ax1.bar(x - 2 * w, lb_1d, width=w, color=SKY, label=r"1D Cauchy–Schwarz")
    ax1.bar(x - w, lb_2d, width=w, color=BLUE, label=r"2D Bessel ($\mathrm{span}\{1,\lambda^2\}$)")
    ax1.bar(x, lb_3d, width=w, color=GREEN, label=r"3D Bessel ($\mathrm{span}\{1,\lambda^2,\lambda^4\}$)")
    ax1.bar(x + w, lb_4d, width=w, color=ORANGE, label=r"4D Bessel ($\mathrm{span}\{1,\dots,\lambda^6\}$)")
    ax1.bar(x + 2 * w, actual_c16, width=w, color=RED, label=r"Actual $C_{16}(G)$")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=6.4)
    ax1.set_ylabel(r"Certified 16-Cycles $C_{16}(G)$", fontsize=10)
    ax1.set_title(r"(a) Theorems 2.4–2.8: 1D $\to$ 4D Chebyshev–Bessel $C_{16}$ Cascade", fontsize=10.5)
    ax1.grid(True, axis="y", ls=":", alpha=0.45)
    ax1.legend(loc="upper right", fontsize=7.6)

    regimes = [
        "Girth $g \\geq 9$\n($r=1$, 1D)",
        "Girth $g \\geq 10$\n($r=2$, 2D)",
        "Girth $g \\geq 10$\n($r=2$, 3D Env)",
        "Girth $g \\geq 12$\n($r=2$, 3D)",
        "Girth $g \\geq 16$\n($r=2$, 4D)",
    ]
    n_max = np.array([131, 90, 250, 290, 1020])
    y_pos = np.arange(len(regimes))
    bars = ax2.barh(y_pos, n_max, color=[SKY, BLUE, GREEN, ORANGE, RED], edgecolor="#222222", height=0.55)
    for bar, val in zip(bars, n_max):
        ax2.text(val + 18, bar.get_y() + bar.get_height() / 2, f"$n \\leq {val}$", va="center", fontsize=9)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(regimes, fontsize=8.8)
    ax2.set_xlim(0, 1180)
    ax2.set_xlabel(r"Certified Maximum Order $n$ with Unconditional $C_{16}(G) > 0$", fontsize=10)
    ax2.set_title(r"(b) Theorems 2.6–2.10: Higher-Girth $C_{16}(G)>0$ Exclusion Horizon", fontsize=10.5)
    ax2.grid(True, axis="x", ls=":", alpha=0.45)

    save_both(fig, "fig6_chebyshev_bessel_cascade_and_discharging")


def main() -> int:
    plot_fig1_moore_bounds()
    plot_fig2_cage_multiplicities()
    plot_fig3_2opt_seesaw()
    plot_fig4_cycle_spectra()
    plot_fig5_phase3_commutator_and_tabu()
    plot_fig6_chebyshev_bessel_cascade()
    print(f"Generated 6 publication figures (PDF + PNG) in {FIG_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
