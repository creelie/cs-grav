"""Figures of the paper.  Run from this directory: python3 make_figures.py
Writes PDF and PNG files into ../../paper/tex/figures/."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "paper", "tex", "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 10, "mathtext.fontset": "cm"})
BLUE, RED, GREY, GREEN = "#2b5d9b", "#b8312f", "#7f7f7f", "#2e8b57"


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_signatures():
    """(s,t,k): the plane s+t-2k=10, the wall k=t, Nishino's points, the families of 27+3."""
    fig = plt.figure(figsize=(6.4, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    # plane s + t - 2k = 10 over the admissible region k <= t <= s
    T, K = np.meshgrid(np.linspace(0, 6, 13), np.linspace(0, 6, 13))
    S = 10 + 2 * K - T
    mask = K > T
    S = np.where(mask, np.nan, S)
    ax.plot_surface(S, T, K, color=BLUE, alpha=0.18, linewidth=0)
    # wall k = t
    Sg, Tg = np.meshgrid(np.linspace(8, 30, 12), np.linspace(0, 6, 12))
    ax.plot_surface(Sg, Tg, Tg, color=GREY, alpha=0.10, linewidth=0)
    # Nishino: (9+m, 1+m, m)
    m = np.arange(0, 6)
    ax.scatter(9 + m, 1 + m, m, color=BLUE, s=28, depthshade=False, label=r"$(9+m,1+m)$, $k=m$")
    # the transplanted families (9+8n+q, 1+q), every admissible k <= 1+q
    for n in (1, 2):
        for q in range(4):
            s, t = 9 + 8 * n + q, 1 + q
            ks = np.arange(0, t + 1)
            ax.plot([s] * len(ks), [t] * len(ks), ks, color=RED, lw=1.0)
            ax.scatter([s] * len(ks), [t] * len(ks), ks, color=RED, s=10, depthshade=False)
    ax.scatter([27], [3], [0], color=RED, s=60, marker="*", depthshade=False, label=r"$27+3$, every $k\leq 3$")
    ax.plot([], [], color=RED, lw=1, label=r"$(9+8n+q,1+q)$, $n=1,2$")
    ax.set_xlabel("$s$"); ax.set_ylabel("$t$"); ax.set_zlabel("$k$")
    ax.set_xlim(8, 30); ax.set_ylim(0, 6); ax.set_zlim(0, 6)
    ax.set_xticks(range(10, 31, 5))
    ax.view_init(elev=20, azim=-62)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.text(14, 4.6, 4.9, r"$s+t-2k=10$", color=BLUE, fontsize=9)
    ax.text(24, 5.5, 5.8, r"$k=t$", color=GREY, fontsize=9)
    save(fig, "fig_signatures")


def fig_spinorcount():
    """Spinor dimension against 2(d-2): they meet only at d = 3, 4, 6, 10."""
    d = np.arange(3, 19)
    spin = np.array([2 ** (int(np.ceil(x / 2)) - 1) for x in d])
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    ax.semilogy(d, spin, "o-", color=RED, ms=4, label=r"least $N$: $2^{\lceil d/2\rceil-1}$")
    ax.semilogy(d, 2 * (d - 2), "s-", color=BLUE, ms=4, label=r"$N=2(d-2)$ (cubic identity)")
    for x in (3, 4, 6, 10):
        ax.axvline(x, color=GREY, lw=0.6, ls=":")
    hits = [(3, 2), (4, 4), (6, 8), (10, 16)]
    ax.scatter([h[0] for h in hits], [h[1] for h in hits], s=70, facecolors="none", edgecolors="k", zorder=5)
    ax.set_xlabel("$d$"); ax.set_ylabel("$N$")
    ax.set_xticks(d)
    ax.legend(frameon=False, fontsize=9)
    save(fig, "fig_spinorcount")


def fig_lightcone():
    """R^{2,1}: null vectors on the cone; two independent null vectors are never orthogonal."""
    fig = plt.figure(figsize=(5.4, 5.0))
    ax = fig.add_subplot(111, projection="3d")
    th = np.linspace(0, 2 * np.pi, 80)
    r = np.linspace(-1.2, 1.2, 30)
    TH, R = np.meshgrid(th, r)
    ax.plot_surface(R * np.cos(TH), R * np.sin(TH), R, color=BLUE, alpha=0.15, linewidth=0)
    a, b = 0.3, 2.1
    for ang, col, lab in [(a, RED, r"$n$"), (b, GREEN, r"$n'$")]:
        v = np.array([np.cos(ang), np.sin(ang), 1.0])
        ax.plot([0, 1.2 * v[0]], [0, 1.2 * v[1]], [0, 1.2 * v[2]], color=col, lw=2)
        ax.text(*(1.3 * v), lab, color=col, fontsize=12)
    # the plane spanned by n, n': it cuts the cone in exactly these two lines
    v1 = np.array([np.cos(a), np.sin(a), 1.0]); v2 = np.array([np.cos(b), np.sin(b), 1.0])
    P = [0 * v1, 1.2 * v1, 1.2 * (v1 + v2) / 2 * 1.0, 1.2 * v2]
    ax.add_collection3d(Poly3DCollection([P], color=GREY, alpha=0.25))
    ax.text(0, 0, -1.55, r"$n\cdot n'=\cos(\theta-\theta')-1<0$", fontsize=9, ha="center")
    ax.set_xlabel("$x^1$"); ax.set_ylabel("$x^2$"); ax.set_zlabel("$x^0$")
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_zlim(-1.2, 1.2)
    ax.view_init(elev=18, azim=-50)
    save(fig, "fig_lightcone")


def fig_multiplets():
    """log2 of the least massless multiplet against log2 of the vector multiplet, D = 10+8n."""
    fig = plt.figure(figsize=(6.0, 4.6))
    ax = fig.add_subplot(111, projection="3d")
    ns = np.arange(0, 4)
    S = 2 ** (4 + 4 * ns)
    least = np.log2(S / 4)                      # log2 log2 of 2^{S/4}
    vec = np.log2(np.log2((8 + 8 * ns) + S / 2))  # log2 log2 of bosons + fermions
    w = 0.35
    ax.bar3d(ns - w, np.zeros_like(ns) , np.zeros_like(ns), w, 0.5, vec, color=BLUE, alpha=0.85, shade=True)
    ax.bar3d(ns, np.ones_like(ns), np.zeros_like(ns), w, 0.5, least, color=RED, alpha=0.85, shade=True)
    ax.set_xticks(ns); ax.set_xticklabels([f"$D={10+8*n}$" for n in ns], fontsize=8)
    ax.set_yticks([0.25, 1.25]); ax.set_yticklabels(["vector", "least"], fontsize=8)
    ax.set_zlabel(r"$\log_2\log_2(\#\mathrm{states})$")
    ax.view_init(elev=22, azim=-58)
    save(fig, "fig_multiplets")


def fig_residuals():
    """Relative residual of the absorption equation for the signatures checked."""
    data = [((9, 1), 0, 3.0e-15), ((10, 2), 1, 1.6e-15), ((11, 3), 2, 1.5e-15), ((12, 4), 3, 1.4e-15),
            ((13, 5), 4, 2.6e-15), ((6, 6), 1, 1.4e-15), ((7, 7), 2, 2.8e-15),
            ((17, 1), 0, 9.8e-1), ((18, 2), 1, 9.9e-1), ((19, 3), 0, 1.0), ((19, 3), 2, 9.9e-1),
            ((20, 4), 1, 1.0), ((20, 4), 3, 9.8e-1)]
    fig = plt.figure(figsize=(6.2, 4.8))
    ax = fig.add_subplot(111, projection="3d")
    for (s, t), k, res in data:
        h = np.log10(res) + 16
        col = BLUE if res < 1e-10 else RED
        ax.bar3d(s - 0.35, k - 0.35, 0, 0.7, 0.7, h, color=col, alpha=0.8, shade=True)
    ax.set_xlabel("$s$"); ax.set_ylabel("$k$")
    ax.set_xticks(range(6, 21, 2)); ax.set_yticks(range(0, 5))
    ax.set_zlabel(r"$16+\log_{10}(\mathrm{residual})$")
    ax.view_init(elev=26, azim=-60)
    save(fig, "fig_residuals")


if __name__ == "__main__":
    fig_signatures(); fig_spinorcount(); fig_lightcone(); fig_multiplets(); fig_residuals()
    print("figures written to", os.path.normpath(OUT))
