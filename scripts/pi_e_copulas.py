import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str((Path("outputs") / ".mplconfig").resolve()),
)

import math

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scripts.plot_io import save_figure
except ModuleNotFoundError:
    from plot_io import save_figure


# -----------------------------
# 1. p-adic valuation utilities
# -----------------------------

def vp_factorial(n, p):
    """Legendre formula: exponent of prime p in n!."""
    total = 0
    while n:
        n //= p
        total += n
    return total


def vp_chudnovsky_block(n, p):
    """
    A_n = (6n)! / ((3n)! * (n!)^3)

    Returns v_p(A_n), without constructing huge integers.
    """
    return (
        vp_factorial(6 * n, p)
        - vp_factorial(3 * n, p)
        - 3 * vp_factorial(n, p)
    )


# -----------------------------
# 2. Rank / empirical copula tools
# -----------------------------

def rankdata_average(x):
    """Average ranks, like scipy.stats.rankdata(method='average')."""
    x = np.asarray(x)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)

    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and x[order[j]] == x[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + j + 1) / 2.0
        i = j

    return ranks


def pearson_corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    aa = a - a.mean()
    bb = b - b.mean()
    den = np.sqrt((aa @ aa) * (bb @ bb))
    return np.nan if den == 0 else (aa @ bb) / den


def spearman_corr(x, y):
    """Spearman correlation = Pearson correlation of ranks."""
    return pearson_corr(rankdata_average(x), rankdata_average(y))


def empirical_tail_dependence(u, v, alpha=0.90):
    """
    Crude upper-tail dependence estimate:
    P(V > alpha | U > alpha)
    """
    mask = u > alpha
    if mask.sum() == 0:
        return np.nan
    return np.mean(v[mask] > alpha)


# -----------------------------
# 3. Main experiment
# -----------------------------

def prime_copula_experiment(p=2, q=3, N=1000, alpha=0.90, seed=123, save_path=None, show=False):
    n_values = np.arange(1, N + 1)

    X = np.array([vp_chudnovsky_block(n, p) for n in n_values])
    Y = np.array([vp_chudnovsky_block(n, q) for n in n_values])

    # Empirical copula coordinates
    U = rankdata_average(X) / (N + 1)
    V = rankdata_average(Y) / (N + 1)

    rho = spearman_corr(X, Y)
    tail = empirical_tail_dependence(U, V, alpha=alpha)

    # Null: destroy dependence but preserve margins
    rng = np.random.default_rng(seed)
    Y_perm = rng.permutation(Y)
    V_perm = rankdata_average(Y_perm) / (N + 1)

    rho_null = spearman_corr(X, Y_perm)
    tail_null = empirical_tail_dependence(U, V_perm, alpha=alpha)

    print(f"Coefficient family: A_n = (6n)! / ((3n)! (n!)^3)")
    print(f"N = {N}")
    print(f"Prime pair: p={p}, q={q}")
    print()
    print(f"Spearman rho observed: {rho:.6f}")
    print(f"Spearman rho null:     {rho_null:.6f}")
    print()
    print(f"Upper-tail dependence observed, alpha={alpha}: {tail:.6f}")
    print(f"Upper-tail dependence null:                    {tail_null:.6f}")

    # Empirical copula scatter
    fig = plt.figure(figsize=(6, 6))
    plt.scatter(U, V, s=8, alpha=0.45)
    plt.xlabel(f"rank-copula U for v_{p}(A_n)")
    plt.ylabel(f"rank-copula V for v_{q}(A_n)")
    plt.title(f"Empirical prime-valuation copula: p={p}, q={q}")
    plt.grid(True, alpha=0.25)

    if save_path is not None:
        save_figure(fig, save_path)
    if show:
        plt.show()
    plt.close(fig)

    return {
        "p": p,
        "q": q,
        "N": N,
        "spearman_observed": rho,
        "spearman_null": rho_null,
        "tail_observed": tail,
        "tail_null": tail_null,
        "X": X,
        "Y": Y,
        "U": U,
        "V": V,
    }


def main():
    pairs = [(2, 3), (2, 5), (3, 5), (5, 7), (7, 11), (2, 101), (61, 1103), (7, 1103)]
    outdir = Path("outputs") / "prime_copulas"

    print("\n--- Prime copula panels ---")
    for p, q in pairs:
        path = outdir / f"p{p}_q{q}.png"
        prime_copula_experiment(p=p, q=q, N=1000, alpha=0.90, save_path=path)
        print(f"Saved {path}")
        print()


if __name__ == "__main__":
    main()
