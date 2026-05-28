import numpy as np
import math


def vp_factorial(n, p):
    total = 0
    while n:
        n //= p
        total += n
    return total


def vp_A(n, p):
    return (
        vp_factorial(6 * n, p)
        - vp_factorial(3 * n, p)
        - 3 * vp_factorial(n, p)
    )


def rank_average(x):
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


def tail_dep(x, y, alpha=0.90):
    n = len(x)
    u = rank_average(x) / (n + 1)
    v = rank_average(y) / (n + 1)

    mask = u > alpha
    if mask.sum() == 0:
        return np.nan

    return np.mean(v[mask] > alpha)


def permutation_test_tail(p=2, q=5, N=2000, alpha=0.90, B=1000, seed=1):
    rng = np.random.default_rng(seed)

    X = np.array([vp_A(n, p) for n in range(1, N + 1)])
    Y = np.array([vp_A(n, q) for n in range(1, N + 1)])

    observed = tail_dep(X, Y, alpha=alpha)

    nulls = []
    for _ in range(B):
        nulls.append(tail_dep(X, rng.permutation(Y), alpha=alpha))

    nulls = np.array(nulls)
    p_value = np.mean(nulls >= observed)

    print(f"p={p}, q={q}, N={N}, alpha={alpha}")
    print(f"Observed tail dependence: {observed:.6f}")
    print(f"Null mean: {nulls.mean():.6f}")
    print(f"Null std:  {nulls.std():.6f}")
    print(f"Permutation p-value: {p_value:.4f}")

    return observed, nulls


for pair in [(2, 3), (2, 5), (3, 5), (5, 7), (5, 11)]:
    permutation_test_tail(*pair, N=2000, alpha=0.90, B=500)
    print()