import numpy as np


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


def tail_dep(x, y, alpha=0.9):
    n = len(x)
    u = rank_average(x) / (n + 1)
    v = rank_average(y) / (n + 1)
    mask = u > alpha
    if mask.sum() == 0:
        return np.nan
    return np.mean(v[mask] > alpha)


def permutation_pvalue(x, y, alpha=0.9, B=1000, seed=1):
    rng = np.random.default_rng(seed)
    obs = tail_dep(x, y, alpha)

    null = np.array([
        tail_dep(x, rng.permutation(y), alpha)
        for _ in range(B)
    ])
    valid_null = null[~np.isnan(null)]

    if np.isnan(obs) or len(valid_null) == 0:
        pval = np.nan
        mu = np.nan
        sd = np.nan
    else:
        pval = np.mean(valid_null >= obs)
        mu = valid_null.mean()
        sd = valid_null.std()

    return obs, mu, sd, pval


def bh_fdr(pvals):
    """
    Benjamini-Hochberg FDR correction.
    Returns adjusted q-values in original order.
    """
    pvals = np.asarray(pvals, dtype=float)
    q = np.full(len(pvals), np.nan)
    finite_mask = np.isfinite(pvals)
    if not finite_mask.any():
        return q

    finite = pvals[finite_mask]
    m = len(finite)
    order = np.argsort(finite)
    ranked = finite[order]

    finite_q = np.empty(m)
    prev = 1.0

    for i in range(m - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * m / rank
        prev = min(prev, val)
        finite_q[order[i]] = prev

    q[finite_mask] = finite_q
    return q


def scan_prime_pairs(
    pairs=((2,3),(2,5),(3,5),(5,7),(5,11)),
    Ns=(500, 1000, 2000, 5000),
    alphas=(0.85, 0.90, 0.95),
    B=500,
    mode="level",
    seed=1
):
    rows = []

    for N in Ns:
        cache = {}

        for p, q in pairs:
            if p not in cache:
                vals = np.array([vp_A(n, p) for n in range(1, N + 1)])
                if mode == "jump":
                    vals = np.diff(vals)
                cache[p] = vals

            if q not in cache:
                vals = np.array([vp_A(n, q) for n in range(1, N + 1)])
                if mode == "jump":
                    vals = np.diff(vals)
                cache[q] = vals

            X = cache[p]
            Y = cache[q]

            for alpha in alphas:
                obs, mu, sd, pval = permutation_pvalue(
                    X, Y, alpha=alpha, B=B, seed=seed
                )

                rows.append({
                    "N": N,
                    "pair": (p, q),
                    "alpha": alpha,
                    "mode": mode,
                    "obs": obs,
                    "null_mean": mu,
                    "null_std": sd,
                    "pval": pval,
                })

    qvals = bh_fdr([r["pval"] for r in rows])
    for r, qv in zip(rows, qvals):
        r["BH_qval"] = qv

    rows = sorted(rows, key=lambda r: (np.isnan(r["pval"]), r["pval"]))

    for r in rows[:20]:
        print(
            f"mode={r['mode']:5s} N={r['N']:5d} "
            f"pair={r['pair']} alpha={r['alpha']:.2f} "
            f"obs={r['obs']:.4f} null={r['null_mean']:.4f} "
            f"p={r['pval']:.4f} BH-q={r['BH_qval']:.4f}"
        )

    return rows


print("=== Raw valuation level test ===")
rows_level = scan_prime_pairs(mode="level")

print("\n=== Jump / carry-change test ===")
rows_jump = scan_prime_pairs(mode="jump")
