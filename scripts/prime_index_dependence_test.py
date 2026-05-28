import argparse
from math import isqrt

import numpy as np


def primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : n + 1 : p] = b"\x00" * (((n - start) // p) + 1)
    return np.flatnonzero(np.frombuffer(sieve, dtype=np.uint8))


def digit_sum_base(n, base):
    total = 0
    while n:
        total += n % base
        n //= base
    return total


def vp_chudnovsky_block(n, p):
    """
    For A_n = (6n)! / ((3n)! * (n!)^3), Legendre's formula simplifies to

        v_p(A_n) = (3 s_p(n) + s_p(3n) - s_p(6n)) / (p - 1),

    so the valuation depends only on base-p digit sums.
    """
    return (
        3 * digit_sum_base(n, p)
        + digit_sum_base(3 * n, p)
        - digit_sum_base(6 * n, p)
    ) // (p - 1)


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


def spearman_corr(x, y):
    rx = rank_average(x)
    ry = rank_average(y)
    rx -= rx.mean()
    ry -= ry.mean()
    return (rx @ ry) / np.sqrt((rx @ rx) * (ry @ ry))


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    qvals = np.empty(m, dtype=float)
    prev = 1.0

    for i in range(m - 1, -1, -1):
        rank = i + 1
        prev = min(prev, ranked[i] * m / rank)
        qvals[order[i]] = prev

    return qvals


def pairwise_corr_vector(matrix):
    corr = np.corrcoef(matrix, rowvar=False)
    upper = np.triu_indices_from(corr, k=1)
    return corr[upper]


def offdiag_corr_mean(matrix):
    return pairwise_corr_vector(matrix).mean()


def precompute_valuations(N, bases):
    vals = np.empty((N + 1, len(bases)), dtype=np.int16)
    vals[0, :] = 0
    for n in range(1, N + 1):
        for j, base in enumerate(bases):
            vals[n, j] = vp_chudnovsky_block(n, base)
    return vals


def sample_residue_matched(rng, pools, counts):
    parts = [
        rng.choice(pools[residue], size=counts[residue], replace=False)
        for residue in sorted(counts)
    ]
    return np.sort(np.concatenate(parts))


def monte_carlo_controls(valuations, all_indices, prime_indices, residue_modulus, B, seed):
    rng = np.random.default_rng(seed)
    K = len(prime_indices)
    residues = prime_indices % residue_modulus
    counts = {
        residue: int(np.sum(residues == residue))
        for residue in np.unique(residues)
    }
    pools = {
        residue: all_indices[all_indices % residue_modulus == residue]
        for residue in counts
    }

    num_bases = valuations.shape[1]
    uniform_means = np.empty((B, num_bases), dtype=float)
    matched_means = np.empty((B, num_bases), dtype=float)
    uniform_corr = np.empty(B, dtype=float)
    matched_corr = np.empty(B, dtype=float)
    uniform_gap = np.empty((B, num_bases), dtype=float)
    matched_gap = np.empty((B, num_bases), dtype=float)

    for b in range(B):
        uniform_idx = np.sort(rng.choice(all_indices, size=K, replace=False))
        matched_idx = sample_residue_matched(rng, pools, counts)

        uniform_vals = valuations[uniform_idx]
        matched_vals = valuations[matched_idx]

        uniform_means[b] = uniform_vals.mean(axis=0)
        matched_means[b] = matched_vals.mean(axis=0)
        uniform_corr[b] = offdiag_corr_mean(uniform_vals)
        matched_corr[b] = offdiag_corr_mean(matched_vals)

        uniform_gaps = np.diff(uniform_idx)
        matched_gaps = np.diff(matched_idx)
        for j in range(num_bases):
            uniform_gap[b, j] = spearman_corr(uniform_vals[:-1, j], uniform_gaps)
            matched_gap[b, j] = spearman_corr(matched_vals[:-1, j], matched_gaps)

    return {
        "counts": counts,
        "uniform_means": uniform_means,
        "matched_means": matched_means,
        "uniform_corr": uniform_corr,
        "matched_corr": matched_corr,
        "uniform_gap": uniform_gap,
        "matched_gap": matched_gap,
    }


def empirical_pvalue(null, observed, alternative="greater"):
    null = np.asarray(null, dtype=float)
    if alternative == "greater":
        extreme = np.sum(null >= observed)
    elif alternative == "less":
        extreme = np.sum(null <= observed)
    else:
        extreme = np.sum(np.abs(null) >= abs(observed))
    return (extreme + 1.0) / (len(null) + 1.0)


def summarize_means(bases, observed, uniform_means, matched_means):
    print("Mean valuations on prime indices versus controls")
    print("base  prime_mean  uniform_mean  uniform_z  matched_mean  matched_z")
    for j, base in enumerate(bases):
        uniform_mu = uniform_means[:, j].mean()
        uniform_sd = uniform_means[:, j].std(ddof=1)
        matched_mu = matched_means[:, j].mean()
        matched_sd = matched_means[:, j].std(ddof=1)
        uniform_z = (observed[j] - uniform_mu) / uniform_sd
        matched_z = (observed[j] - matched_mu) / matched_sd
        print(
            f"{base:>4d}  {observed[j]:10.4f}  {uniform_mu:12.4f}  "
            f"{uniform_z:9.2f}  {matched_mu:12.4f}  {matched_z:9.2f}"
        )
    print()


def summarize_dependence(bases, observed_corr, uniform_corr, matched_corr):
    uniform_mu = uniform_corr.mean()
    uniform_sd = uniform_corr.std(ddof=1)
    matched_mu = matched_corr.mean()
    matched_sd = matched_corr.std(ddof=1)
    p_uniform = empirical_pvalue(uniform_corr, observed_corr, alternative="greater")
    p_matched = empirical_pvalue(matched_corr, observed_corr, alternative="greater")

    print("Cross-base dependence: mean off-diagonal Pearson correlation")
    print(f"Observed prime subsequence: {observed_corr:.6f}")
    print(
        f"Uniform control:          {uniform_mu:.6f} +/- {uniform_sd:.6f} "
        f"(p={p_uniform:.4f})"
    )
    print(
        f"Residue-matched control:  {matched_mu:.6f} +/- {matched_sd:.6f} "
        f"(p={p_matched:.4f})"
    )
    print()


def summarize_gap_correlations(bases, observed_gap, uniform_gap, matched_gap):
    p_uniform = [
        empirical_pvalue(uniform_gap[:, j], observed_gap[j], alternative="greater")
        for j in range(len(bases))
    ]
    p_matched = [
        empirical_pvalue(matched_gap[:, j], observed_gap[j], alternative="greater")
        for j in range(len(bases))
    ]
    q_uniform = bh_fdr(p_uniform)
    q_matched = bh_fdr(p_matched)

    print("Gap-linked dependence: Spearman(value at p_k, next prime gap)")
    print(
        "base  observed   uniform_mu  uniform_p  uniform_q  "
        "matched_mu  matched_p  matched_q"
    )
    for j, base in enumerate(bases):
        print(
            f"{base:>4d}  {observed_gap[j]:8.4f}  "
            f"{uniform_gap[:, j].mean():10.4f}  {p_uniform[j]:9.4f}  "
            f"{q_uniform[j]:9.4f}  {matched_gap[:, j].mean():10.4f}  "
            f"{p_matched[j]:9.4f}  {q_matched[j]:9.4f}"
        )
    print()


def run_experiment(N, warmup, bases, residue_modulus, B, seed):
    primes = primes_upto(N)
    prime_indices = primes[primes > warmup]
    all_indices = np.arange(warmup + 1, N + 1)
    valuations = precompute_valuations(N, bases)

    observed_vals = valuations[prime_indices]
    observed_means = observed_vals.mean(axis=0)
    observed_corr = offdiag_corr_mean(observed_vals)
    observed_gap = np.array(
        [
            spearman_corr(observed_vals[:-1, j], np.diff(prime_indices))
            for j in range(len(bases))
        ]
    )

    controls = monte_carlo_controls(
        valuations=valuations,
        all_indices=all_indices,
        prime_indices=prime_indices,
        residue_modulus=residue_modulus,
        B=B,
        seed=seed,
    )

    residue_counts = {
        int(residue): int(count)
        for residue, count in sorted(controls["counts"].items())
    }

    print(
        f"Prime-indexed Chudnovsky block test: N={N}, warmup>{warmup}, "
        f"{len(prime_indices)} prime indices, B={B}"
    )
    print(f"Valuation bases: {list(bases)}")
    print(f"Prime residue counts mod {residue_modulus}: {residue_counts}")
    print()

    summarize_means(
        bases,
        observed_means,
        controls["uniform_means"],
        controls["matched_means"],
    )
    summarize_dependence(
        bases,
        observed_corr,
        controls["uniform_corr"],
        controls["matched_corr"],
    )
    summarize_gap_correlations(
        bases,
        observed_gap,
        controls["uniform_gap"],
        controls["matched_gap"],
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=50000)
    parser.add_argument("--warmup", type=int, default=100)
    parser.add_argument("--B", type=int, default=300)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--modulus", type=int, default=30)
    parser.add_argument(
        "--bases",
        type=int,
        nargs="+",
        default=[2, 3, 5, 7, 11],
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_experiment(
        N=args.N,
        warmup=args.warmup,
        bases=tuple(args.bases),
        residue_modulus=args.modulus,
        B=args.B,
        seed=args.seed,
    )
