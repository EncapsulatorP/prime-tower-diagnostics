# Ranked Prime-Power Towers, Drift Diagnostics, and Lyapunov Splits

## Scope

This repository is a numerical notebook in script form. Its scope is limited and specific:

- study the prime-power decomposition of a structured coefficient family,
- examine whether the decomposition heights show nontrivial dependence,
- compare additive and multiplicative realizations of the same base quantity `6`,
- use tower maps, cocycles, and a Borwein-like integral as drift diagnostics,
- record the places where the construction may fail or become biased.

It is not a proof repository. Most of what is here should be read as a sequence of controlled experiments, with a few sharper heuristic claims collected in [Lemma.md](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/Lemma.md:1) and the execution status summarized in [results.md](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/results.md:1).

## Central Object: A Ranked Prime-Power Tower

The arithmetic object behind the prime experiments is

```text
A_n = (6n)! / ((3n)! (n!)^3).
```

For each prime `p`, define the `p`-height

```text
h_p(n) = v_p(A_n).
```

Then the vector

```text
(h_2(n), h_3(n), h_5(n), h_7(n), ...)
```

is the prime-power decomposition tower at level `n`. In this repository, “ranked” means two related things:

- we compare the heights across an ordered list of primes, usually small primes first,
- we compare the same height coordinates across different index sets, especially all `n` versus prime indices `n = q_k`.

So the tower height is not one scalar. It is a coordinatewise stack of `p`-adic heights, and the experiments look for dependence among those coordinates.

For this specific `A_n`, Legendre’s formula gives

```text
h_p(n) = v_p((6n)!) - v_p((3n)!) - 3 v_p(n!),
```

and the fast version used in `scripts/prime_index_dependence_test.py` rewrites this in terms of base-`p` digit sums. This matters conceptually: the tower heights are controlled by carry structure, and the factor `6` enters the carry structure directly through `3n` and `6n`.

## Why `6` Is the Structural Pivot

The number `6` appears in two different ways in this project.

First, it is arithmetic. The coefficient family contains `6n`, `3n`, and `n`, so `6 = 2 * 3` is built into the valuation identities. This makes the `2`- and `3`-channels the most immediate places where dependence or carry interactions can accumulate.

Second, it is dynamical. The cocycle toy models compare two ways of realizing the same gross scale:

- additive path: `6 = 3 + 3`, implemented as a doubled base-`3` channel,
- multiplicative path: `6 = 3 * 2`, implemented as a single base-`6` channel.

This is not just notation. The two realizations have different local derivatives, different average log-growth, and therefore different Lyapunov behavior. The project’s main warning signal is that these two realizations may not stay balanced as the tower grows.

## What the Prime Scripts Actually Test

The prime-dependence scripts do not prove a theorem about all primes. They ask narrower questions:

- do the height coordinates `(h_p(n), h_q(n))` exhibit copula or tail dependence?
- does the prime-indexed subsequence `n = q_k` behave differently from random index controls?
- can random subsets reproduce the observed cross-prime dependence?
- after residue matching, is there still a gap-linked effect tied to prime spacing?

The current run summary in [results.md](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/results.md:47) indicates that the prime-indexed subsequence does retain extra dependence beyond a naive random control, and that part of the effect survives a residue-matched control. That should be read as evidence for structure, not as a finished classification.

## Borwein-Like Integration as a Drift Diagnostic

The script [scripts/borwein_integral_logistic_6.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/borwein_integral_logistic_6.py:1) builds a ladder of scales

```text
{2, 3}
{2, 3, 6}
{2, 3, 6, sqrt(6)}
{2, 3, 6, sqrt(6), 1 + sqrt(6)}
{2, 3, 6, sqrt(6), 1 + sqrt(6), 1 / (1 + sqrt(6))}
```

and evaluates

```text
I_L = ∫_0^∞ ∏ sinc(x / a_j) dx
```

for the scale set at level `L`.

In this repository, that integral is not treated as a theorem-generator. It is treated as a drift probe. The quantities

- `|2 I_L / π - 1|`,
- c10 signature distances,
- jumps between successive levels,

are used to see whether adding the `6`-derived scales nudges the construction toward or away from a `π/2`-like reference.

Operationally, the question is:

> when the tower is extended by scales forced by `6`, does the numerical signature stay coherent, or does it jump?

The current run reports the largest c10 jump at `L3 -> L4`, which is consistent with a drift warning rather than a stability certificate.

## The Main Risk: Nonuniform Cocycle Growth

The most important caveat in the whole project is the cocycle split.

The picture in `spinor_split.png` and the argument in [Lemma.md](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/Lemma.md:1) both point in the same direction:

- the additive realization of `6` and the multiplicative realization of `6` share the same base orbit support,
- but the multiplicative path appears to accumulate larger log-expansion on average,
- so the two channels peel apart in height.

That matters because many tower constructions quietly assume that two ways of packaging the same scale will stay interchangeable after many iterates. If the Lyapunov exponents differ, that interchangeability degrades with depth.

In practical terms, a persistent split can negatively affect the construction in at least three ways:

- it can bias tower height toward one channel, so the decomposition no longer reflects a balanced arithmetic superposition,
- it can make glued layers drift apart, meaning that equal-looking boundary data does not remain equal after transport,
- it can corrupt any inference that assumes approximate uniformity of local growth across channels.

This negative impact is not proved for the prime tower itself. The prime tower and the cocycle toy model are not the same object. But the toy model is relevant precisely because it isolates a mechanism that could contaminate the arithmetic construction: one channel grows faster everywhere, not just on rare spikes.

That is why the right panel in `spinor_split.png` is important. The histogram shift is interpreted here as a uniformity failure, not merely as a larger variance. Uniformity failures are harder to average away, so they are more dangerous for a tower construction than isolated outliers.

## Relation to the Lemma

[Lemma.md](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/Lemma.md:1) states the intended sign claim:

```text
λ_mul > λ_add
```

for the base-`6` comparison, under mild empirical assumptions on orbit occupation.

The README should be read conservatively here:

- the lemma is the conceptual spine of the spinor split story,
- the repository provides numerical support for that sign,
- the sign matters because it predicts systematic peeling in tower height,
- but the transfer from the toy cocycle to the prime-power tower is still heuristic.

So the repository’s position is not “the prime tower theorem is proved.” It is:

> if the multiplicative channel really has uniformly larger growth, then any construction that treats additive and multiplicative `6`-channels as neutral substitutes is at risk.

## Why `π` Looks Special in the Tower Plots

The tower maps in [scripts/tower_maps.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/tower_maps.py:1) compare several bases under

```text
T(y) = frac(b^y)
```

and under the Gauss-like map

```text
T(y) = frac((1 / y)^(1 / y)).
```

For the pure positional tower `frac(b^y)`, `π` is visually richer than `sqrt(2)` and `φ`, and usually richer than `e`, for a simple reason: on `y in [0, 1]`, the quantity `π^y` crosses more integer thresholds than bases closer to `1`. More threshold crossings mean more wraps under `frac`, and more wraps mean more visible branches in the delayed embedding.

This does not make `π` “better” in any absolute sense. It makes it less collapsed for this specific positional tower diagnostic.

In the Pisot-vs-`π` comparison, `π` is used as a reference precisely because it keeps a visibly multi-sheeted embedding where several other candidate bases collapse toward thinner supports.

## The Problem with `sqrt(2)` and the Liouville Constant

The weakness of `sqrt(2)` and the Liouville constant in the pure positional tower is not mysterious once the map is fixed.

For `sqrt(2)`:

- the base is only slightly above `1`,
- `b^y` ranges over a narrow interval,
- there are very few integer threshold crossings,
- so `frac(b^y)` has little wrapping and quickly collapses toward a thin attractor.

For the Liouville constant used here:

- the base is less than `1`,
- the map does even less wrapping,
- the orbit tends to stay confined to a very thin set,
- so its famous Diophantine extremality does not translate into rich positional tower geometry in this model.

This is an important methodological point. Exceptional number-theoretic status by itself does not guarantee a rich tower image. The chosen map can suppress that richness almost completely.

## What to Take Seriously and What Not to Overclaim

Take seriously:

- the valuation tower has observable cross-prime structure,
- the prime-indexed subsequence is not obviously random relative to simple controls,
- the multiplicative and additive `6`-channels need not grow at the same rate,
- a Borwein-like integral can be used as a practical drift monitor.

we do not overclaim:

- that the cocycle lemma has already been transferred rigorously to the prime tower,
- that the Borwein-like integral proves the drift rather than diagnosing it,
- that `π` is canonically privileged rather than just less collapsed in these tower maps,
- that `sqrt(2)` or the Liouville constant are intrinsically defective rather than poorly matched to this specific positional model.

## Repository Map

- [scripts/prime_index_dependence_test.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/prime_index_dependence_test.py:1): prime-indexed dependence versus random and residue-matched controls.
- [scripts/5_copula_test.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/5_copula_test.py:1): small copula / tail-dependence tests on valuation pairs.
- [scripts/bonferroni_failed_readressing.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/bonferroni_failed_readressing.py:1): scan across prime pairs, sizes, and tail levels.
- [scripts/borwein_integral_logistic_6.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/borwein_integral_logistic_6.py:1): Borwein-like drift diagnostic for the `6`-scale ladder.
- [scripts/spinor_split_plot.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/spinor_split_plot.py:1): the additive versus multiplicative Lyapunov split visualization.
- [scripts/spinor_cocycles.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/spinor_cocycles.py:1): Newton / Collatz-like / Pisot-like cocycle comparison.
- [scripts/tower_maps.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/tower_maps.py:1): pure positional and Gauss-like tower maps, including `π`, `sqrt(2)`, and the Liouville constant.
- [scripts/cf_recurrence.py](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/scripts/cf_recurrence.py:1): continued-fraction recurrence plots for `sqrt(2)`, `e`, and `π`.

## Running the Scripts

From the repository root:

```bash
python3 scripts/prime_index_dependence_test.py
python3 scripts/borwein_integral_logistic_6.py
python3 scripts/spinor_split_plot.py
python3 scripts/tower_maps.py
python3 scripts/cf_recurrence.py
python3 scripts/pi_e_copulas.py
python3 scripts/spinor_cocycles.py
python3 scripts/5_copula_test.py
python3 scripts/bonferroni_failed_readressing.py
```

Generated outputs and verification logs are summarized in [results.md](/mnt/c/Users/inger/Downloads/logistic_collatz_pisot_newton_spinors/results.md:1) with GPT and Claude tools support.
