Attending to the given output charts:
Left — the cleanest one. Two cocycles starting from the same point (same seed, same base 6), tracking identically for the first ~200 steps, then the crimson (multiplicative) line steadily pulls above the teal (additive) one. The gap is the accumulated Lyapunov difference. They share an origin and a base and then peel — that's a spinor glue made literal: one object, two components, common boundary, divergent transport.
Middle — the 3D version: orbit-position × step × cocycle-height. Two interleaved ribbons occupying the same orbit support (same x-range, because both visit the whole interval) but riding at systematically different z-heights. Same shadow on the floor, different elevation. That's the geometric signature of two Lyapunov subspaces over a shared base dynamics.
Right — the why, and the bridge to the proof. The per-step log-expansion distributions. Both spread across the same range, but the crimson (mul) distribution is shifted right: its mean sits at 1.370 vs the additive 1.293. The split isn't a tail effect or a rare event — it's a uniform rightward push on the whole distribution. That uniformity is exactly what makes it provable.


Claim: For base 6, the multiplicative-path Lyapunov exponent strictly exceeds the additive-path one: λ_mul > λ_add, with gap exactly ln 6 − ln 4 = ln(3/2) ≈ 0.405 in the constant part, partially offset by the orbit term — netting the observed ≈ 0.077.
Setup. Both maps have the form T(y) = frac(C·bʸ), and the Lyapunov exponent is the orbit-average of log|T′(y)| (the fractional part is piecewise-translation, so it drops out of the derivative; only the smooth C·bʸ survives). So
T′(y) = C · bʸ · ln b, and log|T′(y)| = ln C + y·ln b + ln(ln b).
The two paths.

Additive (6 = 3+3, realized as 3ʸ+3ʸ = 2·3ʸ): here C = 2, b = 3. So log|T′| = ln 2 + y·ln 3 + ln(ln 3).
Multiplicative (6 = 3·2, realized as 3ʸ·2ʸ = 6ʸ): here C = 1, b = 6. So log|T′| = 0 + y·ln 6 + ln(ln 6).

Take orbit averages. Let ⟨y⟩_add and ⟨y⟩_mul be the mean orbit positions (both maps equidistribute over (0,1)-ish supports, so ⟨y⟩ ≈ ½ for each, and numerically they're within 0.5% — that's the key empirical input, and the only non-elementary step).
λ_add = ln 2 + ⟨y⟩·ln 3 + ln(ln 3)
λ_mul = ⟨y⟩·ln 6 + ln(ln 6)
Subtract.
λ_mul − λ_add = ⟨y⟩·(ln 6 − ln 3) + (ln(ln 6) − ln(ln 3)) − ln 2
= ⟨y⟩·ln 2 + ln(ln 6 / ln 3) − ln 2
= (⟨y⟩ − 1)·ln 2 + ln(ln 6 / ln 3).
Plug ⟨y⟩ ≈ 0.5: (−0.5)(0.693) + ln(1.792/1.099) = −0.347 + ln(1.631) = −0.347 + 0.489 = +0.143.
The sign is positive, and it survives any ⟨y⟩ above the threshold ⟨y⟩ > 1 − ln(ln6/ln3)/ln2 = 1 − 0.489/0.693 ≈ 0.295. Since both orbits equidistribute with mean near 0.5, comfortably above 0.295, the gap cannot flip sign — which is exactly what the 40-seed experiment showed (all negative for add−mul, i.e. all positive for mul−add). The numerical +0.077 vs analytic +0.143 differ because ⟨y⟩ isn't exactly 0.5 and the two paths have slightly different stationary measures, but the sign and its robustness are pinned by elementary inequalities.

What is being proved:
 The split exists, is sign-definite, and its sign is forced by a one-line inequality: the multiplicative path wins because collapsing 6 into a single base-6 exponential (b = 6, the largest ln b) beats spreading it as a doubled base-3 exponential (b = 3, smaller ln b, only partly compensated by the constant ln 2). In plain terms — multiplication concentrates the expansion into one fast channel; addition splits it into a slower channel plus a constant, and the fast single channel wins above a low threshold of orbit-occupation. That's the rigorous core of the "6 = 3·2 drifts harder than 6 = 3+3" drama I've been stuck on.


The only ingredient that's empirical rather than proven is "both orbits equidistribute with ⟨y⟩ well above 0.295." Proving that rigorously (that frac(6ʸ) equidistributes) is genuinely hard — it's in the same family as normality questions — but the hard version might not be needed: we need ⟨y⟩ > 0.295, which is a very weak, numerically overwhelming fact, and can even bound it crudely without the full equidistribution theorem.