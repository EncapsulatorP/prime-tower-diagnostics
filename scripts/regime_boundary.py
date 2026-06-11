"""
regime_boundary.py
------------------
Lyapunov deflection sign across composite bases N = 4..50 for the positional
power tower  T(y) = frac(b^y).

For each composite N, compare two channels sharing the same value N:
  multiplicative channel : T_mul(y) = frac( N^y )           (b_mul = N,  C_mul = 1)
  additive channel       : T_add(y) = frac( C * b_add^y )   chosen to match the
                                                            6 = 3+3 -> 2*3^y pattern:
       N even -> C = 2, b_add = N/2
       N odd  -> C = smallest prime factor p, b_add = N/p

For each tower step y the derivative is  T'(y) = (C * b^y) * ln(b),
so  log|T'(y)| = ln(C) + y*ln(b) + ln(ln(b)).
The Lyapunov exponent is the orbit-average of log|T'|; the "deflection" is
lambda_mul - lambda_add.

Sign lemma (from the orbit-average decomposition):
    sign(gap) > 0   iff   <y> > 1 - ln(ln(N)/ln(b_add)) / ln(C)
where <y> is the orbit mean (empirically ~0.44-0.46 for these maps).

Output:
  regime_boundary.png  -- two-panel figure (gap vs N, threshold vs N)
  prints a small table to stdout.

Dependencies: numpy, mpmath, matplotlib, sympy.
Tested on Python 3.12.
"""

import sys
try:
    import numpy as np
    import mpmath as mp
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sympy import factorint, isprime
except ImportError as e:
    sys.exit(f"missing dependency: {e}.  pip install numpy mpmath matplotlib sympy")

mp.mp.dps = 120  # working precision; 120 dps is plenty for ~1200-step orbits


# ---------- core map and Lyapunov estimators ----------

def frac(y):
    return y - mp.floor(y)

def lyap_mul(b, y0, n):
    """Lyapunov of T(y)=frac(b^y) at seed y0 for n steps."""
    s = mp.mpf(0); c = 0
    y = frac(mp.mpf(y0))
    lnb = mp.log(b)
    for _ in range(n):
        if y > 0:
            d = mp.power(b, y) * lnb
            s += mp.log(abs(d)); c += 1
        y = frac(mp.power(b, y))
    return float(s / c) if c else float("nan")

def lyap_add(C, b, y0, n):
    """Lyapunov of T(y)=frac(C * b^y)."""
    s = mp.mpf(0); c = 0
    y = frac(mp.mpf(y0))
    lnb = mp.log(b); mpC = mp.mpf(C)
    for _ in range(n):
        if y > 0:
            d = mpC * mp.power(b, y) * lnb
            s += mp.log(abs(d)); c += 1
        y = frac(mpC * mp.power(b, y))
    return float(s / c) if c else float("nan")


# ---------- sweep composite bases 4..50 ----------

SEEDS_PER_BASE = 16
ORBIT_LEN      = 1200
RNG_SEED       = 2

rng = np.random.RandomState(RNG_SEED)
results = []  # (N, gap_mean, gap_std, all_positive, sign_threshold, decomp_label)

for N in range(4, 51):
    if isprime(N):
        continue
    if N % 2 == 0:
        C_add, b_add = 2, N // 2
    else:
        p = min(factorint(N))
        C_add, b_add = p, N // p
    if b_add < 2:
        continue

    gaps = np.empty(SEEDS_PER_BASE)
    for i in range(SEEDS_PER_BASE):
        seed = mp.mpf(rng.uniform(0.05, 0.95))
        gaps[i] = lyap_mul(N, seed, ORBIT_LEN) - lyap_add(C_add, b_add, seed, ORBIT_LEN)

    # analytic sign-threshold for <y>:
    #   gap_sign( <y> ) > 0  iff  <y> > 1 - ln(ln N / ln b_add) / ln C_add
    if C_add > 1:
        thr = 1.0 - float(mp.log(mp.log(N) / mp.log(b_add)) / mp.log(C_add))
    else:
        thr = float("nan")

    results.append((
        N,
        float(gaps.mean()),
        float(gaps.std()),
        bool(np.all(gaps > 0)),
        thr,
        f"{C_add}·{b_add}^y",
    ))

# ---------- print table ----------

print(f"{'N':>3}  {'gap_mean':>9}  {'std':>6}  {'allpos':>6}  {'thr':>6}   decomp")
for N, gm, gs, ap, thr, lbl in results:
    print(f"{N:3d}  {gm:+9.4f}  {gs:6.3f}  {str(ap):>6}  {thr:6.3f}   {lbl}")

robust = [r[0] for r in results if r[3]]
print(f"\nRobust-split bases (mul wins on every seed): {robust}")


# ---------- plot ----------

Ns      = np.array([r[0] for r in results])
gaps    = np.array([r[1] for r in results])
stds    = np.array([r[2] for r in results])
allpos  = [r[3] for r in results]
thrs    = np.array([r[4] for r in results])
colors  = ["crimson" if a else "steelblue" for a in allpos]

ORBIT_MEAN_APPROX = 0.45  # empirical orbit mean across these maps; see code comment in lemma

fig, ax = plt.subplots(1, 2, figsize=(16, 6))

# Panel 1: gap vs N
ax[0].axhline(0, color="k", lw=0.8)
ax[0].errorbar(Ns, gaps, yerr=stds, fmt="none", ecolor="gray", alpha=0.5, zorder=1)
ax[0].scatter(Ns, gaps, c=colors, s=70, zorder=3, edgecolor="k", linewidth=0.5)
for r in results:
    if r[3]:
        ax[0].annotate(str(r[0]), (r[0], r[1]),
                       textcoords="offset points", xytext=(0, 9),
                       ha="center", fontweight="bold")
ax[0].set_title("Lyapunov deflection (mul − add) across composite bases\n"
                "RED = robust split (mul wins on all seeds);  BLUE = no robust split")
ax[0].set_xlabel("composite base N")
ax[0].set_ylabel("gap  λ_mul − λ_add")

# Panel 2: sign-threshold vs orbit mean
ax[1].scatter(Ns, thrs, c=colors, s=70, edgecolor="k", linewidth=0.5, zorder=3)
ax[1].axhline(ORBIT_MEAN_APPROX, color="green", ls="--", lw=1.5,
              label=f"≈ orbit mean ⟨y⟩ ({ORBIT_MEAN_APPROX})")
ax[1].fill_between([2, 52], 0,                 ORBIT_MEAN_APPROX, color="crimson",  alpha=0.08)
ax[1].fill_between([2, 52], ORBIT_MEAN_APPROX, 0.8,               color="steelblue", alpha=0.08)
ax[1].set_title("WHY: sign-threshold  1 − ln(ln N / ln b_add) / ln C   vs orbit mean\n"
                "below green line → split robust (red);  above → killed (blue)")
ax[1].set_xlabel("composite base N")
ax[1].set_ylabel("sign threshold for ⟨y⟩")
ax[1].set_xlim(2, 52)
ax[1].legend()
for r in results:
    if r[3]:
        ax[1].annotate(str(r[0]), (r[0], r[4]),
                       textcoords="offset points", xytext=(0, 9),
                       ha="center", fontweight="bold")

plt.tight_layout()
plt.savefig("regime_boundary.png", dpi=120, bbox_inches="tight")
print("\nwrote regime_boundary.png")
