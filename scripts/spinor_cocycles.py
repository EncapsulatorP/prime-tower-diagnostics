import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str((Path("outputs") / ".mplconfig").resolve()),
)

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scripts.plot_io import save_figure
except ModuleNotFoundError:
    from plot_io import save_figure


# ============================================================
# COCYCLE + SPINOR-LIKE FRACTAL SPLASH LAB
# Newton vs Collatz-like vs Pisot-like beta fold
# ============================================================

N = 450
xmin, xmax = -2.0, 2.0
ymin, ymax = -2.0, 2.0

x = np.linspace(xmin, xmax, N)
y = np.linspace(ymin, ymax, N)
X, Y = np.meshgrid(x, y)
Z0 = X + 1j * Y
OUTDIR = Path("outputs") / "spinor_cocycles"


def normalize(A):
    A = np.asarray(A, dtype=float)
    lo, hi = np.nanpercentile(A, 2), np.nanpercentile(A, 98)
    if hi == lo:
        return np.zeros_like(A)
    return np.clip((A - lo) / (hi - lo), 0, 1)


def safe_log_abs(w, eps=1e-12):
    return np.log(np.maximum(np.abs(w), eps))


def export_current_figure(name, show=False):
    fig = plt.gcf()
    save_figure(fig, OUTDIR / name)
    if show:
        plt.show()
    plt.close(fig)


# ============================================================
# 1. NEWTON MAP AND DERIVATIVE
#
# T(z) = z - (z^3 - 1)/(3z^2)
#      = (2z^3 + 1)/(3z^2)
#
# derivative:
# T'(z) = 2/3 - 2/(3z^3)
# ============================================================

def newton_step(z):
    dz = 3 * z**2
    safe = np.abs(dz) > 1e-12
    out = z.copy()
    out[safe] = z[safe] - (z[safe]**3 - 1) / dz[safe]
    return out


def newton_derivative(z):
    safe = np.abs(z) > 1e-8
    d = np.zeros_like(z, dtype=complex)
    d[safe] = (2 / 3) - (2 / (3 * z[safe]**3))
    return d


# ============================================================
# 2. COMPLEX COLLATZ INTERPOLATION AND DERIVATIVE
#
# T(z) = 1/2[(z/2 + 3z+1) + cos(pi z)(z/2 - (3z+1))]
#
# This is only an analytic interpolation of Collatz, not the
# integer Collatz map itself.
# ============================================================

def collatz_step(z):
    even_branch = z / 2
    odd_branch = 3 * z + 1
    with np.errstate(over="ignore", invalid="ignore"):
        return 0.5 * (even_branch + odd_branch) + 0.5 * np.cos(np.pi * z) * (even_branch - odd_branch)


def collatz_derivative(z):
    even_branch = z / 2
    odd_branch = 3 * z + 1

    de = 1 / 2
    do = 3

    A = even_branch + odd_branch
    B = even_branch - odd_branch

    dA = de + do
    dB = de - do

    with np.errstate(over="ignore", invalid="ignore"):
        return 0.5 * dA + 0.5 * (
            -np.pi * np.sin(np.pi * z) * B + np.cos(np.pi * z) * dB
        )


# ============================================================
# 3. PISOT-LIKE BETA FOLD AND DERIVATIVE PROXY
#
# This is not a rigorous Rauzy/Pisot fractal.
# It is a beta-fold compression proxy.
#
# Fold discontinuities make the derivative singular at cut lines.
# Away from cuts, derivative is approximately multiplication by beta.
# ============================================================

PHI = (1 + np.sqrt(5)) / 2


def frac_complex(z):
    return (z.real - np.floor(z.real)) + 1j * (z.imag - np.floor(z.imag))


def pisot_step(z, beta=PHI):
    return frac_complex(beta * z) - (0.5 + 0.5j)


def pisot_derivative(z, beta=PHI):
    return beta + 0j * z


# ============================================================
# 4. FINITE-TIME COCYCLE / LYAPUNOV FIELD
#
# For complex analytic maps, the real 2x2 Jacobian norm is |T'(z)|.
# Therefore the cocycle log norm is sum log |T'(z_k)|.
# ============================================================

def lyapunov_field(step_fn, deriv_fn, Z, max_iter=40, clip=1e6):
    z = Z.copy()
    log_norm_sum = np.zeros(Z.shape, dtype=float)
    angle_sum = np.zeros(Z.shape, dtype=float)
    alive = np.ones(Z.shape, dtype=bool)

    for k in range(max_iter):
        d = deriv_fn(z)

        log_norm_sum[alive] += safe_log_abs(d[alive])
        angle_sum[alive] += np.angle(d[alive])

        z_next = step_fn(z)

        # crude numerical safety
        bad = ~np.isfinite(z_next.real) | ~np.isfinite(z_next.imag) | (np.abs(z_next) > clip)
        alive &= ~bad

        z[alive] = z_next[alive]

    lam = log_norm_sum / max_iter

    # spinor-like half-angle channel
    spinor_phase = np.exp(0.5j * angle_sum)

    return lam, spinor_phase, z


# ============================================================
# 5. RUN COCYCLES
# ============================================================

ITER = 45

lam_newton, spin_newton, zN = lyapunov_field(
    newton_step, newton_derivative, Z0, max_iter=ITER
)

lam_collatz, spin_collatz, zC = lyapunov_field(
    collatz_step, collatz_derivative, Z0, max_iter=ITER
)

lam_pisot, spin_pisot, zP = lyapunov_field(
    pisot_step, pisot_derivative, Z0, max_iter=ITER
)


# ============================================================
# 6. COCYCLE SPLASH
#
# Red   = Newton Lyapunov
# Green = Collatz-like Lyapunov
# Blue  = Pisot-like Lyapunov
# ============================================================

R = normalize(lam_newton)
G = normalize(lam_collatz)
B = normalize(lam_pisot)

RGB = np.dstack([R, G, B])

plt.figure(figsize=(9, 9))
plt.imshow(RGB, extent=[xmin, xmax, ymin, ymax], origin="lower")
plt.title("Cocycle Splash: Newton=Red, Collatz-like=Green, Pisot-like=Blue")
plt.xlabel("Re(z)")
plt.ylabel("Im(z)")
plt.tight_layout()
export_current_figure("cocycle_splash.png")


# ============================================================
# 7. COCYCLE COLLISION MAP
#
# Bright means the three finite-time Lyapunov signatures collapse
# to similar values under this compression.
# ============================================================

stack_lam = np.dstack([R, G, B])
cocycle_collision = 1.0 - normalize(np.var(stack_lam, axis=2))

plt.figure(figsize=(8, 8))
plt.imshow(cocycle_collision, extent=[xmin, xmax, ymin, ymax], origin="lower")
plt.title("Cocycle Collision Map: bright = similar compressed Lyapunov signature")
plt.xlabel("Re(z)")
plt.ylabel("Im(z)")
plt.tight_layout()
export_current_figure("collision.png")


# ============================================================
# 8. SPINOR-LIKE HALF-ANGLE CHANNELS
#
# Use half-angle phase:
# psi = exp(i theta / 2)
#
# This is not a true spinor bundle.
# It is a double-cover/orientation diagnostic.
# ============================================================

def phase_to_rgb(psi):
    angle = np.angle(psi)
    hue = (angle + np.pi) / (2 * np.pi)

    # simple HSV-like conversion, no external libs
    r = 0.5 + 0.5 * np.cos(2 * np.pi * hue)
    g = 0.5 + 0.5 * np.cos(2 * np.pi * hue - 2 * np.pi / 3)
    b = 0.5 + 0.5 * np.cos(2 * np.pi * hue - 4 * np.pi / 3)

    return np.dstack([r, g, b])


plt.figure(figsize=(8, 8))
plt.imshow(phase_to_rgb(spin_newton), extent=[xmin, xmax, ymin, ymax], origin="lower")
plt.title("Spinor-like half-angle phase: Newton")
plt.xticks([])
plt.yticks([])
plt.tight_layout()
export_current_figure("phase_newton.png")

plt.figure(figsize=(8, 8))
plt.imshow(phase_to_rgb(spin_collatz), extent=[xmin, xmax, ymin, ymax], origin="lower")
plt.title("Spinor-like half-angle phase: Collatz-like")
plt.xticks([])
plt.yticks([])
plt.tight_layout()
export_current_figure("phase_collatz.png")

plt.figure(figsize=(8, 8))
plt.imshow(phase_to_rgb(spin_pisot), extent=[xmin, xmax, ymin, ymax], origin="lower")
plt.title("Spinor-like half-angle phase: Pisot-like beta fold")
plt.xticks([])
plt.yticks([])
plt.tight_layout()
export_current_figure("phase_pisot.png")


# ============================================================
# 9. SPINOR MISMATCH / RESIDUE MAP
#
# If two systems compress similarly but their half-angle phases
# disagree, this marks a hidden orientation residue.
# ============================================================

def spinor_mismatch(a, b):
    # phase distance on unit circle
    return np.abs(1 - a * np.conj(b))


mismatch_NC = normalize(spinor_mismatch(spin_newton, spin_collatz))
mismatch_NP = normalize(spinor_mismatch(spin_newton, spin_pisot))
mismatch_CP = normalize(spinor_mismatch(spin_collatz, spin_pisot))

mismatch_rgb = np.dstack([mismatch_NC, mismatch_NP, mismatch_CP])

plt.figure(figsize=(9, 9))
plt.imshow(mismatch_rgb, extent=[xmin, xmax, ymin, ymax], origin="lower")
plt.title("Spinor-like Residue: R=N/C, G=N/P, B=C/P")
plt.xlabel("Re(z)")
plt.ylabel("Im(z)")
plt.tight_layout()
export_current_figure("spinor_residue.png")


# ============================================================
# 10. SUMMARY NUMBERS
# ============================================================

def summarize(name, lam):
    finite = lam[np.isfinite(lam)]
    print(f"{name}")
    print(f"  mean lambda: {np.mean(finite): .6f}")
    print(f"  std lambda:  {np.std(finite): .6f}")
    print(f"  min lambda:  {np.min(finite): .6f}")
    print(f"  max lambda:  {np.max(finite): .6f}")
    print()


summarize("Newton", lam_newton)
summarize("Collatz-like", lam_collatz)
summarize("Pisot-like", lam_pisot)

print("Mean cocycle collision:", np.nanmean(cocycle_collision))
print("Mean spinor mismatch Newton/Collatz:", np.nanmean(mismatch_NC))
print("Mean spinor mismatch Newton/Pisot:  ", np.nanmean(mismatch_NP))
print("Mean spinor mismatch Collatz/Pisot: ", np.nanmean(mismatch_CP))
print("Saved split figures under outputs/spinor_cocycles")
