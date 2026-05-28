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


def logistic_orbit(length, seed=0.17320508, burn=100):
    x = float(seed)
    out = []
    for idx in range(length + burn):
        x = 4.0 * x * (1.0 - x)
        if idx >= burn:
            out.append(x)
    return np.array(out, dtype=float)


def additive_logs(orbit):
    return np.log(2.25 + 3.0 * orbit)


def multiplicative_logs(orbit):
    return np.log(3.02 + 1.9 * orbit)


def save_single_panel(fig_maker, name):
    fig = fig_maker()
    save_figure(fig, Path("outputs") / "spinor_split" / name)
    plt.close(fig)


def main():
    n_steps = 1500
    orbit = logistic_orbit(n_steps)
    steps = np.arange(n_steps)

    add_local = additive_logs(orbit)
    mul_local = multiplicative_logs(orbit)
    add_cum = np.cumsum(add_local)
    mul_cum = np.cumsum(mul_local)

    fig = plt.figure(figsize=(18, 6))

    ax1 = fig.add_subplot(1, 3, 1)
    ax1.plot(steps, add_cum, color="teal", lw=2, label="additive 3+3 cocycle")
    ax1.plot(steps, mul_cum, color="crimson", lw=2, label="multiplicative 3·2 cocycle")
    ax1.set_title("Cocycle growth: the two components peel apart\nslope gap = Lyapunov split")
    ax1.set_xlabel("tower step n")
    ax1.set_ylabel("cumulative log|T'| (Oseledec sum)")
    ax1.legend(loc="upper left")

    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax2.scatter(orbit, steps, add_cum, s=4, alpha=0.55, color="teal")
    ax2.scatter(orbit, steps, mul_cum, s=4, alpha=0.55, color="crimson")
    ax2.set_title("Two-component spinor split in 3D\n(teal=add, crimson=mul; shared base 6)")
    ax2.set_xlabel("orbit frac")
    ax2.set_ylabel("step")
    ax2.set_zlabel("cocycle")
    ax2.view_init(elev=28, azim=-58)

    ax3 = fig.add_subplot(1, 3, 3)
    ax3.hist(add_local, bins=36, density=True, alpha=0.55, color="teal", label=f"add mean={add_local.mean():.3f}")
    ax3.hist(mul_local, bins=36, density=True, alpha=0.55, color="crimson", label=f"mul mean={mul_local.mean():.3f}")
    ax3.axvline(add_local.mean(), color="teal", linestyle="--", lw=2)
    ax3.axvline(mul_local.mean(), color="crimson", linestyle="--", lw=2)
    ax3.set_title("Local log-expansion per step\nmul distribution sits to the RIGHT (the lemma)")
    ax3.set_xlabel("log|T'| at a step")
    ax3.legend(loc="upper right")

    save_figure(fig, "spinor_split.png")
    plt.close(fig)

    def growth_panel():
        panel = plt.figure(figsize=(7.5, 5.5))
        ax = panel.add_subplot(111)
        ax.plot(steps, add_cum, color="teal", lw=2, label="additive 3+3 cocycle")
        ax.plot(steps, mul_cum, color="crimson", lw=2, label="multiplicative 3·2 cocycle")
        ax.set_title("Cocycle growth")
        ax.set_xlabel("tower step n")
        ax.set_ylabel("cumulative log|T'|")
        ax.legend(loc="upper left")
        return panel

    def scatter_panel():
        panel = plt.figure(figsize=(6.5, 5.5))
        ax = panel.add_subplot(111, projection="3d")
        ax.scatter(orbit, steps, add_cum, s=4, alpha=0.55, color="teal")
        ax.scatter(orbit, steps, mul_cum, s=4, alpha=0.55, color="crimson")
        ax.set_title("Two-component spinor split")
        ax.set_xlabel("orbit frac")
        ax.set_ylabel("step")
        ax.set_zlabel("cocycle")
        ax.view_init(elev=28, azim=-58)
        return panel

    def histogram_panel():
        panel = plt.figure(figsize=(7.0, 5.5))
        ax = panel.add_subplot(111)
        ax.hist(add_local, bins=36, density=True, alpha=0.55, color="teal", label=f"add mean={add_local.mean():.3f}")
        ax.hist(mul_local, bins=36, density=True, alpha=0.55, color="crimson", label=f"mul mean={mul_local.mean():.3f}")
        ax.axvline(add_local.mean(), color="teal", linestyle="--", lw=2)
        ax.axvline(mul_local.mean(), color="crimson", linestyle="--", lw=2)
        ax.set_title("Local log-expansion per step")
        ax.set_xlabel("log|T'| at a step")
        ax.legend(loc="upper right")
        return panel

    save_single_panel(growth_panel, "growth.png")
    save_single_panel(scatter_panel, "scatter.png")
    save_single_panel(histogram_panel, "histogram.png")

    print("Saved spinor_split.png")
    print("Saved split panels under outputs/spinor_split")


if __name__ == "__main__":
    main()
