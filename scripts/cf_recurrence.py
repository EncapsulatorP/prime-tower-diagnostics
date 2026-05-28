import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str((Path("outputs") / ".mplconfig").resolve()),
)

import matplotlib
import mpmath as mp
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scripts.plot_io import save_figure
except ModuleNotFoundError:
    from plot_io import save_figure


mp.mp.dps = 220


def continued_fraction_terms(value, n_terms=300):
    x = mp.mpf(value)
    terms = []
    for _ in range(n_terms):
        a = int(mp.floor(x))
        terms.append(a)
        frac = x - a
        if abs(frac) < mp.mpf("1e-80"):
            break
        x = 1 / frac
    return np.array(terms, dtype=int)


def recurrence_matrix(terms):
    return terms[:, None] == terms[None, :]


def draw_panel(ax, label, terms):
    matrix = recurrence_matrix(terms)
    ax.imshow(
        matrix,
        cmap="gray_r",
        origin="lower",
        interpolation="nearest",
        vmin=0,
        vmax=1,
    )
    ax.set_title(f"{label}: CF recurrence plot", fontsize=20)


def save_split_panel(label, terms, slug):
    fig, ax = plt.subplots(figsize=(6, 6))
    draw_panel(ax, label, terms)
    save_figure(fig, Path("outputs") / "recurrence" / f"{slug}.png")
    plt.close(fig)


def main():
    panels = [
        (r"sqrt2", "sqrt2", mp.sqrt(2)),
        (r"e", "e", mp.e),
        (r"pi", "pi", mp.pi),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, (label, slug, value) in zip(axes, panels):
        terms = continued_fraction_terms(value)
        draw_panel(ax, label, terms)
        save_split_panel(label, terms, slug)

    save_figure(fig, "recurrence.png")
    plt.close(fig)

    print("Saved recurrence.png")
    print("Saved split panels under outputs/recurrence")


if __name__ == "__main__":
    main()
