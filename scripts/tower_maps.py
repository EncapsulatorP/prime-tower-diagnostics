import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str((Path("outputs") / ".mplconfig").resolve()),
)

import math

import matplotlib
import mpmath as mp
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scripts.plot_io import save_figure
except ModuleNotFoundError:
    from plot_io import save_figure


mp.mp.dps = 80


def liouville_constant(kmax=6):
    total = mp.mpf("0")
    for k in range(1, kmax + 1):
        total += mp.power(10, -math.factorial(k))
    return total


def fractional_part(x):
    return x - mp.floor(x)


def delay_embed(values):
    values = np.asarray(values, dtype=float)
    return np.column_stack([values[:-2], values[1:-1], values[2:]])


def approximate_d2(points, max_points=420, seed=0):
    points = np.asarray(points, dtype=float)
    n = len(points)
    if n < 6:
        return 0.0

    if n > max_points:
        rng = np.random.default_rng(seed)
        keep = np.sort(rng.choice(n, size=max_points, replace=False))
        points = points[keep]

    diffs = points[:, None, :] - points[None, :, :]
    dists = np.sqrt(np.sum(diffs * diffs, axis=2))
    upper = dists[np.triu_indices(len(points), k=1)]
    upper = upper[upper > 0]
    if len(upper) < 20:
        return 0.0

    lo = np.quantile(upper, 0.10)
    hi = np.quantile(upper, 0.80)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return 0.0

    radii = np.exp(np.linspace(np.log(lo), np.log(hi), 9))
    corr = np.array([(upper <= r).mean() for r in radii])
    mask = (corr > 0) & (corr < 1)
    if mask.sum() < 4:
        return 0.0

    coeffs = np.polyfit(np.log(radii[mask]), np.log(corr[mask]), 1)
    return float(np.clip(coeffs[0], 0.0, 3.0))


def positional_orbit(base, y0=0.37, burn=120, steps=1500):
    y = float(y0)
    out = []
    for i in range(burn + steps + 2):
        y = (float(base) ** y) % 1.0
        if i >= burn:
            out.append(y)
    return np.array(out, dtype=float)


def gauss_tower_cloud(constant, nseeds=700, burn=8):
    constant = mp.mpf(constant)
    points = []
    for k in range(1, nseeds + 1):
        y = fractional_part(k * constant)
        if y <= 0:
            y = mp.mpf("0.3141592653589793")

        orbit = []
        for _ in range(burn + 3):
            y = fractional_part(mp.power(1 / y, 1 / y))
            if y <= 0:
                y = mp.mpf("0.2718281828459045")
            orbit.append(float(y))

        points.append(orbit[-3:])

    return np.asarray(points, dtype=float)


def plastic_constant():
    roots = np.roots([1.0, 0.0, -1.0, -1.0])
    return float(np.max(roots.real[np.isclose(roots.imag, 0.0)]))


def tribonacci_constant():
    roots = np.roots([1.0, -1.0, -1.0, -1.0])
    return float(np.max(roots.real[np.isclose(roots.imag, 0.0)]))


def quartic_pisot_constant():
    roots = np.roots([1.0, -1.0, 0.0, 0.0, -1.0])
    return float(np.max(roots.real[np.isclose(roots.imag, 0.0)]))


def make_scatter_panel(ax, points, label, cmap="plasma"):
    colors = np.linspace(0.0, 1.0, len(points))
    ax.scatter(
        points[:, 0],
        points[:, 1],
        points[:, 2],
        c=colors,
        cmap=cmap,
        s=6,
        alpha=0.72,
        linewidths=0,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=34, azim=-62)
    ax.set_title(f"{label}\nD2~{approximate_d2(points):.2f}", pad=8)


def render_family(title, panels, point_builder, combined_path, split_dir):
    fig = plt.figure(figsize=(19, 4.4))
    fig.suptitle(title, fontsize=20, y=0.98)

    for idx, panel in enumerate(panels, start=1):
        ax = fig.add_subplot(1, len(panels), idx, projection="3d")
        points = point_builder(panel["value"])
        if np.asarray(points).ndim == 1:
            points = delay_embed(points)
        make_scatter_panel(ax, points, panel["label"])

        solo = plt.figure(figsize=(4.2, 4.2))
        solo_ax = solo.add_subplot(111, projection="3d")
        make_scatter_panel(solo_ax, points, panel["label"])
        save_figure(solo, Path(split_dir) / f"{panel['slug']}.png")
        plt.close(solo)

    save_figure(fig, combined_path)
    plt.close(fig)


def main():
    liouville = float(liouville_constant())
    phi = (1.0 + math.sqrt(5.0)) / 2.0

    pure_positional = [
        {"label": r"$\sqrt{2}$", "slug": "sqrt2", "value": math.sqrt(2.0)},
        {"label": r"$\phi$", "slug": "phi", "value": phi},
        {"label": r"$e$", "slug": "e", "value": math.e},
        {"label": r"$\pi$", "slug": "pi", "value": math.pi},
        {
            "label": r"Liouville $\Sigma 10^{-k!}$",
            "slug": "liouville",
            "value": liouville,
        },
    ]

    gauss_panels = [
        {"label": r"$\sqrt{2}$", "slug": "sqrt2", "value": math.sqrt(2.0)},
        {"label": r"$\phi$", "slug": "phi", "value": phi},
        {"label": r"$e$", "slug": "e", "value": math.e},
        {"label": r"$\pi$", "slug": "pi", "value": math.pi},
        {
            "label": r"Liouville $\Sigma 10^{-k!}$",
            "slug": "liouville",
            "value": liouville,
        },
    ]

    pisot_panels = [
        {
            "label": r"$\log \rho$ (plastic)",
            "slug": "plastic",
            "value": math.log(plastic_constant()),
        },
        {
            "label": r"$\log(\mathrm{tribonacci})$",
            "slug": "tribonacci",
            "value": math.log(tribonacci_constant()),
        },
        {
            "label": r"$\log \phi$",
            "slug": "log_phi",
            "value": math.log(phi),
        },
        {
            "label": r"$\log(\mathrm{quartic\ Pisot})$",
            "slug": "quartic_pisot",
            "value": math.log(quartic_pisot_constant()),
        },
        {
            "label": r"$\pi$ (reference)",
            "slug": "pi_reference",
            "value": math.pi,
        },
    ]

    render_family(
        title=r"Gauss-tower map  T(y)=frac((1/y)^(1/y))",
        panels=gauss_panels,
        point_builder=gauss_tower_cloud,
        combined_path="gt_tower.png",
        split_dir=Path("outputs") / "gt_tower",
    )
    render_family(
        title=r"Pure positional tower  T(y)=frac(b^y), b=the number",
        panels=pure_positional,
        point_builder=positional_orbit,
        combined_path="pt_tower.png",
        split_dir=Path("outputs") / "pt_tower",
    )
    render_family(
        title=r"Pure positional tower  T(y)=frac(b^y)  --  logs of Pisot numbers vs $\pi$",
        panels=pisot_panels,
        point_builder=positional_orbit,
        combined_path="pisot_tower.png",
        split_dir=Path("outputs") / "pisot_tower",
    )

    print("Saved gt_tower.png, pt_tower.png, pisot_tower.png")
    print("Saved split panels under outputs/gt_tower, outputs/pt_tower, outputs/pisot_tower")


if __name__ == "__main__":
    main()
