from pathlib import Path


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def save_figure(fig, path, dpi=200):
    path = Path(path)
    ensure_dir(path.parent)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path
