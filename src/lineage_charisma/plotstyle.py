from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# one palette and one set of rcParams for every figure in the project
SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
MASK_COLORS = {"none": BLUE, "taxonomy": ORANGE, "strict": AQUA}
MASK_NAMES = {"none": "no masking", "taxonomy": "taxonomy masked", "strict": "strict (taxonomy and common names masked)"}
BLUE_RAMP = LinearSegmentedColormap.from_list("blue_seq", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
ORANGE_RAMP = LinearSegmentedColormap.from_list("orange_seq", ["#fde3d6", "#f7b79b", "#f08a5f", "#eb6834", "#c04f22", "#8f3815", "#5e230b"])
HALO = [pe.withStroke(linewidth=2.5, foreground=SURFACE)]
MEDIAN_LINE = dict(color=INK, linewidth=1.6, marker="o", markersize=5, markeredgecolor=SURFACE, markeredgewidth=1.2)


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "text.color": INK,
        "axes.labelcolor": INK_2,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "legend.frameon": False,
    })


def save(fig, path) -> None:
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def combo_label(model: str, rule: str, mask: str | None = None) -> str:
    return f"{model.split('/')[-1]} / {rule}" + (f" / {mask}" if mask else "")
