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


def plot_variants(table, value: str, model_rules, masks, primary, xlabel: str, title: str, path) -> None:
    """One row per model and rule, one dot per mask level, with the primary specification ringed."""
    fig, ax = plt.subplots(figsize=(7.2, 0.62 * len(model_rules) + 1.7), constrained_layout=True)
    rows = {mr: i for i, mr in enumerate(reversed(model_rules))}
    by_combo = table.set_index(["model", "rule", "mask"])
    for mr, row in rows.items():
        values = [by_combo.loc[(*mr, mask), value] for mask in masks]
        ax.hlines(row, min(values), max(values), color=AXIS, linewidth=2, zorder=1)
    for mask in masks:
        ax.scatter([by_combo.loc[(*mr, mask), value] for mr in rows], list(rows.values()), s=64, color=MASK_COLORS[mask],
                   edgecolor=SURFACE, linewidth=1.5, zorder=3, label=MASK_NAMES[mask])
    pr = by_combo.loc[primary, value]
    ax.scatter([pr], [rows[primary[:2]]], s=210, facecolor="none", edgecolor=INK, linewidth=1.4, zorder=4, label="primary specification")
    ax.annotate(f"primary, r = {pr:.3f}", xy=(pr, rows[primary[:2]]), xytext=(0, 13), textcoords="offset points",
                ha="center", color=INK, path_effects=HALO)
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_yticks(list(rows.values()), [combo_label(*mr) for mr in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    fig.legend(loc="outside lower center", ncol=2, handletextpad=0.3, columnspacing=1.6)
    save(fig, path)
