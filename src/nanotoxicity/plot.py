"""Observed-versus-predicted figure for the leave-one-out fits."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt

from nanotoxicity.model import LoocvResult
from nanotoxicity.qdb import OxideTable


def save_figure(
    tables: list[OxideTable],
    fits: list[LoocvResult],
    path: Path,
) -> None:
    if len(tables) != len(fits):
        raise ValueError("Each table needs one fit.")
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(1, len(tables), figsize=(6.3 * len(tables), 6.4))
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.26, top=0.86, wspace=0.32)
    if len(tables) == 1:
        axes = [axes]
    for axis, table, fit in zip(axes, tables, fits):
        _draw_panel(axis, table, fit)
    fig.suptitle(
        "Leave-one-out linear predictions from the deposited descriptors",
        fontsize=12,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _draw_panel(axis, table: OxideTable, fit: LoocvResult) -> None:
    observed = table.y
    predicted = fit.y_pred
    lower = float(min(observed.min(), predicted.min()))
    upper = float(max(observed.max(), predicted.max()))
    pad = 0.1 * (upper - lower)
    lower -= pad
    upper += pad
    axis.plot(
        [lower, upper],
        [lower, upper],
        linestyle="--",
        color="#4d4d4d",
        linewidth=1,
        zorder=1,
    )
    axis.scatter(observed, predicted, s=32, color="#0072B2", zorder=3)
    axis.set_xlim(lower, upper)
    axis.set_ylim(lower, upper)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(f"Observed {table.property_id}")
    axis.set_ylabel(f"Leave-one-out {table.property_id}")
    axis.set_title(f"{table.system}\n{table.citation}, n = {len(observed)}")
    note = (
        f"RMSE  model {fit.model.rmse:.3f}\n"
        f"RMSE  mean  {fit.baseline.rmse:.3f}\n"
        f"MAE   model {fit.model.mae:.3f}\n"
        f"MAE   mean  {fit.baseline.mae:.3f}"
    )
    axis.text(
        0.04,
        0.96,
        note,
        transform=axis.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        family="monospace",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#cccccc"},
        zorder=4,
    )
    halo = [path_effects.withStroke(linewidth=2.4, foreground="white")]
    for index, (x_value, y_value) in enumerate(zip(observed, predicted), start=1):
        # Neighbours in this list are neighbours on the x axis, so alternate sides.
        xytext = (5, 4) if index % 2 else (5, -9)
        axis.annotate(
            str(index),
            (x_value, y_value),
            textcoords="offset points",
            xytext=xytext,
            fontsize=7,
            color="#222222",
            path_effects=halo,
            zorder=5,
        )
    axis.text(
        0.0,
        -0.22,
        _key(table.formulas),
        transform=axis.transAxes,
        va="top",
        ha="left",
        fontsize=7.5,
        family="monospace",
        color="#222222",
    )


def _key(formulas: tuple[str, ...]) -> str:
    pairs = [f"{index} {formula}" for index, formula in enumerate(formulas, start=1)]
    row_width = 6
    rows = [
        "   ".join(pairs[start : start + row_width]) for start in range(0, len(pairs), row_width)
    ]
    return "\n".join(rows)
