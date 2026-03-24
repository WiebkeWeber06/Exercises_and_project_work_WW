# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 12:39:17 2026

@author: wiewe372
"""

import math
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd


REQUIRED_RAW_COLUMNS = ["gene", "genotype", "condition", "time", "replicate", "value"]
REQUIRED_SUMMARY_COLUMNS = ["gene", "genotype", "condition", "time", "mean", "sd"]


def validate_columns(df: pd.DataFrame, required: Sequence[str]) -> None:
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def summarize_replicates(
    df: pd.DataFrame,
    group_cols: Optional[Sequence[str]] = None,
    value_col: str = "value",
) -> pd.DataFrame:
    """
    Summarize replicate-level qPCR data into mean, standard deviation and n.

    Parameters
    ----------
    df : pandas.DataFrame
        Long-format dataframe containing replicate-level measurements.
    group_cols : sequence of str, optional
        Columns used to group replicates. If None, defaults to
        ['gene', 'genotype', 'condition', 'time'].
    value_col : str
        Name of the measurement column.

    Returns
    -------
    pandas.DataFrame
        Summary dataframe with columns:
        group columns + ['mean', 'sd', 'n'].
    """
    validate_columns(df, REQUIRED_RAW_COLUMNS)

    if group_cols is None:
        group_cols = ["gene", "genotype", "condition", "time"]

    summary = (
        df.groupby(list(group_cols), dropna=False)[value_col]
        .agg(mean="mean", sd="std", n="count")
        .reset_index()
    )

    summary["sd"] = summary["sd"].fillna(0.0)
    return summary


def _ordered_unique(values: pd.Series, preferred_order: Optional[Sequence[str]] = None) -> list:
    unique_vals = list(pd.unique(values))
    if preferred_order is None:
        return unique_vals
    ordered = [x for x in preferred_order if x in unique_vals]
    remaining = [x for x in unique_vals if x not in ordered]
    return ordered + remaining


def plot_expression(
    df_summary: pd.DataFrame,
    *,
    facet_rows: Optional[str] = "gene",
    facet_cols: Optional[str] = "condition",
    line_group: str = "genotype",
    x: str = "time",
    y: str = "mean",
    yerr: Optional[str] = "sd",
    title: str = "qPCR expression plot",
    ylabel: str = "Relative expression",
    xlabel: Optional[str] = None,
    row_order: Optional[Sequence[str]] = None,
    col_order: Optional[Sequence[str]] = None,
    x_order: Optional[Sequence[str]] = None,
    line_order: Optional[Sequence[str]] = None,
    colors: Optional[dict] = None,
    markers: Optional[dict] = None,
    linestyles: Optional[dict] = None,
    sharey: str | bool = "row",
    figsize: Optional[tuple] = None,
    legend: bool = True,
    legend_title: Optional[str] = None,
    legend_anchor_x: float = 0.845,
    right_margin: float = 0.82,
    top_margin: float = 0.88,
    bottom_margin: float = 0.14,
    left_margin: float = 0.08,
    wspace: float = 0.28,
    hspace: float = 0.28,
    rotate_xticks: int = 45,
    grid: bool = True,
    yscale: str = "linear",
    ylim: Optional[tuple] = None,
    xlim: Optional[tuple] = None,
    show_points: bool = True,
    capsize: float = 4,
    linewidth: float = 2,
    markersize: float = 6,
) -> tuple:
    """
    Plot summarized qPCR-like data from a long-format summary dataframe.

    Parameters
    ----------
    df_summary : pandas.DataFrame
        Summary dataframe with columns such as gene, genotype, condition, time,
        mean, sd, n.
    facet_rows : str or None
        Column used to define subplot rows.
    facet_cols : str or None
        Column used to define subplot columns.
    line_group : str
        Column used to distinguish plotted lines within each panel.
    x : str
        Column used for x-axis categories.
    y : str
        Column containing plotted values.
    yerr : str or None
        Column containing error bar values.
    title : str
        Figure title.
    ylabel : str
        Y-axis label.
    xlabel : str or None
        X-axis label. If None, uses `x`.
    row_order, col_order, x_order, line_order : sequence of str, optional
        Explicit plotting order.
    colors, markers, linestyles : dict, optional
        Style dictionaries keyed by values of `line_group`.
    sharey : {'row', 'col', 'all', True, False}
        Y-axis sharing behavior.
    yscale : str
        Y-axis scale, e.g. 'linear' or 'log'.
    ylim, xlim : tuple, optional
        Axis limits.

    Returns
    -------
    fig, axes
        Matplotlib figure and axes.
    """
    validate_columns(df_summary, REQUIRED_SUMMARY_COLUMNS)

    plot_df = df_summary.copy()

    row_levels = [None] if facet_rows is None else _ordered_unique(plot_df[facet_rows], row_order)
    col_levels = [None] if facet_cols is None else _ordered_unique(plot_df[facet_cols], col_order)
    x_levels = _ordered_unique(plot_df[x], x_order)
    line_levels = _ordered_unique(plot_df[line_group], line_order)

    nrows = len(row_levels)
    ncols = len(col_levels)

    if figsize is None:
        figsize = (5.5 * ncols + 2.5, 3.8 * nrows + 1.2)

    sharey_map = {"all": True, "row": "row", "col": "col", True: True, False: False}
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize,
        squeeze=False,
        sharey=sharey_map.get(sharey, sharey),
    )

    x_pos = {label: i for i, label in enumerate(x_levels)}

    for i, row_val in enumerate(row_levels):
        for j, col_val in enumerate(col_levels):
            ax = axes[i, j]
            panel_df = plot_df.copy()

            if facet_rows is not None:
                panel_df = panel_df[panel_df[facet_rows] == row_val]
            if facet_cols is not None:
                panel_df = panel_df[panel_df[facet_cols] == col_val]

            for line_val in line_levels:
                sub = panel_df[panel_df[line_group] == line_val].copy()
                if sub.empty:
                    continue

                sub["_xpos"] = sub[x].map(x_pos)
                sub = sub.sort_values("_xpos")

                color = None if colors is None else colors.get(line_val)
                marker = "o" if markers is None else markers.get(line_val, "o")
                linestyle = "-" if linestyles is None else linestyles.get(line_val, "-")

                ax.errorbar(
                    sub["_xpos"],
                    sub[y],
                    yerr=None if yerr is None else sub[yerr],
                    label=str(line_val),
                    color=color,
                    marker=marker if show_points else None,
                    linestyle=linestyle,
                    linewidth=linewidth,
                    markersize=markersize,
                    capsize=capsize,
                )

            panel_title_parts = []
            if facet_rows is not None:
                panel_title_parts.append(str(row_val))
            if facet_cols is not None:
                panel_title_parts.append(str(col_val))
            if panel_title_parts:
                ax.set_title(" | ".join(panel_title_parts))

            ax.set_xticks(range(len(x_levels)))
            ax.set_xticklabels(x_levels, rotation=rotate_xticks, ha="right")
            ax.set_yscale(yscale)

            if ylim is not None:
                ax.set_ylim(*ylim)
            if xlim is not None:
                ax.set_xlim(*xlim)

            if grid:
                ax.grid(alpha=0.3)

            if j == 0:
                if facet_rows is not None:
                    ax.set_ylabel(f"{row_val}\n{ylabel}")
                else:
                    ax.set_ylabel(ylabel)

            if i == nrows - 1:
                ax.set_xlabel(xlabel if xlabel is not None else x)

    handles, labels = axes[0, 0].get_legend_handles_labels()

    fig.suptitle(title, fontsize=18)
    fig.subplots_adjust(
        left=left_margin,
        right=right_margin,
        top=top_margin,
        bottom=bottom_margin,
        wspace=wspace,
        hspace=hspace,
    )

    if legend and handles:
        fig.legend(
            handles,
            labels,
            title=legend_title if legend_title is not None else line_group,
            loc="center left",
            bbox_to_anchor=(legend_anchor_x, 0.5),
            bbox_transform=fig.transFigure,
            frameon=False,
        )

    return fig, axes