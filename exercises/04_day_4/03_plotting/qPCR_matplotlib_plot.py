# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 12:09:28 2026

@author: wiewe372
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

# Reproducibility
np.random.seed(42)

# True numeric x positions (correct distances)
x = np.array([1, 6, 11, 13, 18, 23])

# Labels shown on the axis
x_labels = ["1", "6", "11", "1", "6", "11"]


# Base expression patterns (synthetic, biology-inspired)
base_data = {
    "White": {
        "WT": {
            "GAPR4": np.array([1.2, 1.8, 2.4, 2.1, 1.4, 1.0]),
            "GAP1":  np.array([0.9, 1.0, 1.2, 1.5, 1.8, 1.3]),
        },
        "gapr4 mutant": {
            "GAPR4": np.array([0.15, 0.10, 0.12, 0.10, 0.14, 0.11]),
            "GAP1":  np.array([0.8, 0.9, 1.0, 1.2, 1.4, 1.1]),
        },
    },
    "Red": {
        "WT": {
            "GAPR4": np.array([1.0, 1.4, 1.9, 1.7, 1.2, 0.9]),
            "GAP1":  np.array([1.1, 1.2, 1.3, 1.4, 1.6, 1.2]),
        },
        "gapr4 mutant": {
            "GAPR4": np.array([0.10, 0.08, 0.10, 0.09, 0.11, 0.10]),
            "GAP1":  np.array([0.9, 1.0, 1.0, 1.1, 1.2, 1.0]),
        },
    },
    "Blue": {
        "WT": {
            "GAPR4": np.array([1.5, 2.2, 3.0, 2.6, 1.8, 1.3]),
            "GAP1":  np.array([1.0, 1.1, 1.3, 1.7, 2.0, 1.5]),
        },
        "gapr4 mutant": {
            "GAPR4": np.array([0.12, 0.10, 0.09, 0.08, 0.10, 0.11]),
            "GAP1":  np.array([0.85, 0.95, 1.05, 1.25, 1.5, 1.2]),
        },
    },
}

# Generate 3 biological replicates
replicate_data = {}
n_reps = 3

for light in base_data:
    replicate_data[light] = {}
    for genotype in base_data[light]:
        replicate_data[light][genotype] = {}
        for gene in base_data[light][genotype]:
            base = base_data[light][genotype][gene]

            # Add noise: a bit proportional to the signal, plus a small baseline
            noise_sd = 0.08 * base + 0.03
            reps = np.array([
                np.clip(base + np.random.normal(0, noise_sd), 0, None)
                for _ in range(n_reps)
            ])

            replicate_data[light][genotype][gene] = reps

# Compute mean and standard deviation
mean_data = {}
sd_data = {}

for light in replicate_data:
    mean_data[light] = {}
    sd_data[light] = {}
    for genotype in replicate_data[light]:
        mean_data[light][genotype] = {}
        sd_data[light][genotype] = {}
        for gene in replicate_data[light][genotype]:
            reps = replicate_data[light][genotype][gene]
            mean_data[light][genotype][gene] = reps.mean(axis=0)
            sd_data[light][genotype][gene] = reps.std(axis=0, ddof=1)

# -------------------------------------------------------
# togehter (1 x 3)


# Fixed layout settings for consistency
RIGHT_MARGIN = 0.82
LEGEND_X = 0.845

fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

colors = {
    "GAPR4 (WT)": "tab:blue",
    "GAP1 (WT)": "tab:orange",
    "GAPR4 (gapr4 mutant)": "tab:green",
    "GAP1 (gapr4 mutant)": "tab:red",
}

styles = {
    "GAPR4 (WT)": "-o",
    "GAP1 (WT)": "-s",
    "GAPR4 (gapr4 mutant)": "--o",
    "GAP1 (gapr4 mutant)": "--s",
}

for ax, light in zip(axes, ["White", "Red", "Blue"]):
    series = {
        "GAPR4 (WT)": mean_data[light]["WT"]["GAPR4"],
        "GAP1 (WT)": mean_data[light]["WT"]["GAP1"],
        "GAPR4 (gapr4 mutant)": mean_data[light]["gapr4 mutant"]["GAPR4"],
        "GAP1 (gapr4 mutant)": mean_data[light]["gapr4 mutant"]["GAP1"],
    }
    errors = {
        "GAPR4 (WT)": sd_data[light]["WT"]["GAPR4"],
        "GAP1 (WT)": sd_data[light]["WT"]["GAP1"],
        "GAPR4 (gapr4 mutant)": sd_data[light]["gapr4 mutant"]["GAPR4"],
        "GAP1 (gapr4 mutant)": sd_data[light]["gapr4 mutant"]["GAP1"],
    }

    for label in series:
        linestyle = "--" if "mutant" in label else "-"
        marker = "o" if "GAPR4" in label else "s"

        ax.errorbar(
            x, series[label], yerr=errors[label],
            color=colors[label],
            linestyle=linestyle,
            marker=marker,
            linewidth=2,
            markersize=6,
            capsize=4,
            label=label
        )

    ax.set_title(f"{light} light")
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=0, ha="right")
    ax.set_xlabel("Time")
    ax.grid(alpha=0.3)

axes[0].set_ylabel("Relative expression")
fig.suptitle("Synthetic qPCR data for GAPR4 and GAP1 (mean ± SD, n=3)", fontsize=20)

handles, labels = axes[0].get_legend_handles_labels()

fig.subplots_adjust(
    left=0.08, right=RIGHT_MARGIN, top=0.82, bottom=0.23,
    wspace=0.28
)

fig.legend(
    handles, labels,
    loc="center left",
    bbox_to_anchor=(LEGEND_X, 0.5),
    bbox_transform=fig.transFigure,
    frameon=False
)
plt.savefig("synthetic_qpcr_replicates.png", dpi=300, bbox_inches="tight")
plt.show()


# ------------------------------------------------------------------
# seperate (2 x 3)

# Fixed layout settings for consistency
RIGHT_MARGIN = 0.82
LEGEND_X = 0.845

lights = ["White", "Red", "Blue"]
genes = ["GAPR4", "GAP1"]
markers = {"WT": "o", "gapr4 mutant": "s"}
linestyles = {"WT": "-", "gapr4 mutant": "--"}

# Colors for condition bars
light_phase_colors = {
    "White": "lightgray",   # or "white" but gray is visible
    "Red": "orangered",
    "Blue": "royalblue",
}

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey="row")

for col, light in enumerate(lights):
    for row, gene in enumerate(genes):
        ax = axes[row, col]

        for genotype in ["WT", "gapr4 mutant"]:
            y = mean_data[light][genotype][gene]
            yerr = sd_data[light][genotype][gene]

            ax.errorbar(
                x, y, yerr=yerr,
                marker=markers[genotype],
                linestyle=linestyles[genotype],
                linewidth=2,
                markersize=6,
                capsize=4,
                label=genotype
            )
            
        ax.set_title(f"{light} light", pad=25)
        ax.set_xlim(0, 24)
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels)
        ax.set_xlabel("Time")
        ax.grid(alpha=0.3)

        # Add top phase bar using data coords for x and axes coords for y
        trans = blended_transform_factory(ax.transData, ax.transAxes)
        bar_y = 1.02
        bar_h = 0.06

        ax.add_patch(Rectangle(
            (0, bar_y), 12, bar_h,
            transform=trans,
            facecolor=light_phase_colors[light],
            edgecolor="none",
            clip_on=False,
            zorder=10
        ))

        ax.add_patch(Rectangle(
            (12, bar_y), 12, bar_h,
            transform=trans,
            facecolor="black",
            edgecolor="none",
            clip_on=False,
            zorder=10
        ))

        if col == 0:
            ax.set_ylabel(f"{gene}\nRelative expression")

        # Optional phase labels only on bottom row
        if row == 1:
            ax.text(6, -0.18, "light phase",
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=9)
            ax.text(18, -0.18, "dark phase",
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=9)

# Shared legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.92, 0.5), frameon=False)

fig.suptitle("Synthetic qPCR data for GAPR4 and GAP1 (mean ± SD, n=3)", fontsize=18)
fig.subplots_adjust(left=0.08, right=0.9, top=0.86, bottom=0.14, wspace=0.35, hspace=0.35)
#plt.tight_layout(rect=[0, 0.03, 0.9, 0.95])
plt.savefig("synthetic_qpcr_2x3_panels.png", dpi=300, bbox_inches="tight")
plt.show()


# only top plot has the bars 
# --------------------------------------------

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey="row")

for col, light in enumerate(lights):
    for row, gene in enumerate(genes):
        ax = axes[row, col]

        for genotype in ["WT", "gapr4 mutant"]:
            y = mean_data[light][genotype][gene]
            yerr = sd_data[light][genotype][gene]

            ax.errorbar(
                x, y, yerr=yerr,
                marker=markers[genotype],
                linestyle=linestyles[genotype],
                linewidth=2,
                markersize=6,
                capsize=4,
                label=genotype
            )
        
        ax.grid(alpha=0.3)
        
        # Optional phase labels only on bottom row
        if row == 1:
            ax.text(6, -0.18, "light phase",
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=9)
            ax.text(18, -0.18, "dark phase",
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=9)
            ax.set_xlabel("Time", labelpad=15)
            
            
        if row == 0:
            ax.set_title(f"{light} light", pad=25, x = 0.25)
            ax.set_xlim(0, 24)
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels)
            
    
            # Add top phase bar using data coords for x and axes coords for y
            trans = blended_transform_factory(ax.transData, ax.transAxes)
            bar_y = 1.02
            bar_h = 0.06
    
            ax.add_patch(Rectangle(
                (0, bar_y), 12, bar_h,
                transform=trans,
                facecolor=light_phase_colors[light],
                edgecolor="none",
                clip_on=False,
                zorder=10
            ))
    
            ax.add_patch(Rectangle(
                (12, bar_y), 12, bar_h,
                transform=trans,
                facecolor="black",
                edgecolor="none",
                clip_on=False,
                zorder=10
            ))

        if col == 0:
            ax.set_ylabel(f"{gene}\nRelative expression")

        
# Shared legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.92, 0.5), frameon=False)

fig.suptitle("Synthetic qPCR data for GAPR4 and GAP1 (mean ± SD, n=3)", fontsize=18)
fig.subplots_adjust(left=0.08, right=0.9, top=0.86, bottom=0.14, wspace=0.20, hspace=0.1)
#plt.tight_layout(rect=[0, 0.03, 0.9, 0.95])
plt.savefig("synthetic_qpcr_2x3_panels_one_bar.png", dpi=300, bbox_inches="tight")
plt.show()


# ----------------------------------------
# circadian test



# time points
x = np.array([1, 6, 11, 13, 18, 23])
theta = 2 * np.pi * x / 24

# example data (replace with yours)
y = np.array([1.2, 1.8, 2.3, 2.2, 1.4, 1.0])

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111, polar=True)

# plot line
ax.plot(theta, y, marker='o')

# close the circle (important!)
ax.plot(np.append(theta, theta[0]),
        np.append(y, y[0]))

# set labels
ax.set_xticks(2 * np.pi * np.array([0, 6, 12, 18]) / 24)
ax.set_xticklabels(["0h", "6h", "12h", "18h"])

ax.set_title("Circadian expression (GAPR4)")

# light phase (0–12h)
ax.bar(
    x=0,
    height=max(y)*1.1,
    width=np.pi,
    bottom=0,
    color="yellow",
    alpha=0.1,
    align="edge"
)

# dark phase (12–24h)
ax.bar(
    x=np.pi,
    height=max(y)*1.1,
    width=np.pi,
    bottom=0,
    color="black",
    alpha=0.1,
    align="edge"
# %%
)


plt.show()


