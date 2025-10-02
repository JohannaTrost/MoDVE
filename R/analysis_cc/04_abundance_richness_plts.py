import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
matplotlib.use("MacOSX")

# Set overall font scale
sns.set_context("talk", font_scale=1.2)

scenarios = ["climdata_era5_cmip6_1981-2100_ssp245_no_cc", "climdata_era5_cmip6_1981-2100_ssp245"]
base_dir = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua") # /a5/forest0/
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc")

# Load data
data = pd.read_csv(base_dir / "a5_species_distribution_cc_vs_no_cc.csv")

# Compute abundance and richness
div = (
    data.groupby(["Scenario", "Year", "ForestID", "SpeciesPool"])
    .agg(
        Richness=("SpeciesID", "nunique"),   # unique species count
        Abundance=("SpeciesID", "count")     # number of rows = individuals
    )
    .reset_index()
)

div.to_csv(base_dir / "a5_diversity_cc_vs_no_cc.csv")

# ---- Plot individual replicates ----
# Richness over time
g = sns.FacetGrid(
    div,
    col="ForestID",
    row="SpeciesPool",
    hue="Scenario",
    margin_titles=True,
    height=3,
    aspect=1.5
)
g.map_dataframe(sns.lineplot, x="Year", y="Richness")
g.add_legend()
plt.show()
plt.savefig(DirectoryPlots / f"richness_facet.pdf")
# Abundance plot
g = sns.FacetGrid(
    div,
    col="ForestID",
    row="SpeciesPool",
    hue="Scenario",
    margin_titles=True,
    height=3,
    aspect=1.5
)
g.map_dataframe(sns.lineplot, x="Year", y="Abundance")
g.add_legend()
plt.show()
plt.savefig(DirectoryPlots / f"abundance_facet.pdf")

# ---- Plot each Species pool ----

sns.set_context("talk", font_scale=1)

### Abundance ###

# - Step 1: average within each forest (if replicates exist) ---
forest_means = (
    div.groupby(["Year", "ForestID", "SpeciesPool", "Scenario"])
    .agg(mean_abundance=("Abundance", "mean"))
    .reset_index()
)

# - Step 2: compute across-forest mean + CI ---
summary = (
    forest_means.groupby(["Year", "SpeciesPool", "Scenario"])
    .agg(
        mean_abundance=("mean_abundance", "mean"),
        sem=("mean_abundance", "sem")
    )
    .reset_index()
)

# 95% confidence intervals
summary["ci_lower"] = summary["mean_abundance"] - 1.96 * summary["sem"]
summary["ci_upper"] = summary["mean_abundance"] + 1.96 * summary["sem"]

# - Step 3: Plot ---
species_pools = summary["SpeciesPool"].unique()
n_pools = len(species_pools)

fig, axes = plt.subplots(
    nrows=4, ncols=3,
    figsize=(12, 12),
    sharex=True, sharey=True
)

axes = axes.flatten()

# define colors per scenario
scenarios = summary["Scenario"].unique()
colors = plt.cm.tab10(np.arange(len(scenarios)))

for ax, (pool, df_pool) in zip(axes, summary.groupby("SpeciesPool")):
    handles = []
    for color, (scenario, df_scen) in zip(colors, df_pool.groupby("Scenario")):
        # line
        line, = ax.plot(
            df_scen["Year"], df_scen["mean_abundance"],
            color=color, label=scenario
        )
        # CI shading
        patch = ax.fill_between(
            df_scen["Year"],
            df_scen["ci_lower"],
            df_scen["ci_upper"],
            color=color, alpha=0.3
        )
        # store handles for legend (line + patch)
        handles.append((line, Patch(facecolor=color, alpha=0.3)))

    ax.set_title(f"Species Pool: {pool}")
    ax.set_ylabel("Abundance")

# Custom legend: one entry with both line + patch
handles = [h[0] for h in handles] + [h[1] for h in handles]
labels = ([f"{s} (mean)" for s in df_pool["Scenario"].unique()] +
          [f"{s} (95% CI)" for s in df_pool["Scenario"].unique()])
fig.legend(handles, labels, loc='lower right', bbox_to_anchor=(0.8, 0.07))

# Remove empty axes if fewer than 10 pools
for ax in axes[n_pools:]:
    ax.axis("off")

axes[-2].set_xlabel("Year")
axes[-1].set_xlabel("Year")
plt.tight_layout()
plt.savefig(DirectoryPlots / "abundance_mean_across_forests_CI.pdf")
plt.show()

### Richness ###

# - Step 1: average within each forest (if replicates exist) ---
forest_means = (
    div.groupby(["Year", "ForestID", "SpeciesPool", "Scenario"])
    .agg(mean_richness=("Richness", "mean"))
    .reset_index()
)

# - Step 2: compute across-forest mean + CI ---
summary = (
    forest_means.groupby(["Year", "SpeciesPool", "Scenario"])
    .agg(
        mean_richness=("mean_richness", "mean"),
        sem=("mean_richness", "sem")
    )
    .reset_index()
)

# 95% confidence intervals
summary["ci_lower"] = summary["mean_richness"] - 1.96 * summary["sem"]
summary["ci_upper"] = summary["mean_richness"] + 1.96 * summary["sem"]

# - Step 3: Plot ---
species_pools = summary["SpeciesPool"].unique()
n_pools = len(species_pools)

fig, axes = plt.subplots(
    nrows=4, ncols=3,
    figsize=(12, 12),
    sharex=True, sharey=True
)

axes = axes.flatten()

# define colors per scenario
scenarios = summary["Scenario"].unique()
colors = plt.cm.tab10(np.arange(len(scenarios)))

for ax, (pool, df_pool) in zip(axes, summary.groupby("SpeciesPool")):
    handles = []
    for color, (scenario, df_scen) in zip(colors, df_pool.groupby("Scenario")):
        # line
        line, = ax.plot(
            df_scen["Year"], df_scen["mean_richness"],
            color=color, label=scenario
        )
        # CI shading
        patch = ax.fill_between(
            df_scen["Year"],
            df_scen["ci_lower"],
            df_scen["ci_upper"],
            color=color, alpha=0.3
        )
        # store handles for legend (line + patch)
        handles.append((line, Patch(facecolor=color, alpha=0.3)))

    ax.set_title(f"Species Pool: {pool}")
    ax.set_ylabel("Richness")

# Custom legend: one entry with both line + patch
handles = [h[0] for h in handles] + [h[1] for h in handles]
labels = ([f"{s} (mean)" for s in df_pool["Scenario"].unique()] +
          [f"{s} (95% CI)" for s in df_pool["Scenario"].unique()])
fig.legend(handles, labels, loc='lower right', bbox_to_anchor=(0.8, 0.07))

# Remove empty axes if fewer than 10 pools
for ax in axes[n_pools:]:
    ax.axis("off")

axes[-2].set_xlabel("Year")
axes[-1].set_xlabel("Year")
plt.tight_layout()
plt.savefig(DirectoryPlots / "richness_mean_across_forests_CI.pdf")
plt.show()

# ---- Plot each Forest ----

### Abundance ###

# - Step 1: average within each forest (if replicates exist) ---
forest_means = (
    div.groupby(["Year", "SpeciesPool", "ForestID", "Scenario"])
    .agg(mean_abundance=("Abundance", "mean"))
    .reset_index()
)

# - Step 2: compute across-forest mean + CI ---
summary = (
    forest_means.groupby(["Year", "ForestID", "Scenario"])
    .agg(
        mean_abundance=("mean_abundance", "mean"),
        sem=("mean_abundance", "sem")
    )
    .reset_index()
)

# 95% confidence intervals
summary["ci_lower"] = summary["mean_abundance"] - 1.96 * summary["sem"]
summary["ci_upper"] = summary["mean_abundance"] + 1.96 * summary["sem"]

# - Step 3: Plot ---
species_pools = summary["ForestID"].unique()
n_pools = len(species_pools)

fig, axes = plt.subplots(
    nrows=3, ncols=1,
    figsize=(7, 10),
    sharex=True, sharey=True
)

axes = axes.flatten()

# define colors per scenario
scenarios = summary["Scenario"].unique()
colors = plt.cm.tab10(np.arange(len(scenarios)))

for ax, (pool, df_pool) in zip(axes, summary.groupby("ForestID")):
    handles = []
    for color, (scenario, df_scen) in zip(colors, df_pool.groupby("Scenario")):
        # line
        line, = ax.plot(
            df_scen["Year"], df_scen["mean_abundance"],
            color=color, label=scenario
        )
        # CI shading
        patch = ax.fill_between(
            df_scen["Year"],
            df_scen["ci_lower"],
            df_scen["ci_upper"],
            color=color, alpha=0.3
        )
        # store handles for legend (line + patch)
        handles.append((line, Patch(facecolor=color, alpha=0.3)))

    ax.set_title(f"Forest: {pool}")
    ax.set_ylabel("Abundance")

# Custom legend: one entry with both line + patch
handles = [h[0] for h in handles] + [h[1] for h in handles]
labels = ([f"{s} (mean)" for s in df_pool["Scenario"].unique()] +
          [f"{s} (95% CI)" for s in df_pool["Scenario"].unique()])
fig.legend(handles, labels, loc='lower right', bbox_to_anchor=(0.8, 0.07))

# Remove empty axes if fewer than 10 pools
for ax in axes[n_pools:]:
    ax.axis("off")

axes[-2].set_xlabel("Year")
axes[-1].set_xlabel("Year")
plt.tight_layout()
plt.savefig(DirectoryPlots / "abundance_mean_across_species_pools_CI.pdf")
plt.show()


### Richness ###

# - Step 1: average within each forest (if replicates exist) ---
forest_means = (
    div.groupby(["Year", "SpeciesPool", "ForestID", "Scenario"])
    .agg(mean_richness=("Richness", "mean"))
    .reset_index()
)

# - Step 2: compute across-forest mean + CI ---
summary = (
    forest_means.groupby(["Year", "ForestID", "Scenario"])
    .agg(
        mean_richness=("mean_richness", "mean"),
        sem=("mean_richness", "sem")
    )
    .reset_index()
)

# 95% confidence intervals
summary["ci_lower"] = summary["mean_richness"] - 1.96 * summary["sem"]
summary["ci_upper"] = summary["mean_richness"] + 1.96 * summary["sem"]

# - Step 3: Plot ---
species_pools = summary["ForestID"].unique()
n_pools = len(species_pools)

fig, axes = plt.subplots(
    nrows=3, ncols=1,
    figsize=(7, 10),
    sharex=True, sharey=True
)

axes = axes.flatten()

# define colors per scenario
scenarios = summary["Scenario"].unique()
colors = plt.cm.tab10(np.arange(len(scenarios)))

for ax, (pool, df_pool) in zip(axes, summary.groupby("ForestID")):
    handles = []
    for color, (scenario, df_scen) in zip(colors, df_pool.groupby("Scenario")):
        # line
        line, = ax.plot(
            df_scen["Year"], df_scen["mean_richness"],
            color=color, label=scenario
        )
        # CI shading
        patch = ax.fill_between(
            df_scen["Year"],
            df_scen["ci_lower"],
            df_scen["ci_upper"],
            color=color, alpha=0.3
        )
        # store handles for legend (line + patch)
        handles.append((line, Patch(facecolor=color, alpha=0.3)))

    ax.set_title(f"Forest: {pool}")
    ax.set_ylabel("Richness")

# Custom legend: one entry with both line + patch
handles = [h[0] for h in handles] + [h[1] for h in handles]
labels = ([f"{s} (mean)" for s in df_pool["Scenario"].unique()] +
          [f"{s} (95% CI)" for s in df_pool["Scenario"].unique()])
fig.legend(handles, labels, loc='lower right', bbox_to_anchor=(0.8, 0.07))

# Remove empty axes if fewer than 10 pools
for ax in axes[n_pools:]:
    ax.axis("off")

axes[-2].set_xlabel("Year")
axes[-1].set_xlabel("Year")
plt.tight_layout()
plt.savefig(DirectoryPlots / "richness_mean_across_species_pools_CI.pdf")
plt.show()
