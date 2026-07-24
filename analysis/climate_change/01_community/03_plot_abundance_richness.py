# -----
# Visualize species richness and abundance comparing climate change and baseline

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

matplotlib.use("MacOSX")

# Set overall font scale
sns.set_context("talk", font_scale=1.2)

# Directories
scenarios = ["climdata_era5_cmip6_1981-2100_ssp245_no_cc", "climdata_era5_cmip6_1981-2100_ssp245"]
base_dir = Path("../modve_data_zenodo/modve_output/regua") # /a5/forest0/
DirectoryPlots = Path("../modve_figs/climate_change")

DirectoryPlots.mkdir(exist_ok=True)

# Load data
data = pd.read_csv(base_dir / "species_distribution_cc_vs_no_cc.csv")

# Compute abundance and richness
div = (
    data.groupby(["Scenario", "Year", "ForestID", "SpeciesPool"])
    .agg(
        Richness=("SpeciesID", "nunique"),   # unique species count
        Abundance=("SpeciesID", "count")     # number of rows = individuals
    )
    .reset_index()
)

div.to_csv(base_dir / "diversity_cc_vs_no_cc.csv")

# ---- Plot all replicates in one ----

# boxplot for 2030
plt.figure(figsize=(4, 6))
sns.boxplot(
    data=div[div["Year"] == 2030],
    x="Scenario",
    y="Richness",
    hue="Scenario",
    notch=True,
    palette={"CC": "#f7766e", "No CC": "#004aad"},
    medianprops=dict(linewidth=1.2, color='black'),
    flierprops=dict(marker='o', markeredgecolor='gray')
)

plt.xlabel("", fontsize=12)
plt.ylabel("", fontsize=13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=25)
#plt.grid(alpha=0.3)
ax = plt.gca()
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
for spine in plt.gca().spines.values(): spine.set_visible(False)
plt.tight_layout()
plt.savefig(DirectoryPlots / f"richness_2030_boxplot_v2.png", dpi=700)

plt.close("all")
sns.boxplot(
    data=div[div["Year"] == 2050],
    x="Scenario",
    y="Abundance",
    hue="Scenario",
    notch=True,
    palette={"CC": "#f7766e", "No CC": "#004aad"},
    medianprops=dict(linewidth=1.2, color='black'),
    flierprops=dict(marker='o', markeredgecolor='gray')
)

plt.xlabel("Scenario", fontsize=12)
plt.ylabel("Abundance", fontsize=13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(alpha=0.3)
for spine in plt.gca().spines.values(): spine.set_visible(False)
plt.savefig(DirectoryPlots / f"abundance_2050_boxplot.png", dpi=700)

# --- Ts summary

summary = (
    div.groupby(["Year", "Scenario"])
    .agg(
        mean_richness = ("Richness", "mean"),
        mean_abundance=("Abundance", "mean"),
        sem_abundance=("Abundance", "sem"),
        sem_richness=("Richness", "sem")
    )
    .reset_index()
)

# 95% confidence intervals
summary["ci_lower_abundance"] = summary["mean_abundance"] - 1.96 * summary["sem_abundance"]
summary["ci_upper_abundance"] = summary["mean_abundance"] + 1.96 * summary["sem_abundance"]
summary["ci_lower_richness"] = summary["mean_richness"] - 1.96 * summary["sem_richness"]
summary["ci_upper_richness"] = summary["mean_richness"] + 1.96 * summary["sem_richness"]

plt.rcParams.update({'font.size': 12})

# Define colors
colors = {'CC': '#f7766e', 'No CC': '#004aad'}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharex=True)

# ---- Abundance subplot ----
for scenario, group in summary.groupby('Scenario'):
    ax1.plot(group['Year'], group['mean_abundance'],
             color=colors[scenario], label=scenario, linewidth=2)
    ax1.fill_between(group['Year'], group['ci_lower_abundance'], group['ci_upper_abundance'],
                     color=colors[scenario], alpha=0.2)

ax1.set_xlabel('')
ax1.set_ylabel('Abundance', fontsize=30)
xticks = range(summary['Year'].min() - 1, summary['Year'].max() + 1, 10)
ax1.set_xticks(xticks)
ax1.set_xticklabels(xticks, rotation=25)
ax1.set_xlim(summary['Year'].min(), summary['Year'].max())
ax1.grid(alpha=0.3)
for spine in ax1.spines.values(): spine.set_visible(False)

# ---- Richness subplot ----
for scenario, group in summary.groupby('Scenario'):
    ax2.plot(group['Year'], group['mean_richness'],
             color=colors[scenario], label=scenario, linewidth=2)
    ax2.fill_between(group['Year'], group['ci_lower_richness'], group['ci_upper_richness'],
                     color=colors[scenario], alpha=0.2)

ax2.set_xlabel('')
ax2.set_ylabel('Richness', fontsize=30)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xticks, rotation=25)
ax2.set_xlim(summary['Year'].min(), summary['Year'].max())
ax2.grid(alpha=0.3)
for spine in ax2.spines.values(): spine.set_visible(False)

# Legend (only one, centered below)
handles, labels = ax1.get_legend_handles_labels()
# Legend (only one, centered below)
handles, labels = ax1.get_legend_handles_labels()
labels = ["Climate change", "Baseline"]
fig.legend(handles, labels, title=None,
           loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0))

# Increase bottom margin so the legend fits
plt.tight_layout(rect=[0, 0.12, 1, 1])
plt.savefig(DirectoryPlots / "abundance_richness_ts.pdf")

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


# --- Compute relative percent differences ---

# 1. Filter for year 2030
df_2030 = div[div['Year'] == 2050]

# 2. Split by scenario
cc = df_2030[df_2030['Scenario'] == 'CC']
nocc = df_2030[df_2030['Scenario'] == 'No CC']

# 3. Merge CC and No CC rows by ForestID + SpeciesPool (or other IDs you need)
merged = pd.merge(
    cc, nocc,
    on=['ForestID', 'SpeciesPool'],
    suffixes=('_CC', '_NoCC')
)

# 4. Compute percent change in richness
merged['pct_change_richness'] = (
    (merged['Richness_CC'] - merged['Richness_NoCC']) /
     merged['Richness_NoCC'] * 100
)
merged['pct_change_abundance'] = (
    (merged['Abundance_CC'] - merged['Abundance_NoCC']) /
     merged['Abundance_NoCC'] * 100
)

# 5. Average percent change
average_pct_change = merged['pct_change_richness'].mean()
avg_pct_change_abundance = merged['pct_change_abundance'].mean()

# CI
print(merged['pct_change_richness'].sem() * 1.96)
print(merged['pct_change_abundance'].sem() * 1.96)