import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from matplotlib.lines import Line2D
from setuptools.command.rotate import rotate

matplotlib.use("MacOSX")

# Set overall font scale
sns.set_context("talk", font_scale=1.2)

scenarios = ["climdata_era5_cmip6_1981-2100_ssp245_no_cc", "climdata_era5_cmip6_1981-2100_ssp245"]
base_dir = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua") # /a5/forest0/
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc")

# ----- Overall abundance and richness

div = pd.read_csv(base_dir / "a5_diversity_cc_vs_no_cc.csv", index_col=0)

# Convert to longer df
div_long = pd.melt(
    div,
    id_vars=['Year', 'Scenario', 'ForestID', 'SpeciesPool'],
    value_vars=['Richness', 'Abundance'],
    var_name='metric_type',
    value_name='value'
)

# Ts summary
summary = (
    div_long.groupby(["Year", "Scenario", "metric_type"])
    .agg(
        mean = ("value", "mean"),
        sem=("value", "sem")
    )
    .reset_index()
)
# 95% confidence intervals
summary["ci_lower"] = summary["mean"] - 1.96 * summary["sem"]
summary["ci_upper"] = summary["mean"] + 1.96 * summary["sem"]

# ----- Height of vertical maximum abundance and richness

vertical_div = pd.read_csv(base_dir / "a5_vertical_diversity_cc_vs_no_cc.csv", index_col=0)

div_peak = (
    vertical_div.loc[
        vertical_div.groupby(['Scenario', 'Year', 'ForestID', 'SpeciesPool', 'level_5'])['Count'].idxmax(),
        ['Scenario', 'Year', 'ForestID', 'SpeciesPool', 'level_5', 'Count', 'Height']
    ]
    .reset_index(drop=True)
)

# Filter groups with at least 2 species in Richness
valid_groups = (
    div_peak.query("level_5 == 'Richness' and Count >= 2")
    [['Scenario', 'Year', 'ForestID', 'SpeciesPool']]
)
filtered_div_peak = div_peak.merge(valid_groups, on=['Scenario', 'Year', 'ForestID', 'SpeciesPool'])
filtered_div_peak = div_peak[div_peak["level_5"].isin(["Abundance", "Richness"])]

# Ts summary
summary_vertical = (
    filtered_div_peak.groupby(["Year", "Scenario", "level_5"])
    .agg(
        mean = ("Height", "mean"),
        sem=("Height", "sem")
    )
    .reset_index()
)
# 95% confidence intervals
summary_vertical["ci_lower"] = summary_vertical["mean"] - 1.96 * summary_vertical["sem"]
summary_vertical["ci_upper"] = summary_vertical["mean"] + 1.96 * summary_vertical["sem"]

# Create consisten metric type column for merging
summary_vertical["metric_type"] = summary_vertical["level_5"].map({
    "Abundance": "Max. abundance\nheight (m)",
     "Richness": "Max. richness\nheight (m)"
})
# Remove orifinal level_5 column
summary_vertical.drop(["level_5"], axis=1, inplace=True)

# ----- Get species range and species position

species_distr = pd.read_csv(base_dir / "a5_species_distribution_cc_vs_no_cc.csv")

# This computes for each species (within each pool/forest/year/scenario) the mean height and the within-species range
per_species = (
    species_distr
    .groupby(['SpeciesID', 'SpeciesPool', 'ForestID', 'Year', 'Scenario'], observed=True)
    .agg(
        Position = ('Height', 'mean'),
        Range = ('Height', lambda x: x.max() - x.min()),
    )
    .reset_index()
)

# Convert to long format
per_species_long = pd.melt(
    per_species,
    id_vars=['SpeciesID', 'SpeciesPool', 'ForestID', 'Year', 'Scenario'],
    value_vars=['Position', 'Range'],
    var_name='metric_type',
    value_name='value'
)

# Summarize for ts
summary_species = (
        per_species_long
        .groupby(['Scenario', 'Year', 'metric_type'], observed=True)["value"]
        .agg(['mean', 'sem'])
        .reset_index()
)
# 95% CI
summary_species['ci_lower'] = summary_species['mean'] - 1.96 * summary_species['sem']
summary_species['ci_upper'] = summary_species['mean'] + 1.96 * summary_species['sem']

# Merge overall and vertical summaries
plot_data = pd.concat([summary, summary_vertical, summary_species], ignore_index=True)

plot_data.dropna()

# ----- Generate plot

#plt.rcParams.update({'font.size': 12})

# Define colors
colors = {'CC': '#f7766e', 'No CC': '#004aad'}

xticks = range(1980, 2100 + 1, 10)

# Create a figure with 2 columns and 3 rows
fig, axes = plt.subplots(3, 2, figsize=(19, 14), sharex=True)

# Flatten axes for easy iteration
axes = axes.flatten()

# Define the order of metrics for subplots
metrics = ['Richness', 'Abundance', 'Max. richness\nheight (m)',
           'Max. abundance\nheight (m)', 'Position', 'Range']

# Plot each metric in its respective subplot
for i, metric in enumerate(metrics):
    ax = axes[i]
    for scenario, group in plot_data[plot_data['metric_type'] == metric].groupby('Scenario'):
        ax.plot(group['Year'], group['mean'], color=colors[scenario], label=scenario, linewidth=2)
        ax.fill_between(group['Year'], group['ci_lower'], group['ci_upper'], color=colors[scenario], alpha=0.2)

    ax.set_ylabel(metric, fontsize=30)
    ax.grid(alpha=0.3)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation=45)

# Create a single legend
handles, labels = axes[0].get_legend_handles_labels()
labels = ["Climate change", "Baseline"]
fig.legend(handles, labels, title=None, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0))

# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0.05, 1, 1])

# Save the figure
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc")
plt.savefig(DirectoryPlots / 'a_species_comm_resp_v2.pdf')

