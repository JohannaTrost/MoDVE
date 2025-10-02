import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
from matplotlib.pyplot import ylabel
from scipy.stats import entropy
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
matplotlib.use("MacOSX")

from scipy.signal import savgol_filter


def apply_savgol_filter(group, window_length=5, polyorder=2,
        column_to_filter='Height'):
    """
    Apply Savitzky-Golay filter to a specific column within each group.

    Parameters:
    - group: pandas DataFrame group
    - window_length: int, length of the filter window (must be odd)
    - polyorder: int, order of polynomial to fit
    - column_to_filter: str, name of column to apply filter to

    Returns:
    - group with new column containing filtered values
    """
    # Make a copy to avoid modifying original data
    group = group.copy()

    # Sort by height to ensure proper ordering for the filter
    group = group.sort_values(column_to_filter)

    # Check if we have enough points for the filter
    if len(group) >= window_length:
        # Apply Savitzky-Golay filter
        filtered_values = savgol_filter(group[column_to_filter],
                                        window_length=window_length,
                                        polyorder=polyorder)
        group[f'{column_to_filter}_savgol'] = filtered_values
    else:
        # If not enough points, just copy original values
        print(
            f"Warning: Group has {len(group)} points, less than window_length "
            f"{window_length}. Using original values.")
        group[f'{column_to_filter}_savgol'] = group[column_to_filter]

    return group


scenarios = ["climdata_era5_cmip6_1981-2100_ssp245_no_cc", "climdata_era5_cmip6_1981-2100_ssp245"]
base_dir = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua") # /a5/forest0/
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc")

# Load data
data = pd.read_csv(base_dir / "a5_species_distribution_cc_vs_no_cc.csv")

species_counts = data.groupby(["Scenario", "Year", "ForestID", "SpeciesPool",
                               "SpeciesID", "Height"]).size().reset_index(name='Count')

# Function to calculate alpha diversity metrics
def diversity(counts):
    counts_array = np.array(counts)

    # Handle edge cases
    if len(counts_array) == 0 or counts_array.sum() == 0:
        return pd.Series({
            "Richness": 0,
            "Shannon": 0,
            "Simpson": 0,
            "Evenness": 0,
            "Abundance": 0
        })

    # Filter out zero counts for calculations
    non_zero_counts = counts_array[counts_array > 0]
    richness = len(non_zero_counts)

    # Shannon diversity
    shannon = entropy(counts_array, base=np.e)

    # Simpson index
    proportions = counts_array / counts_array.sum()
    simpson = 1 - np.sum(proportions ** 2)

    # Evenness (Pielou's evenness)
    if richness > 1:
        evenness = shannon / np.log(richness)
    elif richness == 1:
        evenness = 1.0  # Perfect evenness with one species
    else:
        evenness = 0

    return pd.Series({
        "Richness": richness,
        "Shannon": shannon,
        "Simpson": simpson,
        "Evenness": evenness,
        "Abundance": int(counts_array.sum())
    })

# Compute alpha diversity per forest
div = species_counts.groupby(["Scenario", "Year", "ForestID", "SpeciesPool", "Height"])['Count'].apply(diversity).reset_index()

# - Plot vertical diversity profiles averaged last 10 years (showing mean and CI), facet over forests and species pools
div_last10 = div[div['Year'] >= 2079]

def plot_div_profiles(div_data, div_metric):

    div_metric_last10 = div_data[div_data["level_5"] == div_metric]
    div_metric_median = div_metric_last10.groupby(
        ["Scenario", "ForestID", "SpeciesPool", "Height"]).agg(
        Median=('Count', 'median'),
        SD=('Count', 'std')
    ).reset_index()

    # Calculate confidence intervals
    div_metric_median = div_metric_median.copy()
    div_metric_median['CI_lower'] = div_metric_median['Median'] - 1.96 * \
                                  div_metric_median['SD']
    div_metric_median['CI_upper'] = div_metric_median['Median'] + 1.96 * \
                                  div_metric_median['SD']

    div_metric_median_smooth = div_metric_median.groupby(
        ['ForestID', 'SpeciesPool', 'Scenario']).apply(
        lambda x: apply_savgol_filter(x, window_length=5, polyorder=2,
                                      column_to_filter='Median')
    ).reset_index(drop=True)

    # Create FacetGrid
    g = sns.FacetGrid(div_metric_median_smooth,
                      col='ForestID', row='SpeciesPool',
                      height=4, aspect=0.8,
                      margin_titles=True)

    # Plot lines and confidence intervals
    def plot_profile(x, y, **kwargs):
        data = kwargs.pop('data')
        scenarios = data['Scenario'].unique()
        colors = {'CC': '#1f77b4', 'No CC': '#ff7f0e'}

        for scenario in scenarios:
            scenario_data = data[data['Scenario'] == scenario].sort_values(
                'Height')
            if not scenario_data.empty:
                plt.plot(scenario_data['Median'],
                         scenario_data['Height'],
                         color=colors.get(scenario, 'black'),
                         label=scenario, linewidth=2, marker='.', markersize=4)
                # Plot smoothed line
                plt.plot(scenario_data['Median_savgol'],
                         scenario_data['Height'],
                         color=colors.get(scenario, 'black'),
                         label=scenario, linewidth=1.5, linestyle='--')

                # Add confidence interval
                plt.fill_betweenx(scenario_data['Height'],
                                  scenario_data['CI_lower'],
                                  scenario_data['CI_upper'],
                                  color=colors.get(scenario, 'gray'), alpha=0.2)

    g.map_dataframe(plot_profile, 'Median', 'Height')
    g.set_axis_labels(f'{div_metric} Median', 'Height (m)')
    g.add_legend(title='Scenario')

    # Add grid
    for ax in g.axes.flat:
        ax.grid(True, alpha=0.3)

    return g

# plot and save for each metric
metrics = ["Richness", "Shannon", "Simpson", "Evenness", "Abundance"]
for metric in metrics:
    plt.close("all")
    g = plot_div_profiles(div_last10, metric)
    plt.show()
    plt.savefig(DirectoryPlots / f"Vertival_{metric}_facet_last20.pdf")


# Analyse peak in each time step, where at least 2 species are present

div_peak = (
    div.loc[
        div.groupby(['Scenario', 'Year', 'ForestID', 'SpeciesPool', 'level_5'])['Count'].idxmax(),
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

# ------ Facet over forests and species pools ------ #

# --- Peak height over time, facet over forests and species pools
for metric in metrics:
    plt.close("all")
    g = sns.FacetGrid(
        filtered_div_peak[filtered_div_peak["level_5"] == metric],
        col="ForestID",
        row="SpeciesPool",
        hue="Scenario",
        margin_titles=True,
        height=3,
        aspect=1.5
    )
    g.map_dataframe(sns.lineplot, x="Year", y="Count")
    g.set_ylabels(f"{metric} peak")
    g.add_legend()
    plt.show()
    plt.savefig(DirectoryPlots / f"vertical_{metric.lower()}_peak_facet_ts.pdf")

# --- Peak value over time
for metric in metrics:
    plt.close("all")
    g = sns.FacetGrid(
        filtered_div_peak[filtered_div_peak["level_5"] == metric],
        col="ForestID",
        row="SpeciesPool",
        hue="Scenario",
        margin_titles=True,
        height=3,
        aspect=1.5
    )
    g.map_dataframe(sns.lineplot, x="Year", y="Height")
    g.set_ylabels(f"{metric} peak\n Height (m)")
    g.add_legend()
    plt.show()
    plt.savefig(DirectoryPlots / f"vertical_{metric.lower()}_peak_height_facet_ts.pdf")

# ------ Facet over species pools only ------ #

# - Step 1: average within each forest (if replicates exist) ---
forest_means = (
    filtered_div_peak.groupby(["Year", "ForestID", "SpeciesPool", "Scenario", "level_5"])
    .agg(mean_div=("Count", "mean"),
         mean_div_hgt=("Height", "mean"))
    .reset_index()
)

# - Step 2: compute across-forest mean + CI ---
summary = (
    forest_means.groupby(["Year", "SpeciesPool", "Scenario", "level_5"])
    .agg(
        mean_div=("mean_div", "mean"),
        sem_div=("mean_div", "sem"),
        mean_div_hgt=("mean_div_hgt", "mean"),
        sem_div_hgt=("mean_div_hgt", "sem")
    )
    .reset_index()
)

# 95% confidence intervals
summary["div_ci_lower"] = summary["mean_div"] - 1.96 * summary["sem_div"]
summary["div_ci_upper"] = summary["mean_div"] + 1.96 * summary["sem_div"]
summary["div_hgt_ci_lower"] = summary["mean_div_hgt"] - 1.96 * summary["sem_div_hgt"]
summary["div_hgt_ci_upper"] = summary["mean_div_hgt"] + 1.96 * summary["sem_div_hgt"]

# - Step 3: Plot ---
species_pools = summary["SpeciesPool"].unique()
n_pools = len(species_pools)

# define colors per scenario
scenarios = summary["Scenario"].unique()
colors = plt.cm.tab10(np.arange(len(scenarios)))

# Plot peak value over time
for metric in metrics:
    plt.close("all")

    metric_summary = summary[summary["level_5"] == metric]

    fig, axes = plt.subplots(
        nrows=4, ncols=3,
        figsize=(12, 12),
        sharex=True, sharey=True
    )

    axes = axes.flatten()

    for ax, (pool, df_pool) in zip(axes, metric_summary.groupby("SpeciesPool")):
        handles = []
        for color, (scenario, df_scen) in zip(colors, df_pool.groupby("Scenario")):
            # line
            line, = ax.plot(
                df_scen["Year"], df_scen["mean_div"],
                color=color, label=scenario
            )
            # CI shading
            patch = ax.fill_between(
                df_scen["Year"],
                df_scen["div_ci_lower"],
                df_scen["div_ci_upper"],
                color=color, alpha=0.3
            )
            # store handles for legend (line + patch)
            handles.append((line, Patch(facecolor=color, alpha=0.3)))

        ax.set_title(f"Species Pool: {pool}")
        ax.set_ylabel(f"{metric} peak")

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
    plt.savefig(DirectoryPlots / f"{metric.lower()}_peak_mean_across_forests_CI.pdf")
    plt.show()


# Plot peak height over time
for metric in metrics:
    plt.close("all")

    metric_summary = summary[summary["level_5"] == metric]

    fig, axes = plt.subplots(
        nrows=4, ncols=3,
        figsize=(12, 12),
        sharex=True, sharey=True
    )

    axes = axes.flatten()

    for ax, (pool, df_pool) in zip(axes, metric_summary.groupby("SpeciesPool")):
        handles = []
        for color, (scenario, df_scen) in zip(colors, df_pool.groupby("Scenario")):
            # line
            line, = ax.plot(
                df_scen["Year"], df_scen["mean_div_hgt"],
                color=color, label=scenario
            )
            # CI shading
            patch = ax.fill_between(
                df_scen["Year"],
                df_scen["div_hgt_ci_lower"],
                df_scen["div_hgt_ci_upper"],
                color=color, alpha=0.3
            )
            # store handles for legend (line + patch)
            handles.append((line, Patch(facecolor=color, alpha=0.3)))

        ax.set_title(f"Species Pool: {pool}")
        ax.set_ylabel(f"{metric} peak\nheight (m)")

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
    plt.savefig(DirectoryPlots / f"{metric.lower()}_peak_hgt_mean_across_forests_CI.pdf")
    plt.show()

# ------ Facet over forests only ------ #

# - Step 2: compute across-forest mean + CI ---
summary = (
    forest_means.groupby(["Year", "ForestID", "Scenario", "level_5"])
    .agg(
        mean_div=("mean_div", "mean"),
        sem_div = ("mean_div", "sem"),
        mean_div_hgt = ("mean_div_hgt", "mean"),
        sem_div_hgt = ("mean_div_hgt", "sem")
    )
    .reset_index()
)

# 95% confidence intervals
summary["div_ci_lower"] = summary["mean_div"] - 1.96 * summary["sem_div"]
summary["div_ci_upper"] = summary["mean_div"] + 1.96 * summary["sem_div"]
summary["div_hgt_ci_lower"] = summary["mean_div_hgt"] - 1.96 * summary["sem_div_hgt"]
summary["div_hgt_ci_upper"] = summary["mean_div_hgt"] + 1.96 * summary["sem_div_hgt"]

# - Step 3: Plot ---
species_pools = summary["ForestID"].unique()
n_pools = len(species_pools)

# define colors per scenario
scenarios = summary["Scenario"].unique()
colors = plt.cm.tab10(np.arange(len(scenarios)))

# Plot peak value over time
for metric in metrics:
    plt.close("all")

    metric_summary = summary[summary["level_5"] == metric]

    fig, axes = plt.subplots(
        nrows=3, ncols=1,
        figsize=(7, 10),
        sharex=True, sharey=True
    )

    axes = axes.flatten()

    for ax, (pool, df_pool) in zip(axes, metric_summary.groupby("ForestID")):
        handles = []
        for color, (scenario, df_scen) in zip(colors, df_pool.groupby("Scenario")):
            # line
            line, = ax.plot(
                df_scen["Year"], df_scen["mean_div"],
                color=color, label=scenario
            )
            # CI shading
            patch = ax.fill_between(
                df_scen["Year"],
                df_scen["div_ci_lower"],
                df_scen["div_ci_upper"],
                color=color, alpha=0.3
            )
            # store handles for legend (line + patch)
            handles.append((line, Patch(facecolor=color, alpha=0.3)))

        ax.set_title(f"Forest: {pool}")
        ax.set_ylabel(f"{metric} peak")

    # Custom legend: one entry with both line + patch
    handles = [h[0] for h in handles] + [h[1] for h in handles]
    labels = ([f"{s} (mean)" for s in df_pool["Scenario"].unique()] +
              [f"{s} (95% CI)" for s in df_pool["Scenario"].unique()])
    fig.legend(
        handles, labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),  # push legend below plots
        ncol=2,                       # two-column layout
        frameon=False
    )

    # Remove empty axes if fewer than 10 pools
    for ax in axes[n_pools:]:
        ax.axis("off")

    axes[-2].set_xlabel("Year")
    axes[-1].set_xlabel("Year")
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # leave space at bottom
    plt.savefig(DirectoryPlots / f"{metric}_peak_mean_across_species_pools_CI.pdf")
    plt.show()

# Plot peak height over time
for metric in metrics:
    plt.close("all")

    metric_summary = summary[summary["level_5"] == metric]

    fig, axes = plt.subplots(
        nrows=3, ncols=1,
        figsize=(7, 10),
        sharex=True, sharey=True
    )

    axes = axes.flatten()

    for ax, (pool, df_pool) in zip(axes, metric_summary.groupby("ForestID")):
        handles = []
        for color, (scenario, df_scen) in zip(colors, df_pool.groupby("Scenario")):
            # line
            line, = ax.plot(
                df_scen["Year"], df_scen["mean_div_hgt"],
                color=color, label=scenario
            )
            # CI shading
            patch = ax.fill_between(
                df_scen["Year"],
                df_scen["div_hgt_ci_lower"],
                df_scen["div_hgt_ci_upper"],
                color=color, alpha=0.3
            )
            # store handles for legend (line + patch)
            handles.append((line, Patch(facecolor=color, alpha=0.3)))

        ax.set_title(f"Forest: {pool}")
        ax.set_ylabel(f"{metric} peak\nheight (m)")

    # Custom legend: one entry with both line + patch
    handles = [h[0] for h in handles] + [h[1] for h in handles]
    labels = ([f"{s} (mean)" for s in df_pool["Scenario"].unique()] +
              [f"{s} (95% CI)" for s in df_pool["Scenario"].unique()])
    fig.legend(
        handles, labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),  # push legend below plots
        ncol=2,                       # two-column layout
        frameon=False
    )

    # Remove empty axes if fewer than 10 pools
    for ax in axes[n_pools:]:
        ax.axis("off")

    axes[-2].set_xlabel("Year")
    axes[-1].set_xlabel("Year")
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # leave space at bottom
    plt.savefig(DirectoryPlots / f"{metric}_peak_height_mean_across_species_pools_CI.pdf")
    plt.show()
 
# Save processed data
div.to_csv(base_dir / "a5_vertical_diversity_cc_vs_no_cc.csv")
div_peak.to_csv(base_dir / "a5_vertical_diversity_peak_cc_vs_no_cc.csv")