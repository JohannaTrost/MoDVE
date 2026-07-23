# -----
# Visualize vertical diversity as smooth vertical profile and over time
# comparing climate change and baseline

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
from scipy.stats import entropy
import seaborn as sns
from matplotlib.patches import Patch
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import imageio

matplotlib.use("MacOSX")

from scipy.signal import savgol_filter

plt.rcParams.update({'font.size': 12})

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
data = pd.read_csv(base_dir / "species_distribution_cc_vs_no_cc.csv")

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
div_last10 = div[div['Year'] == 2050]

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
    plt.savefig(DirectoryPlots / f"Vertival_{metric}_facet_2030.pdf")


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
div.to_csv(base_dir / "vertical_diversity_cc_vs_no_cc.csv")
div_peak.to_csv(base_dir / "vertical_diversity_peak_cc_vs_no_cc.csv")

# ------ Overall summary ------ #

overall_peak = div_peak.groupby(["Scenario", "level_5", "Year"]).agg(
        mean_div=("Count", "mean"),
        sem_div = ("Count", "sem"),
        mean_div_hgt= ("Height", "mean"),
        sem_div_hgt = ("Height", "sem")
    ).reset_index()

print(overall_peak[overall_peak["Year"] == 2050])

# ------ Compute percent change in 2050 ------ #

# 1. Filter for year 2030
df_2050 = div_peak[div_peak['Year'] == 2080]

# 2. Split by scenario
cc = df_2050[df_2050['Scenario'] == 'CC']
nocc = df_2050[df_2050['Scenario'] == 'No CC']

# 3. Merge CC and No CC rows by ForestID + SpeciesPool (or other IDs you need)
merged = pd.merge(
    cc, nocc,
    on=['ForestID', 'SpeciesPool', 'Year', 'level_5'],
    suffixes=('_CC', '_NoCC')
)

# 4. Compute percent change in richness
merged['pct_change_div'] = (
    (merged['Count_CC'] - merged['Count_NoCC']) /
     merged['Count_NoCC'] * 100
)
merged['pct_change_height'] = (
    (merged['Height_CC'] - merged['Height_NoCC']) /
     merged['Height_NoCC'] * 100
)

# 5. Average percent change
metric = "Richness"
avg_pct_change_div = merged.loc[merged['level_5'] == metric, 'pct_change_div'].mean()
avg_pct_change_height = merged.loc[merged['level_5'] == metric, 'pct_change_height'].mean()
sem_pct_change_div = merged.loc[merged['level_5'] == metric, 'pct_change_div'].sem()
sem_pct_change_height = merged.loc[merged['level_5'] == metric, 'pct_change_height'].sem()
print(f"Average percent change in diversity (2050): {avg_pct_change_div:.1f} % ± {1.96 * sem_pct_change_div:.1f} %")
print(f"Average percent change in peak height (2050): {avg_pct_change_height:.1f} % ± {1.96 * sem_pct_change_height:.1f} %")

# ------ Overall vertical div all replicates in one ------ #

# - Compute across-forest mean + CI ---
div2030 = div[(div['Year'] == 2030)  & (div['SpeciesPool'] == 3) & (div['ForestID'] == 0)]
summary = (
    div2030.groupby(["Scenario", "level_5", "Height"])
    .agg(
        mean_div=("Count", "mean"),
        sem_div = ("Count", "sem")
    )
    .reset_index()
)

# 95% confidence intervals
summary["div_ci_lower"] = summary["mean_div"] - 1.96 * summary["sem_div"]
summary["div_ci_upper"] = summary["mean_div"] + 1.96 * summary["sem_div"]

data = summary.copy()

# Plot lines and confidence intervals
scenarios = data['Scenario'].unique()
colors = {'CC': '#f7766e', 'No CC': '#004aad'}

plt.close("all")

metric = "Richness"
sigma = 2  # controls smoothing strength

plt.figure(figsize=(6, 8))

for scenario in scenarios:
    scenario_data = data[
        (data['Scenario'] == scenario) &
        (data["level_5"] == metric)
        ].sort_values('Height')

    if not scenario_data.empty:
        # Apply Gaussian smoothing
        smoothed_mean = gaussian_filter1d(scenario_data['mean_div'],
                                          sigma=sigma)
        smoothed_lower = gaussian_filter1d(scenario_data['div_ci_lower'],
                                           sigma=sigma)
        smoothed_upper = gaussian_filter1d(scenario_data['div_ci_upper'],
                                           sigma=sigma)

        color = colors.get(scenario, 'black')

        # --- Plot smoothed mean ---
        plt.plot(smoothed_mean, scenario_data['Height'],
                 color=color, linewidth=2, label=f"{scenario} (smoothed)")

        # --- Plot smoothed CI ---
        plt.fill_betweenx(scenario_data['Height'],
                          smoothed_lower, smoothed_upper,
                          color=color, alpha=0.2)

        # --- Plot raw (non-smoothed) mean for transparency ---
        plt.plot(scenario_data['mean_div'], scenario_data['Height'],
                 color=color, alpha=0.3, linestyle='--', linewidth=1)

plt.xlabel("Species richness")
plt.ylabel("Height (m)")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
           fancybox=True, shadow=False, ncol=1, frameon=False)
plt.grid(alpha=0.3)
for spine in plt.gca().spines.values(): spine.set_visible(False)
plt.tight_layout()
plt.show()
plt.savefig(DirectoryPlots / "vertical_richness_smooth_sp3_f0.png", dpi=700)

# ---- Plot vertical diversity across replicates for each year ---- #

plt.close("all")

# ---- Make GIF for richness over time ---- #

years = np.asarray(sorted(div['Year'].unique()))
years = years[years >= 2030]
metrics = ["Richness", "Abundance"]

# Determine global plotting limits
global_limits = {'Richness':
                     {'xmin': 1.0, 'xmax': 6.0, 'ymin': 0.5, 'ymax': 51.5},
                 'Abundance':
                     {'xmin': 1.0, 'xmax': 50.0, 'ymin': 0.5, 'ymax': 51.5},
                 'Shannon':
                     {'xmin': 0.0, 'xmax': 1.3, 'ymin': 0.5, 'ymax': 51.5}
                 }

out_dir = Path(DirectoryPlots)
out_dir.mkdir(exist_ok=True)

sigma = 2
gif_files = {m: [] for m in metrics}

for metric in metrics:
    for year in years:

        fig, ax = plt.subplots(figsize=(7, 8))

        div_year = div[div['Year'] == year]

        # Filter: richness >=2 combos still applied
        rich = div_year[div_year["level_5"] == "Richness"]
        valid = rich[rich["Count"] >= 2][["SpeciesPool","Scenario"]]
        div_filtered = div_year.merge(valid, on=["SpeciesPool","Scenario"], how="inner")

        summary = (
            div_filtered.groupby(["Scenario", "level_5", "Height"])
            .agg(mean_div=("Count", "mean"), sem_div=("Count", "sem"))
            .reset_index()
        )

        summary["div_ci_lower"] = summary["mean_div"] - 1.96 * summary["sem_div"]
        summary["div_ci_upper"] = summary["mean_div"] + 1.96 * summary["sem_div"]

        data = summary.copy()

        for scenario in scenarios:
            df_s = data[(data["Scenario"] == scenario) & (data["level_5"] == metric)]
            df_s = df_s.sort_values("Height")

            if df_s.empty:
                continue

            sm_mean = gaussian_filter1d(df_s["mean_div"], sigma=sigma)
            sm_low  = gaussian_filter1d(df_s["div_ci_lower"], sigma=sigma)
            sm_up   = gaussian_filter1d(df_s["div_ci_upper"], sigma=sigma)

            color = colors.get(scenario, 'black')

            # Smooth curve
            ax.plot(sm_mean, df_s["Height"], color=color, linewidth=2)

            # CI band
            ax.fill_betweenx(df_s["Height"], sm_low, sm_up, color=color, alpha=0.2)

        ## --- fixed axes for all frames ---
        ax.set_xlim(global_limits[metric]["xmin"], global_limits[metric]["xmax"])
        ax.set_ylim(global_limits[metric]["ymin"], global_limits[metric]["ymax"])

        ax.set_title(f"Shannon diversity index - Year {year}", fontsize=22)
        ax.set_xlabel(metric, fontsize=24)
        ax.set_ylabel("Height (m)", fontsize=24)
        ax.grid(alpha=0.3)

        # remove spines
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Save temporary PNG per frame
        frame_path = out_dir / f"{metric}_{year}.png"
        plt.tight_layout()
        plt.savefig(frame_path, dpi=120)
        gif_files[metric].append(frame_path)

        plt.close(fig)

# ---- Build GIFs ----
for metric in metrics:
    frames = [imageio.imread(p) for p in gif_files[metric]]
    imageio.mimsave(out_dir / f"{metric}_animation.gif", frames, duration=0.5)  # 0.5s per frame

# -------------------- Time series of max and max height ---- #

colors = {'CC': '#f7766e', 'No CC': '#004aad'}  # keep same colors as your example

# --- Aggregate across Species to get time series (mean ± 95% CI) for position and range ---
# For each Year and Scenario we compute the mean of mean_height across species, sem, and 95% CI
def summarize_across_species(df, value_col, prefix):
    summary = (
        df
        .groupby(['Scenario', 'Year', 'level_5'], observed=True)[value_col]
        .agg(['mean', 'std', 'count'])
        .rename(columns={'mean': f'mean_{prefix}', 'std': f'std_{prefix}', 'count': f'n_{prefix}'})
        .reset_index()
    )
    # SEM and 95% CI (approx normal, z=1.96)
    summary[f'sem_{prefix}'] = summary[f'std_{prefix}'] / np.sqrt(summary[f'n_{prefix}'])
    summary[f'ci_lower_{prefix}'] = summary[f'mean_{prefix}'] - 1.96 * summary[f'sem_{prefix}']
    summary[f'ci_upper_{prefix}'] = summary[f'mean_{prefix}'] + 1.96 * summary[f'sem_{prefix}']
    return summary

summary_max_val = summarize_across_species(filtered_div_peak, 'Count', 'Count')
summary_max_hgt = summarize_across_species(filtered_div_peak, 'Height', 'Height')

# For plotting convenience, merge the two summaries on Scenario+Year so we can pick xticks from range of years
# But plotting will iterate over each independently (like your example)
years_min = min(summary_max_val['Year'].min(), summary_max_hgt['Year'].min())
years_max = max(summary_max_val['Year'].max(), summary_max_hgt['Year'].max())

# --- Plotting: 3 metrics × 2 columns = 6 subplots ---

metrics = [
    ("Abundance", "Vertical max.\nabundance", "Max. abundance\nheight (m)"),
    ("Shannon", "Vertical max.\nShannon index", "Max. Shannon index\nheight (m)"),
    ("Richness", "Vertical max.\nrichness", "Max. richness\nheight (m)")
]
yrs = [2050, 2080, 2090]

plt.rcParams.update({'font.size': 20})

fig, axes = plt.subplots(
    nrows=3, ncols=3, figsize=(15, 15)
)

xticks = range(years_min - 1, years_max + 1, 10)

for col, (metric, ylabel_top, ylabel_middle) in enumerate(metrics):

    # ---------------- TOP ROW: maximum value ----------------
    ax_top = axes[0, col]
    data_max = summary_max_val[summary_max_val["level_5"] == metric]

    for scenario, group in data_max.groupby("Scenario"):
        ax_top.plot(group["Year"], group["mean_Count"],
                    color=colors.get(scenario, None),
                    linewidth=2, label=scenario)
        ax_top.fill_between(group["Year"],
                            group["ci_lower_Count"],
                            group["ci_upper_Count"],
                            color=colors.get(scenario, None),
                            alpha=0.2)

    ax_top.set_ylabel(ylabel_top, fontsize=25)
    ax_top.set_xticks(xticks)
    ax_top.set_xticklabels([], )
    ax_top.set_xlim(years_min, years_max)
    ax_top.grid(alpha=0.3)
    for spine in ax_top.spines.values():
        spine.set_visible(False)
    ax_top.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    # ---------------- MIDDLE ROW: height of maximum ----------------
    ax_middle = axes[1, col]   # <--- FIXED
    data_hgt = summary_max_hgt[summary_max_hgt["level_5"] == metric]

    for scenario, group in data_hgt.groupby("Scenario"):
        ax_middle.plot(group["Year"], group["mean_Height"],
                       color=colors.get(scenario, None),
                       linewidth=2)
        ax_middle.fill_between(group["Year"],
                               group["ci_lower_Height"],
                               group["ci_upper_Height"],
                               color=colors.get(scenario, None),
                               alpha=0.2)

    ax_middle.set_ylabel(ylabel_middle, fontsize=25)
    ax_middle.set_xticks(xticks)
    ax_middle.set_xticklabels(xticks, rotation=90, fontsize=18)
    ax_middle.set_xlim(years_min, years_max)
    ax_middle.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    # ---------------- BOTTOM ROW: vertical profile ----------------
    ax_bottom = axes[2, col]
    year = yrs[col]

    div2030 = div[(div['Year'] == year)]

    richness = div2030[div2030['level_5'] == 'Richness']
    valid_combinations = richness[richness['Count'] >= 2][
        ['SpeciesPool', 'Scenario']]

    div_filtered = div2030.merge(valid_combinations,
                                 on=['SpeciesPool', 'Scenario'], how='inner')

    summary = (
        div_filtered.groupby(["Scenario", "level_5", "Height"])
        .agg(
            mean_div=("Count", "mean"),
            sem_div=("Count", "sem")
        )
        .reset_index()
    )

    summary["div_ci_lower"] = summary["mean_div"] - 1.96 * summary["sem_div"]
    summary["div_ci_upper"] = summary["mean_div"] + 1.96 * summary["sem_div"]
    summary["div_ci_lower"] = summary["div_ci_lower"].clip(lower=0)

    data = summary.copy()
    sigma = 2

    for scenario in ["CC", "No CC"]:

        scenario_data = data[
            (data['Scenario'] == scenario) &
            (data["level_5"] == metric)
        ].sort_values('Height')

        if not scenario_data.empty:

            smoothed_mean = gaussian_filter1d(scenario_data['mean_div'], sigma)
            smoothed_lower = gaussian_filter1d(scenario_data['div_ci_lower'], sigma)
            smoothed_upper = gaussian_filter1d(scenario_data['div_ci_upper'], sigma)

            color = colors.get(scenario, 'black')

            ax_bottom.plot(smoothed_mean, scenario_data['Height'],
                           color=color, linewidth=2)

            ax_bottom.fill_betweenx(scenario_data['Height'],
                                    smoothed_lower, smoothed_upper,
                                    color=color, alpha=0.2)

            ax_bottom.plot(scenario_data['mean_div'], scenario_data['Height'],
                           color=color, alpha=0.3, linestyle='--', linewidth=1)

    ax_bottom.set_title(str(year), fontsize=25, loc='center')
    ax_bottom.set_xlabel(metric if metric != "Shannon" else "Shannon index",
                         fontsize=25)

# Shared Y label for bottom-left plot
axes[2, 0].set_ylabel("Height (m)", fontsize=25)

for ax in axes.ravel():
    ax.grid(alpha=0.3)
    for spine in ax.spines.values():
        spine.set_visible(False)

handles, labels = axes[0, 0].get_legend_handles_labels()
label_map = {'CC': 'Climate change', 'No CC': 'Baseline'}
friendly_labels = [label_map.get(l, l) for l in labels]

fig.legend(handles, friendly_labels,
           loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.02))

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(DirectoryPlots / "vertical_div_ts_profile_v2.pdf")



