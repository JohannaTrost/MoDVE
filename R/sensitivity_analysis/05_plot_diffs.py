import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import matplotlib.ticker as mticker

matplotlib.use("MacOSX")

DirectoryPlots = Path("/Users/johanna/Uni/masterarbeit/figs/sensitivity_analysis/")
DirectoryModelResults = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios")

# ----- Species distribution ----- #

fileName = "a5_Rep_1_climdata_era5_cmip6_1906-2024_ssp245_119ts_SpeciesVertical_mcGradients.csv"
verticalDistr = pd.read_csv(DirectoryModelResults / fileName)

filtered_distr = verticalDistr[(verticalDistr["timeStep"] >= 187) & (verticalDistr["timeStep"] <= 197)]

agg_distr = (
    filtered_distr
        .groupby(["speciesPool", "SpeciesID", "scenario"])
        .agg(
            mean_iqrHeight=("iqrHeight", "mean"),
            mean_meanHeight=("meanHeight", "mean")
        )
        .reset_index()
)

#steepness_map = {1.0: "1.0", 0.5: "0.5", 0.0: "0.0", 1.5: "1.5"}
steepness_map = {1.0: "1.0", 0.5: "0.5", 0.0: "0.0", 1.5: "1.5"}
agg_distr["Steepness (s)"] = agg_distr["scenario"].map(steepness_map)

# ----- Diversity data ----- #

agg_richness = pd.read_csv(DirectoryModelResults / "a5_vertical_richness_steepness.csv")

# Use names for scenarios instead of numbers
agg_richness["Steepness (s)"] = agg_richness["Scenario"].map(steepness_map)

cols = {
    "1.0": "#D4DDE0",  # move first
    "0.0":     "#E3E9EA",
    "0.5":  "#C5D2D5",
    "1.5":    '#88A2AA'
}
order = ["1.0", "0.0", "0.5", "1.5"]

# ----- Plotting ----- #

fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(11, 3), sharey=True)

sns.boxplot(
    data=agg_distr,
    x="mean_meanHeight",
    y="Steepness (s)",
    hue="Steepness (s)",
    orient='h',
    ax=axs[0],
    notch=True,
    palette=cols,
    medianprops=dict(linewidth=1.2, color='black'),
    flierprops=dict(marker='o', markeredgecolor='gray'),
    showmeans=True,
    meanprops=dict(
        marker='o',
        markerfacecolor='#F7766E',
        markeredgecolor='#F7766E',
        markersize=4
    ),
    order=order,
    hue_order=order,
    legend=False
)

sns.boxplot(
    data=agg_distr,
    x="mean_iqrHeight",
    y="Steepness (s)",
    hue="Steepness (s)",
    orient='h',
    ax=axs[1],
    notch=True,
    palette=cols,
    medianprops=dict(linewidth=1.2, color='black'),
    flierprops=dict(marker='o', markeredgecolor='gray'),
    showmeans=True,
    meanprops=dict(
        marker='o',
        markerfacecolor='#F7766E',
        markeredgecolor='#F7766E',
        markersize=4
    ),
    order=order,
    hue_order=order,
    legend=False
)

sns.boxplot(
    data=agg_richness,
    x="maxRichness",
    y="Steepness (s)",
    hue="Steepness (s)",
    orient='h',
    ax=axs[2],
    notch=True,
    palette=cols,
    medianprops=dict(linewidth=1.2, color='black'),
    flierprops=dict(marker='o', markeredgecolor='gray'),
    showmeans=True,
    meanprops=dict(
        marker='o',
        markerfacecolor='#F7766E',
        markeredgecolor='#F7766E',
        markersize=4
    ),
    order=order,
    hue_order=order,
    legend=False
)

sns.boxplot(
    data=agg_richness,
    x="maxHeight",
    y="Steepness (s)",
    hue="Steepness (s)",
    orient='h',
    ax=axs[3],
    notch=True,
    palette=cols,
    medianprops=dict(linewidth=1.2, color='black'),
    flierprops=dict(marker='o', markeredgecolor='gray'),
    showmeans=True,
    meanprops=dict(
        marker='o',
        markerfacecolor='#F7766E',
        markeredgecolor='#F7766E',
        markersize=4
    ),
    order=order,
    hue_order=order,
    legend=False
)

axs[0].set_xlabel("Species \nposition (m)", fontsize=18)
axs[1].set_xlabel("Species \nrange (IQR) (m)", fontsize=18)
axs[2].set_xlabel("Max. richness", fontsize=18)
axs[3].set_xlabel("Max. richness\n height (m)", fontsize=18)

for ax in axs:
    ax.set_ylabel("Steepness (s)", fontsize=18)
    ax.tick_params(axis='x', labelsize=13)
    ax.tick_params(axis='y', labelsize=13)
    ax.grid(alpha=0.3, axis='x')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

plt.tight_layout()
plt.savefig(DirectoryPlots / "spec_distr_richness_last10yrs_boxplot_v3.pdf")

plt.close("all")
