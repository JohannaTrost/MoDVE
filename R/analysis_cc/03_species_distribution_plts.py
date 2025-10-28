import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import alpha

matplotlib.use("MacOSX")

scenarios = ["climdata_era5_cmip6_1981-2100_ssp245_no_cc", "climdata_era5_cmip6_1981-2100_ssp245"]
base_dir = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua") # /a5/forest0/
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc")

DirectoryPlots.mkdir(exist_ok=True)

# --- Load data --- #
forests = np.arange(3)
species_pools = np.arange(1, 11)
rep = 1
ts_start = 80
ts_end = 199

vars = ["SpeciesID", "IndividualID", 'Status', 'Mass', 'Z', 'X', 'Y']

species_distr = None
for scenario in scenarios:
    for forest in forests:
        for sp in species_pools:
            print(f"Scenario: {scenario}, Forest: {forest}, SpeciesPool: {sp}")
            for ts in range(ts_start, ts_end + 1):
                path = (base_dir / scenario / "a5" / f"forest{forest}" /
                        f"ID_SpeciesP_{sp}_Rep_{rep}" / f"IndividualMatrixTimeStep{ts}.csv")
                curr_species_distr = pd.read_csv(path, usecols=vars)
                curr_species_distr["Scenario"] = "CC" if scenario == "climdata_era5_cmip6_1981-2100_ssp245" else "No CC"
                curr_species_distr["ForestID"] = forest
                curr_species_distr["SpeciesPool"] = sp
                curr_species_distr["TimeStep"] = ts

                species_distr = pd.concat((curr_species_distr, species_distr)) if species_distr is not None else curr_species_distr

# Map time steps from 80-199 to 1981-2100
species_distr["Year"] = species_distr["TimeStep"] + 1901
species_distr["Height"] = species_distr["Z"] - 0.5
# Filter individuals that are alive
species_distr = species_distr[species_distr["Status"] == 1]

species_distr.to_csv(base_dir / "a5_species_distribution_cc_vs_no_cc.csv", index=False)

species_distr = pd.read_csv(base_dir / "a5_species_distribution_cc_vs_no_cc.csv")

# --- Plot position and range for each replicate (sp, forest) --- #

# For each Forest and species pool, make a boxplot of the height of each species in each scenario (exclude all years before 2020)
for forest in forests:
    for sp in species_pools:
        df_plot = species_distr[(species_distr["ForestID"] == forest) &
                                (species_distr["SpeciesPool"] == sp) &
                                (species_distr["Year"] >= 2020)]
        # Take average height of individual
        df_plot = df_plot.groupby(
            ["Scenario", "ForestID", "SpeciesPool", "SpeciesID", "IndividualID"]
        )["Height"].mean().reset_index()

        # Sort the DataFrame by avg species height
        df_plot["AvgHeight"] = df_plot.groupby(["SpeciesID"])["Height"].transform(
            "mean")
        df_plot = df_plot.sort_values(by="AvgHeight").reset_index(drop=True)

        # First, calculate species order based on AvgHeight
        species_order = (
            df_plot.groupby("SpeciesID")["AvgHeight"]
            .mean()
            .sort_values()
            .index
        )

        # Create the boxplot
        plt.figure(figsize=(9, 3))
        sns.boxplot(
            data=df_plot,
            x="SpeciesID",
            y="Height",
            hue="Scenario",
            order=species_order,
            medianprops=dict(linewidth=1.2),
            flierprops=dict(marker='o', markeredgecolor='gray'),
        )

        # Style
        plt.xticks(rotation=45)
        plt.xlabel("Species ID")
        plt.ylabel("Height (m)")
        plt.legend(title="Scenario")
        plt.tight_layout()
        plt.show()
        plt.savefig(DirectoryPlots / f"Species_Height_Distribution_Forest{forest}_SP{sp}.png")

# --- Plot position and range comparing CC and no CC --- #

# 1. Overall last 10 years (2090-2100)
df_plot = species_distr[species_distr["Year"] >= 2089]
# Take average height of individual
df_plot = df_plot.groupby(
    ["Scenario", "ForestID", "SpeciesPool", "SpeciesID", "IndividualID"]
)["Height"].mean().reset_index()

# 1.1
plt.figure(figsize=(4, 4))
sns.violinplot(
    data=df_plot,
    x="Scenario",
    y="Height",
    hue="Scenario"
)
plt.tight_layout()
plt.savefig(DirectoryPlots / f"cc_vs_no_cc_position_violin.png")
# 1.2
plt.figure(figsize=(4, 4))
sns.boxplot(
    data=df_plot,
    x="Scenario",
    y="Height",
    hue="Scenario",
    notch=True,
    medianprops=dict(linewidth=1.2),
    flierprops=dict(marker='o', markeredgecolor='gray')
)
plt.tight_layout()
plt.savefig(DirectoryPlots / f"cc_vs_no_cc_position_box.png")
# 1.3 Species-wise position
df_plot_sp = df_plot.groupby(
    ["Scenario", "ForestID", "SpeciesPool", "SpeciesID"]
)["Height"].mean().reset_index()
plt.figure(figsize=(4, 4))
sns.boxplot(
    data=df_plot_sp,
    x="Scenario",
    y="Height",
    hue="Scenario",
    notch=True,
    medianprops=dict(linewidth=1.2),
    flierprops=dict(marker='o', markeredgecolor='gray')
)
plt.tight_layout()
plt.savefig(DirectoryPlots / f"cc_vs_no_cc_species_position_box.png")

plt.figure(figsize=(4, 4))
sns.violinplot(
    data=df_plot_sp,
    x="Scenario",
    y="Height",
    hue="Scenario"
)
plt.tight_layout()
plt.savefig(DirectoryPlots / f"cc_vs_no_cc_species_position_violin.png")

# 1.4 Species-wise IQR
df_plot_sp_iqr = df_plot.groupby(
    ["Scenario", "ForestID", "SpeciesPool", "SpeciesID"]
)["Height"].agg(lambda x: np.percentile(x, 75) - np.percentile(x, 25)).reset_index(name="IQR")

plt.figure(figsize=(4, 4))
sns.boxplot(
    data=df_plot_sp_iqr,
    x="Scenario",
    y="IQR",
    hue="Scenario",
    notch=True,
    medianprops=dict(linewidth=1.2),
    flierprops=dict(marker='o', markeredgecolor='gray')
)
plt.tight_layout()
plt.savefig(DirectoryPlots / f"cc_vs_no_cc_species_range_box.png")

plt.figure(figsize=(5, 5))
sns.violinplot(
    data=df_plot_sp_iqr,
    x="Scenario",
    y="IQR",
    hue="Scenario"
)
plt.tight_layout()
plt.savefig(DirectoryPlots / f"cc_vs_no_cc_species_range_violin.png")

# --- 2. Grouped by forests for last 10 years (2090-2100)

# - 2.1 Individual position
df_plot = species_distr[species_distr["Year"] >= 2089]
# Take average height of individual
df_plot = df_plot.groupby(
    ["Scenario", "ForestID", "SpeciesPool", "SpeciesID", "IndividualID"]
)["Height"].mean().reset_index()

plt.figure(figsize=(4, 4))
sns.violinplot(
    data=df_plot,
    x="ForestID",
    y="Height",
    hue="Scenario"
)
plt.legend(ncol=2, title="Scenario")
plt.tight_layout()
plt.savefig(DirectoryPlots / f"forest_cc_vs_no_cc_position_violin.pdf")
plt.figure(figsize=(4, 4))
sns.boxplot(
    data=df_plot,
    x="ForestID",
    y="Height",
    hue="Scenario",
    notch=True,
    medianprops=dict(linewidth=1.2),
    flierprops=dict(marker='o', markeredgecolor='gray')
)
plt.legend(ncol=2, title="Scenario", loc='upper right')
plt.tight_layout()
plt.savefig(DirectoryPlots / f"forest_cc_vs_no_cc_position_box.pdf")

# - 2.2 Species-wise position

plt.figure(figsize=(4, 4))
sns.violinplot(
    data=df_plot_sp,
    x="ForestID",
    y="Height",
    hue="Scenario"
)
plt.legend(ncol=2, title="Scenario")
plt.tight_layout()
plt.savefig(DirectoryPlots / f"forest_cc_vs_no_cc_species_position_violin.pdf")
plt.figure(figsize=(4, 4))
sns.boxplot(
    data=df_plot_sp,
    x="ForestID",
    y="Height",
    hue="Scenario",
    notch=True,
    medianprops=dict(linewidth=1.2),
    flierprops=dict(marker='o', markeredgecolor='gray')
)
plt.legend(ncol=2, title="Scenario", loc='upper right')
plt.tight_layout()
plt.savefig(DirectoryPlots / f"forest_cc_vs_no_cc_species_position_box_last20yrs.pdf")

# - 2.3 Species-wise IQR

plt.figure(figsize=(4, 4))
sns.violinplot(
    data=df_plot_sp_iqr,
    x="ForestID",
    y="IQR",
    hue="Scenario"
)
plt.legend(ncol=2, title="Scenario")
plt.tight_layout()
plt.savefig(DirectoryPlots / f"forest_cc_vs_no_cc_species_range_violin_last20yrs.pdf")
plt.figure(figsize=(4, 4))
sns.boxplot(
    data=df_plot_sp_iqr,
    x="ForestID",
    y="IQR",
    hue="Scenario",
    notch=True,
    medianprops=dict(linewidth=1.2),
    flierprops=dict(marker='o', markeredgecolor='gray')
)
plt.legend(ncol=2, title="Scenario", loc='upper right')
plt.tight_layout()
plt.savefig(DirectoryPlots / f"forest_cc_vs_no_cc_species_range_box_last20yrs.pdf")

# --- 3. Grouped by species pools and forest for last 10 years (2090-2100)

def plot_forest_species_position_violin(data, fname_save, y="Height"):
    # Get all unique species pools
    species_pools = data["SpeciesPool"].unique()

    # Create subplots: 5 rows, 2 columns
    fig, axes = plt.subplots(5, 2, figsize=(8, 12), sharey=True, sharex=True)
    axes = axes.flatten()  # Flatten to 1D array for easy looping

    # Loop through each species pool and create a violin plot
    for i, pool in enumerate(species_pools):
        ax = axes[i]
        subset = data[data["SpeciesPool"] == pool]

        sns.violinplot(
            data=subset,
            x="ForestID",
            y=y,
            hue="Scenario",
            ax=ax
        )

        ax.set_title(f"Species Pool: {pool}")
        ax.legend_.remove()  # Remove individual legends to avoid clutter

    # Put one legend for all plots at the bottom
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, title="Scenario", ncol=2, loc="lower center")

    plt.tight_layout(rect=[0, 0.05, 1, 1])  # Leave space at bottom for legend
    plt.savefig(DirectoryPlots / fname_save) #"sp_forest_cc_vs_no_cc_position_violin.pdf")
    plt.show()

def plot_forest_species_position_box(data, fname_save, y="Height"):

    # Create subplots: 5 rows, 2 columns
    fig, axes = plt.subplots(5, 2, figsize=(8, 12), sharey=True, sharex=True)
    axes = axes.flatten()  # Flatten to 1D array for easy looping

    # Loop through each species pool and create a violin plot
    for i, pool in enumerate(species_pools):
        ax = axes[i]
        subset = data[data["SpeciesPool"] == pool]

        sns.boxplot(
            data=subset,
            x="ForestID",
            y=y,
            hue="Scenario",
            ax=ax,
            notch=True,
            medianprops=dict(linewidth=1.2),
            flierprops=dict(marker='o', markeredgecolor='gray')
        )

        ax.set_title(f"Species Pool: {pool}")
        ax.legend_.remove()  # Remove individual legends to avoid clutter

    # Put one legend for all plots at the bottom
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, title="Scenario", ncol=2, loc="lower center")

    plt.tight_layout(rect=[0, 0.05, 1, 1])  # Leave space at bottom for legend
    plt.savefig(DirectoryPlots / fname_save)
    plt.show()

def plot_forest_species_position_swarm(data, fname_save, y="Height"):

    # Create subplots: 5 rows, 2 columns
    fig, axes = plt.subplots(5, 2, figsize=(8, 12), sharey=True, sharex=True)
    axes = axes.flatten()  # Flatten to 1D array for easy looping

    # Loop through each species pool and create a violin plot
    for i, pool in enumerate(species_pools):
        ax = axes[i]
        subset = data[data["SpeciesPool"] == pool]

        sns.swarmplot(
            data=subset,
            x="ForestID",
            y=y,
            hue="Scenario",
            ax=ax
        )

        ax.set_title(f"Species Pool: {pool}")
        ax.legend_.remove()  # Remove individual legends to avoid clutter

    # Put one legend for all plots at the bottom
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, title="Scenario", ncol=2, loc="lower center")

    plt.tight_layout(rect=[0, 0.05, 1, 1])  # Leave space at bottom for legend
    plt.savefig(DirectoryPlots / fname_save)
    plt.show()


# - 3.1 Individual position
plot_forest_species_position_violin(df_plot, "sp_forest_cc_vs_no_cc_position_violin.pdf")
plot_forest_species_position_box(df_plot, "sp_forest_cc_vs_no_cc_position_box.pdf")

# - 3.2 Species-wise position
plot_forest_species_position_violin(df_plot_sp, "sp_forest_cc_vs_no_cc_species_position_violin.pdf")
plot_forest_species_position_box(df_plot_sp, "sp_forest_cc_vs_no_cc_species_position_box.pdf")
plot_forest_species_position_swarm(df_plot_sp, "sp_forest_cc_vs_no_cc_species_position_swarm.pdf")

# - 3.3 Species-wise IQR
plot_forest_species_position_violin(df_plot_sp_iqr, "sp_forest_cc_vs_no_cc_species_position_violin_last20yrs.pdf", "IQR")
plot_forest_species_position_box(df_plot_sp_iqr, "sp_forest_cc_vs_no_cc_species_position_box_last20yrs.pdf", "IQR")
plot_forest_species_position_swarm(df_plot_sp_iqr, "sp_forest_cc_vs_no_cc_species_position_swarm_last20yrs.pdf", "IQR")

# --- Plot time series of avg position and range for each replicate (sp, forest) --- #

# - 1. Position
# For each Forest and species pool, make a boxplot of the height of each species in each scenario (exclude all years before 2020)
for forest in forests:
    for sp in species_pools:
        df_plot = species_distr[(species_distr["ForestID"] == forest) &
                                (species_distr["SpeciesPool"] == sp)]
        # Take average height of individual
        df_plot = df_plot.groupby(
            ["Scenario", "ForestID", "SpeciesPool", "SpeciesID", "TimeStep"]
        )["Height"].mean().reset_index()

        df_summary = (
            df_plot.groupby(["Scenario", "TimeStep"])["Height"]
            .agg(["mean", "std"])
            .reset_index()
        )
        # Handle NaN std (if only one data point)
        df_summary.loc[np.isnan(df_summary["std"]), "std"] = 0
        df_summary["mean"] = np.asarray(df_summary["mean"]).astype(float)

        # Make the plot
        plt.figure(figsize=(10, 6))

        # Plot mean line with ribbon for each scenario
        for scenario, df_scen in df_summary.groupby("Scenario"):
            sns.lineplot(
                data=df_scen,
                x="TimeStep",
                y="mean",
                label=scenario
            )
            plt.fill_between(
                df_scen["TimeStep"],
                df_scen["mean"] - df_scen["std"],
                df_scen["mean"] + df_scen["std"],
                alpha=0.2
            )

        plt.xlabel("TimeStep")
        plt.ylabel("Height (m)")
        plt.legend(title="Scenario")
        plt.tight_layout()
        plt.show()
        plt.savefig(DirectoryPlots / f"Species_Position_ts_Forest{forest}_SP{sp}.pdf")

# 2. Range (IQR) TODO

# --- Plot time series of avg difference position and range for each replicate (sp, forest) --- #

# - 1. Position
scenarios_str = ["No CC", "CC"]
for forest in forests:
    for sp in species_pools:

        plt.close("all")

        df_plot = species_distr[(species_distr["ForestID"] == forest) &
                                (species_distr["SpeciesPool"] == sp)]
        # Take average height of individual
        df_plot = df_plot.groupby(
            ["Scenario", "ForestID", "SpeciesPool", "SpeciesID", "Year"]
        )["Height"].mean().reset_index()

        # Pivot the data so we have species heights per scenario per time step
        df_pivot = df_plot.pivot_table(
            index=["SpeciesID", "Year"],
            columns="Scenario",
            values="Height"
        ).reset_index()

        # Remove rows where either scenario is missing (species not present in both scenarios at this timestep)
        df_pivot = df_pivot.dropna(subset=scenarios_str)

        # Calculate difference: "No CC" - "CC"
        df_pivot["diff"] = df_pivot["CC"] - df_pivot["No CC"]

        # Now summarize: mean and std of the differences at each time step
        df_diff_summary = df_pivot.groupby("Year")["diff"].agg(
            ["mean", "std"]).reset_index()

        # Handle NaN std if needed (e.g., only one species present at that timestep)
        df_diff_summary["std"] = df_diff_summary["std"].fillna(0)
        df_diff_summary["mean"] = np.asarray(df_diff_summary["mean"]).astype(float)

        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(df_diff_summary["Year"], df_diff_summary["mean"],
                 label="Mean Difference", color="grey")
        plt.fill_between(
            df_diff_summary["Year"],
            df_diff_summary["mean"] - df_diff_summary["std"],
            df_diff_summary["mean"] + df_diff_summary["std"],
            alpha=0.2, color="grey"
        )
        plt.hlines(y=0, xmin=ts_start + 1901, xmax=ts_end + 1901,
                   color='darkgrey', linestyles='--', label='No Difference')
        plt.xlabel("Year")
        plt.ylabel("Absolute height difference (CC - No CC) (m)")
        plt.title("Mean ± SD of species height differences between scenarios")
        plt.legend()
        plt.tight_layout()
        plt.show()

        # Save figure
        plt.savefig(
            DirectoryPlots / f"Species_HeightDiff_ts_Forest{forest}_SP{sp}.pdf")

# Species-wise Difference mean forests
for sp in species_pools:

    plt.close("all")

    df_plot = species_distr[species_distr["SpeciesPool"] == sp]

    # Pivot the data so we have species heights per scenario per time step
    df_pivot = df_plot.pivot_table(
        index=["SpeciesID", "Year", "SpeciesPool", "ForestID"],
        columns="Scenario",
        values="Height"
    ).reset_index()

    # Remove rows where either scenario is missing (species not present in both scenarios at this timestep)
    df_pivot = df_pivot.dropna(subset=scenarios_str)

    # Calculate difference: "No CC" - "CC"
    df_pivot["diff"] = df_pivot["CC"] - df_pivot["No CC"]

    # Take average height of individual
    df_plot = df_pivot.groupby(
        ["SpeciesPool", "SpeciesID", "Year"]
    )["diff"].mean().reset_index()

    # Restrict to 2001–2100
    subset_range = df_plot[
        (df_plot["Year"] >= 2030) & (df_plot["Year"] <= 2100)]

    # Compute median diff per species
    medians = subset_range.groupby("SpeciesID")["diff"].median().reset_index(
        name="median_diff")

    # Split into negative and positive median groups
    neg_species = medians[medians["median_diff"] < 0]["SpeciesID"]
    pos_species = medians[medians["median_diff"] >= 0]["SpeciesID"]

    # Last year per species (for annotation)
    last_years = df_plot.groupby("SpeciesID")["Year"].max().reset_index()
    last_values = df_plot.merge(last_years, on=["SpeciesID", "Year"],
                                how="inner")

    # Build one shared palette for all species
    all_species = pd.concat([neg_species, pos_species]).unique()
    palette = sns.color_palette("colorblind", n_colors=len(all_species))

    # Map SpeciesID -> unique color
    color_map = dict(zip(all_species, palette))

    # Define plotting function
    def plot_species(ax, species_list, title):
        data = df_plot[df_plot["SpeciesID"].isin(species_list)]

        # Map species in this group to their unique color
        group_palette = {sp: color_map[sp] for sp in species_list}

        sns.lineplot(
            data=data,
            x="Year",
            y="diff",
            hue="SpeciesID",
            palette=group_palette,
            ax=ax,
            legend=False,
            alpha=0.6
        )

        ax.hlines(
            y=0,
            xmin=df_plot["Year"].min(),
            xmax=df_plot["Year"].max(),
            color="darkgrey",
            linestyles="--",
            alpha=0.6
        )
        ax.set_title(title)
        ax.set_ylabel("Species position shift (CC - No CC) (m)")

        # annotate each species at its last year
        for sp in species_list:
            sp_last = last_values[last_values["SpeciesID"] == sp].iloc[0]
            ax.text(
                sp_last["Year"] + 1,
                sp_last["diff"],
                str(sp),
                fontsize=14,
                va="center",
                color=color_map[sp]  # match annotation color with line
            )


    # --- Step 7: Make subplots
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    plot_species(
        axes[0],
        neg_species,
        "Downward shift with climate change (2030–2100)"
    )
    plot_species(
        axes[1],
        pos_species,
        "Upward shift with climate change (2030–2100)"
    )

    axes[1].set_xlabel("Year")
    plt.tight_layout()
    plt.show()

    # Save figure
    plt.savefig(
        DirectoryPlots / f"Species_HeightShift_ts_SP{sp}.pdf")

# - 2. Range (IQR)
for forest in forests:
    for sp in species_pools:

        plt.close("all")

        df_plot = species_distr[(species_distr["ForestID"] == forest) &
                                (species_distr["SpeciesPool"] == sp)]
        # iqr
        df_plot = df_plot.groupby(
            ["Scenario", "ForestID", "SpeciesPool", "SpeciesID", "Year"]
        )["Height"].agg(
            lambda x: np.percentile(x, 75) - np.percentile(x, 25)).reset_index(
            name="IQR")

        # Pivot the data so we have species heights per scenario per time step
        df_pivot = df_plot.pivot_table(
            index=["SpeciesID", "Year"],
            columns="Scenario",
            values="IQR"
        ).reset_index()

        # Remove rows where either scenario is missing (species not present in both scenarios at this timestep)
        df_pivot = df_pivot.dropna(subset=scenarios_str)

        # Calculate difference: "No CC" - "CC"
        df_pivot["diff"] = df_pivot["CC"] - df_pivot["No CC"]

        # Now summarize: mean and std of the differences at each time step
        df_diff_summary = df_pivot.groupby("Year")["diff"].agg(
            ["mean", "std"]).reset_index()

        # Handle NaN std if needed (e.g., only one species present at that timestep)
        df_diff_summary["std"] = df_diff_summary["std"].fillna(0)
        df_diff_summary["mean"] = np.asarray(df_diff_summary["mean"]).astype(float)

        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(df_diff_summary["Year"], df_diff_summary["mean"],
                 label="Mean Difference", color="grey")
        plt.fill_between(
            df_diff_summary["Year"],
            df_diff_summary["mean"] - df_diff_summary["std"],
            df_diff_summary["mean"] + df_diff_summary["std"],
            alpha=0.2, color="grey"
        )
        plt.hlines(y=0, xmin=ts_start + 1901, xmax=ts_end + 1901,
                   color='darkgrey', linestyles='--', label='No Difference')
        plt.xlabel("Year")
        plt.ylabel("Absolute range difference (CC - No CC) (m)")
        plt.title("Mean ± SD of species range differences between scenarios")
        plt.legend()
        plt.tight_layout()
        plt.show()

        # Save figure
        plt.savefig(
            DirectoryPlots / f"Species_RangeDiff_ts_Forest{forest}_SP{sp}.pdf")

# --- Plot position vs range for each replicate (sp, forest) and for CC vs No CC --- #

df_plot = species_distr[species_distr["Year"] >= 2020]
# Take average height of individual
df_plot = df_plot.groupby(
    ["Scenario", "ForestID", "SpeciesPool", "SpeciesID", "IndividualID"]
)["Height"].mean().reset_index()

# Species-wise position
df_plot_sp = df_plot.groupby(
    ["Scenario", "ForestID", "SpeciesPool", "SpeciesID"]
)["Height"].mean().reset_index()
# Species-wise IQR
df_plot_sp["IQR"] = df_plot.groupby(
            ["Scenario", "ForestID", "SpeciesPool", "SpeciesID"]
        )["Height"].agg(
            lambda x: np.percentile(x, 75) - np.percentile(x, 25)).reset_index(
            name="IQR")["IQR"]

# ALl reoplicates in one plot
plt.figure(figsize=(5, 5))
sns.scatterplot(data=df_plot_sp, x="Height", y="IQR", hue="Scenario")
plt.savefig(DirectoryPlots / f"Species_Position_vs_Range_scatter_all.pdf")

# Individual replicates
g = sns.relplot(
    data=df_plot_sp,
    x="Height",
    y="IQR",
    hue="Scenario",
    col="ForestID",
    row="SpeciesPool",
    kind="scatter",
    height=2.5,      # height of each facet in inches
    aspect=1.2,    # width = aspect * height
    facet_kws={"sharex": False, "sharey": False}  # optional, can adjust depending on scale
)
plt.tight_layout()
plt.savefig(DirectoryPlots / f"Species_Position_vs_Range_scatter_facet.pdf")


# ---------------- PLot Preds Pos ---------------- #


preds = pd.read_csv(base_dir / "a5_pred_pos_cc_vs_no_cc.csv")

g = sns.FacetGrid(
    preds,
    col="ForestID",
    row="SpeciesPool",
    hue="Scenario",
    margin_titles=True,
    height=3,
    aspect=1.5
)

# observed "Position" lines
g.map_dataframe(sns.lineplot, x="Year", y="Position", linestyle="--",
                alpha=0.6)

# predicted "fit" line
g.map_dataframe(sns.lineplot, x="Year", y="fit")

# CI ribbon (upr / lwr)
def add_ci(data, color, **kwargs):
    plt.fill_between(
        data["Year"], data["lwr"], data["upr"],
        color=color, alpha=0.2, **kwargs
    )

g.map_dataframe(add_ci)

g.add_legend()
plt.show()
plt.savefig(DirectoryPlots / "pred_position_facet.pdf")
