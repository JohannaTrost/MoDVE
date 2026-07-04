import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker

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

# - Species position for specific year
df_plot = species_distr[species_distr["Year"] == 2030]



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
            DirectoryPlots / f"Species_height_diff_ts_Forest{forest}_SP{sp}.pdf")

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

# --- Plot position shift on species level:

DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc/position_shift")
species_distr = pd.read_csv(base_dir / "a5_species_distribution_cc_vs_no_cc.csv")

# Assuming species_distr is a pandas DataFrame
species_distr_stats = (
    species_distr
    .groupby(["Scenario", "SpeciesPool", "SpeciesID", "ForestID", "TimeStep", "Year"], dropna=False)
    .agg(
        Position=("Height", lambda x: np.nanmean(x)),
        Mass=("Mass", lambda x: np.nanmean(x)),
        IQR=("Height", lambda x: np.nanpercentile(x.dropna(), 75) - np.nanpercentile(x.dropna(), 25) if x.notna().any() else np.nan),
        Range=("Height", lambda x: np.nanmax(x) - np.nanmin(x) if x.notna().any() else np.nan)
    )
    .reset_index()
)

# Ensure Scenario has categorical order like factor(levels = c("No CC", "CC"))
species_distr_stats["Scenario"] = pd.Categorical(
    species_distr_stats["Scenario"], categories=["No CC", "CC"], ordered=True
)

# Arrange (sort) by Scenario
species_distr_stats = species_distr_stats.sort_values("Scenario")

# Select specific columns (like dplyr::select)
species_distr_stats = species_distr_stats[
    ["Scenario", "SpeciesPool", "SpeciesID", "ForestID", "TimeStep", "Year", "Position"]
]

# Pivot wider (similar to tidyr::pivot_wider)
species_distr_stats = species_distr_stats.pivot(
    index=["SpeciesPool", "SpeciesID", "ForestID", "TimeStep", "Year"],
    columns="Scenario",
    values="Position"
).reset_index()

# Calculate diff = `CC` - `No CC`
species_distr_stats["diff"] = species_distr_stats["CC"] - species_distr_stats["No CC"]
species_distr_stats.dropna(subset = ["diff"], inplace = True)

species_shift = species_distr_stats[species_distr_stats["Year"] >= 2080]
species_shift = species_shift.groupby(["SpeciesPool", "SpeciesID"], dropna=False).agg(
        AvgDiff=("diff", lambda x: np.nanmean(x)),
        SdDiff=("diff", lambda x: np.nanstd(x) / np.sqrt(len(x)) if len(x) > 0 else np.nan)
).reset_index()

# -- Plot the shift per species

# Sort by AvgDiff
df_sorted = species_shift.sort_values("AvgDiff")

# Compute error margins
error = 1.96 * df_sorted["SdDiff"]

# Plot
plt.figure(figsize=(8, 5))

plt.errorbar(
    x=range(len(df_sorted)),
    y=df_sorted["AvgDiff"],
    yerr=error,
    fmt='o',
    color='#886d7a',         # point color
    ecolor='#886d7a',        # error bar color
    elinewidth=1,            # thinner error bars
    capsize=4,
    markersize=4
)

# Labels
plt.xlabel("Species", fontsize=14)
plt.ylabel("Species shift with CC (m)", fontsize=14)

# Increase tick font sizes
plt.xticks([], fontsize=12)
plt.yticks(fontsize=12)

# Add horizontal line at y=0
plt.axhline(0, color='gray', linestyle='--', linewidth=1)

# Remove frame (spines)
for spine in plt.gca().spines.values():
    spine.set_visible(False)

# Grid and layout
plt.grid(True, linestyle='-', alpha=0.5)
plt.tight_layout()
plt.show()
plt.savefig(DirectoryPlots / "species_shift.png", dpi=700)

# --- Overall trend of species shift -- #

plt.rcParams.update({'font.size': 20})
colors = {'CC': '#f7766e', 'No CC': '#004aad'}  # keep same colors as your example

# --- 1) Aggregate per SpeciesID / SpeciesPool / ForestID / Year / Scenario ---
# This computes for each species (within each pool/forest/year/scenario) the mean height and the within-species range
per_species = (
    species_distr
    .groupby(['SpeciesID', 'SpeciesPool', 'ForestID', 'Year', 'Scenario'], observed=True)
    .agg(
        mean_height = ('Height', 'mean'),
        range_height = ('Height', lambda x: x.max() - x.min()),
        n_individuals = ('Height', 'count')
    )
    .reset_index()
)

# --- 2) Aggregate across Species to get time series (mean ± 95% CI) for position and range ---
# For each Year and Scenario we compute the mean of mean_height across species, sem, and 95% CI
def summarize_across_species(df, value_col, prefix):
    summary = (
        df
        .groupby(['Scenario', 'Year'], observed=True)[value_col]
        .agg(['mean', 'std', 'count'])
        .rename(columns={'mean': f'mean_{prefix}', 'std': f'std_{prefix}', 'count': f'n_{prefix}'})
        .reset_index()
    )
    # SEM and 95% CI (approx normal, z=1.96)
    summary[f'sem_{prefix}'] = summary[f'std_{prefix}'] / np.sqrt(summary[f'n_{prefix}'])
    summary[f'ci_lower_{prefix}'] = summary[f'mean_{prefix}'] - 1.96 * summary[f'sem_{prefix}']
    summary[f'ci_upper_{prefix}'] = summary[f'mean_{prefix}'] + 1.96 * summary[f'sem_{prefix}']
    return summary

summary_pos = summarize_across_species(per_species, 'mean_height', 'position')
summary_range = summarize_across_species(per_species, 'range_height', 'range')

# For plotting convenience, merge the two summaries on Scenario+Year so we can pick xticks from range of years
# But plotting will iterate over each independently (like your example)
years_min = min(summary_pos['Year'].min(), summary_range['Year'].min())
years_max = max(summary_pos['Year'].max(), summary_range['Year'].max())

# --- 3) Plotting - make two side-by-side subplots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharex=True)

# ---- Position subplot (left) ----
for scenario, group in summary_pos.groupby('Scenario'):
    ax1.plot(group['Year'], group['mean_position'],
             color=colors.get(scenario, None), label=scenario, linewidth=2)
    ax1.fill_between(group['Year'],
                     group['ci_lower_position'],
                     group['ci_upper_position'],
                     color=colors.get(scenario, None), alpha=0.2)

ax1.set_xlabel('')
ax1.set_ylabel('Species position (m)', fontsize=30)
xticks = range(years_min - 1, years_max + 1, 10)
ax1.set_xticks(xticks)
ax1.set_xticklabels(xticks, rotation=25)
ax1.set_xlim(years_min, years_max)
ax1.grid(alpha=0.3)
for spine in ax1.spines.values():
    spine.set_visible(False)

# ---- Range subplot (right) ----
for scenario, group in summary_range.groupby('Scenario'):
    ax2.plot(group['Year'], group['mean_range'],
             color=colors.get(scenario, None), label=scenario, linewidth=2)
    ax2.fill_between(group['Year'],
                     group['ci_lower_range'],
                     group['ci_upper_range'],
                     color=colors.get(scenario, None), alpha=0.2)

ax2.set_xlabel('')
ax2.set_ylabel('Species range (m)', fontsize=30)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xticks, rotation=25)
ax2.set_xlim(years_min, years_max)
ax2.grid(alpha=0.3)
for spine in ax2.spines.values():
    spine.set_visible(False)

# Force integers
ax1.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
ax2.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# Legend (one legend centered below)
handles, labels = ax1.get_legend_handles_labels()
# If you want friendly labels instead of 'CC'/'No CC' adapt here:
label_map = {'CC': 'Climate change', 'No CC': 'Baseline'}
friendly_labels = [label_map.get(l, l) for l in labels]
fig.legend(handles, friendly_labels, title=None,
           loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0))

plt.tight_layout(rect=[0, 0.12, 1, 1])

outpath = DirectoryPlots / "position_range_ts.pdf"
plt.savefig(outpath)

# --- Compute position/ range shifts --- #

# ===========================================================
# 1. POSITION SHIFT  (already done before, but included for completeness)
# ===========================================================

# Pivot for position
pos_pivot = (
    per_species.pivot_table(
        index=["SpeciesID", "SpeciesPool", "ForestID", "Year"],
        columns="Scenario",
        values="mean_height",
        observed=True
    )
    .reset_index()
)

# Keep only species present in both scenarios
pos_pivot = pos_pivot.dropna(subset=["CC", "No CC"])

# Compute position shift
pos_pivot["position_shift"] = pos_pivot["CC"] - pos_pivot["No CC"]

position_shift = pos_pivot[
    ["SpeciesID", "SpeciesPool", "ForestID", "Year",
     "No CC", "CC", "position_shift"]
].rename(columns={"No CC": "position_NoCC", "CC": "position_CC"})

# ----- Relative change in position ---- #
pos_pivot["rel_position_shift"] = (
    (pos_pivot["CC"] - pos_pivot["No CC"]) / pos_pivot["No CC"]
) * 100
# Average relatie shiftin 2030
rel_shift_2030 = pos_pivot[pos_pivot["Year"] >= 2030
]["rel_position_shift"].mean()

# ===========================================================
# 2. RANGE SHIFT
# ===========================================================

# Pivot for range
range_pivot = (
    per_species.pivot_table(
        index=["SpeciesID", "SpeciesPool", "ForestID", "Year"],
        columns="Scenario",
        values="range_height",
        observed=True
    )
    .reset_index()
)

# Keep only species present in both scenarios
range_pivot = range_pivot.dropna(subset=["CC", "No CC"])

# Compute range shift
range_pivot["range_shift"] = range_pivot["CC"] - range_pivot["No CC"]

range_shift = range_pivot[
    ["SpeciesID", "SpeciesPool", "ForestID", "Year",
     "No CC", "CC", "range_shift"]
].rename(columns={"No CC": "range_NoCC", "CC": "range_CC"})


# ===========================================================
# 3. AVERAGE SHIFT PER YEAR (mean, sem, CI)
# ===========================================================

def summarize_shift(df, value_col, prefix):
    out = (
        df.groupby("Year")[value_col]
        .agg(["mean", "std", "count"])
        .rename(columns={"mean": f"mean_{prefix}", "std": f"std_{prefix}", "count": f"n_{prefix}"})
        .reset_index()
    )

    out[f"sem_{prefix}"] = out[f"std_{prefix}"] / np.sqrt(out[f"n_{prefix}"])
    out[f"ci_lower_{prefix}"] = out[f"mean_{prefix}"] - 1.96 * out[f"sem_{prefix}"]
    out[f"ci_upper_{prefix}"] = out[f"mean_{prefix}"] + 1.96 * out[f"sem_{prefix}"]

    return out

summary_pos_shift = summarize_shift(position_shift, "position_shift", "pos_shift")
summary_range_shift = summarize_shift(range_shift, "range_shift", "range_shift")

# How many species are moving up between 2080 and 2100?
shift_subset = position_shift[
    (position_shift["Year"] >= 2080) & (position_shift["Year"] <= 2100)
]

# Count the number of species moving up each year
shift_counts = (
    shift_subset
    .groupby("Year")
    .agg(
        n_species_up=("position_shift", lambda x: (x > 0).sum()),
        n_species_down=("position_shift", lambda x: (x < 0).sum())
    )
    .reset_index()
)

perc_up = shift_counts["n_species_up"] / (shift_counts["n_species_up"] + shift_counts["n_species_down"]) * 100
print(round(perc_up.mean(), 1))
print(round(perc_up.sem() * 1.96, 1))  # 95% CI