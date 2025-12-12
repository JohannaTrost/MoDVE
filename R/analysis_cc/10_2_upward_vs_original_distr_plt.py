import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker
from scipy.stats import alpha

matplotlib.use("MacOSX")

scenarios = ["climdata_era5_cmip6_1981-2100_ssp245_no_cc", "climdata_era5_cmip6_1981-2100_ssp245"]
base_dir = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua") # /a5/forest0/
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc/position_shift")

DirectoryPlots.mkdir(exist_ok=True)

# ------------------ Original community

species_shift = pd.read_csv(base_dir / "a5_species_shift_cc_vs_no_cc_last20_yrs_past_2050.csv")
species_shift.rename(columns={"diff": "AvgDiff", "diff_sem": "SdDiff"}, inplace=True)
species_shift.drop(columns=["initial_diff", "last_year_alive"], inplace=True)
species_shift = species_shift[species_shift["AvgDiff"] >= species_shift["AvgDiff"].quantile(0.75)]

# ------------------ Upward community

# Process upward shift of species of upward species simulation
upward_shift_spec = pd.read_csv(base_dir / "a5_upward_shifted_species_distribution_cc_vs_no_cc.csv")

# Assuming species_distr is a pandas DataFrame
upward_shift_spec_stats = (
    upward_shift_spec
    .groupby(["Scenario", "SpeciesPool", "SpeciesID", "ForestID", "TimeStep", "Year"], dropna=False)
    .agg(
        Position=("Height", lambda x: np.nanmean(x)),
        Mass=("Mass", lambda x: np.nanmean(x)),
        IQR=("Height", lambda x: np.nanpercentile(x.dropna(), 75) - np.nanpercentile(x.dropna(), 20) if x.notna().any() else np.nan),
        Range=("Height", lambda x: np.nanmax(x) - np.nanmin(x) if x.notna().any() else np.nan)
    )
    .reset_index()
)

# Ensure Scenario has categorical order like factor(levels = c("No CC", "CC"))
upward_shift_spec_stats["Scenario"] = pd.Categorical(
    upward_shift_spec_stats["Scenario"], categories=["No CC", "CC"], ordered=True
)

# Arrange (sort) by Scenario
upward_shift_spec_stats = upward_shift_spec_stats.sort_values("Scenario")

# Select specific columns (like dplyr::select)
upward_shift_spec_stats = upward_shift_spec_stats[
    ["Scenario", "SpeciesPool", "SpeciesID", "ForestID", "TimeStep", "Year", "Position"]
]

# Pivot wider (similar to tidyr::pivot_wider)
upward_shift_spec_stats = upward_shift_spec_stats.pivot(
    index=["SpeciesPool", "SpeciesID", "ForestID", "TimeStep", "Year"],
    columns="Scenario",
    values="Position"
).reset_index()

# Calculate diff = `CC` - `No CC`
upward_shift_spec_stats["diff_upward"] = upward_shift_spec_stats["CC"] - upward_shift_spec_stats["No CC"]
upward_shift_spec_stats.dropna(subset = ["diff_upward"], inplace = True)

# Take last 20 years for each species
species_shift_up = (
    upward_shift_spec_stats
    .sort_values(["SpeciesPool", "SpeciesID", "Year"])
    .groupby(["SpeciesPool", "SpeciesID"], group_keys=False)
    .apply(lambda g: g.tail(20))
)

species_shift_up = species_shift_up.groupby(["SpeciesPool", "SpeciesID"], dropna=False).agg(
        AvgDiff=("diff_upward", lambda x: np.nanmean(x)),
        SdDiff=("diff_upward", lambda x: np.nanstd(x) / np.sqrt(len(x)) if len(x) > 0 else np.nan)
).reset_index()

# --- Merge with original simulation shifts

species_shift_up["Simulation"] = "Upward shifted species"
species_shift["Simulation"] = "Entire community"
species_shift_comb = pd.concat([species_shift_up, species_shift])

# --- Plot shifts

# Sort
species_order = species_shift.sort_values('AvgDiff', ascending=True)[['SpeciesPool', 'SpeciesID']]
species_order = [(row["SpeciesPool"], row["SpeciesID"]) for _, row in species_order.iterrows()]

# mapping from SpeciesID -> x position
pos_map = {species: i for i, species in enumerate(species_order)}

# Filter rows whose (SpeciesPool, SpeciesID) pair is in species_order
df_plot = species_shift_comb.copy()
df_plot['pair'] = list(zip(df_plot['SpeciesPool'], df_plot['SpeciesID']))
df_plot = df_plot[df_plot['pair'].isin(species_order)]
# Map x position
df_plot['x'] = df_plot['pair'].map(pos_map)

# compute error margins
df_plot['error'] = 1.96 * df_plot['SdDiff']

# get simulation groups (preserve original order of appearance)
sims = df_plot['Simulation'].unique()
palette = ["#214E34", "#88A2AA"]

plt.rcParams.update({'font.size': 18})
plt.figure(figsize=(35, 5))

dodge_amount = 0.15  # adjust as needed
offsets = [-dodge_amount, dodge_amount]  # for 2 simulations

for i, (sim, color) in enumerate(zip(sims, palette)):
    sub = df_plot[df_plot['Simulation'] == sim].sort_values('x')
    plt.errorbar(
        x=sub['x'] + offsets[i],  # apply offset
        y=sub['AvgDiff'],
        yerr=sub['error'],
        fmt='o',
        color=color,
        ecolor=color,
        elinewidth=3,
        capsize=5,
        markersize=7,
        label=sim
    )

# Labels
plt.xlabel("Species", fontsize=23)
plt.ylabel("Species shift with CC (m)", fontsize=23)

plt.xticks([])

plt.axhline(0, color='gray', linestyle='--', linewidth=1)

# Add vertical gridlines at each species position
species_positions = df_plot['x'].unique()
for x_pos in species_positions:
    plt.axvline(x_pos, color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

for spine in plt.gca().spines.values():
    spine.set_visible(False)
plt.grid(alpha=0.3, axis='y')  # Only horizontal gridlines from grid()

plt.legend(title="Simulation", fontsize=18, title_fontsize=20, loc='lower right')

plt.tight_layout()

# Save then show
plt.savefig(DirectoryPlots / "species_shift_comp_v2.pdf")
