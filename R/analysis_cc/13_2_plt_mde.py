
import pandas as pd
from pathlib import Path
import matplotlib
from matplotlib.pyplot import ylabel
from scipy.stats import entropy
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter1d
from matplotlib.lines import Line2D
matplotlib.use("MacOSX")

from scipy.signal import savgol_filter

base_dir = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua")
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc")

# Laod MDE data
mde = pd.read_csv(base_dir / "climdata_era5_cmip6_1981-2100_ssp245" /
                  "a5_mde_richness.csv")
mde.rename(columns={"bin_mid": "Height",
                    "mean_richness": "mean_mde_richness",
                    "sem_richness": "sem_mde_richness"}, inplace=True)
mde["sem_mde_richness"] *= 1.96 # CI

# Load simulated vertical richness
div = pd.read_csv(base_dir / "a5_vertical_diversity_cc_vs_no_cc.csv")

# - Prepare data for plotting

colors = {'CC': '#f7766e', 'No CC': '#004aad'}
year = 2050

div2030 = div[(div['Year'] == year)]

# # Identify the combinations to keep (Richness >= 2)
# richness = div2030[div2030['level_5'] == 'Richness']
# valid_combinations = richness[richness['Count'] >= 2][['SpeciesPool', 'Scenario']]
#
# # Merge back to keep only valid rows
# div_filtered = div2030.merge(valid_combinations, on=['SpeciesPool', 'Scenario'], how='inner')

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

sigma = 2  # smoothing strength

scenario_data = data[
    (data['Scenario'] == "CC") &
    (data["level_5"] == "Richness")
].sort_values('Height')

# Merge with MDE data
scenario_data = pd.merge(scenario_data, mde, on="Height")

# Gaussian smoothing
smoothed_mean = gaussian_filter1d(scenario_data['mean_div'], sigma=sigma)
smoothed_lower = gaussian_filter1d(scenario_data['div_ci_lower'], sigma=sigma)
smoothed_upper = gaussian_filter1d(scenario_data['div_ci_upper'], sigma=sigma)

color = colors.get("CC", 'black')

plt.close("all")

fig, ax = plt.subplots(figsize=(7, 6), sharey=True)

# Plot (smoothed) MoDVE richness
ax.plot(smoothed_mean, scenario_data['Height'],
        color="#f7766e", linewidth=2, label=f"MoDVE")
ax.fill_betweenx(scenario_data['Height'],
                 smoothed_lower, smoothed_upper,
                 color="#f7766e", alpha=0.2)
ax.plot(scenario_data['mean_div'], scenario_data['Height'],
        color="#f7766e", alpha=0.3, linestyle='--', linewidth=1)

# Plot MDE richness
ax.plot(scenario_data["mean_mde_richness"], scenario_data['Height'],
        color="#88A2AA", linewidth=2, label=f"MDE")
ax.fill_betweenx(scenario_data['Height'],
                 scenario_data["mean_mde_richness"] - scenario_data["sem_mde_richness"],
                 scenario_data["mean_mde_richness"] + scenario_data["sem_mde_richness"],
                 color="#88A2AA", alpha=0.2)

# Axis labels, grid, title
ax.set_xlabel("Richness", fontsize=30)
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
ax.grid(alpha=0.3)
for spine in ax.spines.values():
    spine.set_visible(False)

# Shared Y label only on left-most subplot
ax.set_ylabel("Height (m)", fontsize=30)

# --- One global legend ---
#plt.rcParams.update({'font.size': 22})
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', ncol=1, frameon=False,
           fontsize=22) #, bbox_to_anchor=(0.5, 1.0))

plt.tight_layout()  # leave space for global legend
plt.show()

fig.savefig(DirectoryPlots / f"mde_vs_modve_vertical_richness_smooth_2050.pdf")