import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
import scipy.stats as stats

matplotlib.use("MacOSX")

base_dir = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/pirineus/scenarios") # /a5/forest0/
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc")

senarios = {
    "climdata_era5_cmip6_1906-2024_ssp245_119ts_no_mc_grad" : 0,
    "climdata_era5_cmip6_1906-2024_ssp245_119ts_0.5_mc_grad" : 0.5,
    "climdata_era5_cmip6_1906-2024_ssp245_119ts" : 1,
    "climdata_era5_cmip6_1906-2024_ssp245_119ts_1.5_mc_grad" : 1.5
}


DirectoryPlots.mkdir(exist_ok=True)

# --- Load data --- #
forests = np.arange(3)
species_pools = np.arange(1, 11)
rep = 1
ts_start = 80
ts_end = 197

vars = ["SpeciesID", "IndividualID", 'Status', 'Mass', 'Z', 'X', 'Y']

species_distr = None
for scenario, s in senarios.items():
        for sp in species_pools:
            print(f"Scenario: {scenario}, SpeciesPool: {sp}")
            for ts in range(ts_start, ts_end + 1):
                path = (base_dir / scenario / "a5" / f"spec200_rand" /
                        f"ID_SpeciesP_{sp}_Rep_1" / f"IndividualMatrixTimeStep{ts}.csv")

                curr_species_distr = pd.read_csv(path, usecols=vars)
                curr_species_distr["Scenario"] = s
                curr_species_distr["SpeciesPool"] = sp
                curr_species_distr["TimeStep"] = ts

                species_distr = pd.concat((curr_species_distr, species_distr)) if species_distr is not None else curr_species_distr

# Map time steps from 80-197 to 1906-2024
species_distr["Year"] = species_distr["TimeStep"] + 1826
species_distr["Height"] = species_distr["Z"] - 0.5
# Filter individuals that are alive
species_distr = species_distr[species_distr["Status"] == 1]

species_distr.to_csv(base_dir / "a5_species_distribution_steepness.csv", index=False)

# --- Compute richness

species_distr = pd.read_csv(base_dir / "a5_species_distribution_steepness.csv")
species_counts = species_distr.groupby(["Scenario", "Year", "SpeciesPool",
                                        "SpeciesID", "Height"]).size().reset_index(name='Count')

# - Overall richness
# Filter for years 2013 to 2023
filtered = species_counts[(species_counts['Year'] >= 2013) & (species_counts['Year'] <= 2023)]

# Group by Scenario and SpeciesPool, and count unique SpeciesID
unique_species_counts = (
    filtered.groupby(['Scenario', 'SpeciesPool'])['SpeciesID']
    .nunique()
    .reset_index(name='UniqueSpecies')
)
unique_species_counts

# Assuming `unique_species_counts` is already computed as before

# Group by Scenario to compute mean and SEM across SpeciesPools
summary_stats = unique_species_counts.groupby('Scenario')['UniqueSpecies'].agg(
    Average='mean',
    SEM=lambda x: stats.sem(x, nan_policy='omit')  # SEM, ignoring NaNs if any
).reset_index()

summary_stats

# per height
richness = species_counts.groupby(["Scenario", "Year", "SpeciesPool", "Height"])['Count'].apply(lambda count: len(count > 0)).reset_index()

# Filter last 10 time steps with richness >= 1
filtered_richness = richness[(richness["Year"] >= 2013) & (richness["Year"] <= 2023)]

# Aggregate across last ten time steps
agg_richness = (
    filtered_richness.loc[
        filtered_richness.groupby(['Scenario', 'Year', 'SpeciesPool'])['Count'].idxmax(),
        ['Scenario', 'Year', 'SpeciesPool', 'Count', 'Height']
    ]
    .reset_index(drop=True)
).rename(columns={
    "Count": "maxRichness",
    "Height": "maxHeight"
})

agg_richness.to_csv(base_dir / "a5_vertical_richness_steepness.csv", index=False)
