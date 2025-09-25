import pandas as pd
import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("MacOSX")

scenarios = ["climdata_era5_cmip6_1981-2100_ssp245_no_cc", "climdata_era5_cmip6_1981-2100_ssp245"]
base_dir = Path("/Users/johanna/Uni/masterarbeit/data/modve_output/regua") # /a5/forest0/
DirectoryPlots = Path("../../figs/a5_plots_test/cc_vs_no_cc")

DirectoryPlots.mkdir(exist_ok=True)

forests = np.arange(3)
species_pools = np.arange(1, 11)
rep = 1
ts_start = 100
ts_end = 199

vars = ['timeStep', 'Recruits', 'MortalityBranchFall', 'MortalityLight', 'MortalityCompetition',
        'MortalityNatural', 'MortalityHum', 'MortalityTemp']

community = None
for scenario in scenarios:
    for forest in forests:
        for sp in species_pools:
            path = base_dir / scenario / "a5" / f"forest{forest}" / f"ID_SpeciesP_{sp}_Rep_{rep}" / "CommunitySummary.csv"
            curr_community = pd.read_csv(path, usecols=vars)
            curr_community["Scenario"] = "CC" if scenario == "climdata_era5_cmip6_1981-2100_ssp245" else "No CC"
            curr_community["ForestID"] = forest
            curr_community["SpeciesPool"] = sp

            community = pd.concat((curr_community, community)) if community is not None else curr_community

# Map time steps from 80-199 to 1981-2100
community["Year"] = community["timeStep"] + 1901

# --- Plot on forest level for each scenario

# Columns to plot
cols_to_plot = ['Recruits', 'MortalityBranchFall', 'MortalityLight',
                'MortalityCompetition', 'MortalityNatural', 'MortalityHum',
                'MortalityTemp']

# Scenarios to plot
scenarios_str = ['CC', 'No CC']

# Set font size
plt.rcParams.update({'font.size': 14})

plt_shape = (3, 3)

for scenario in scenarios_str:
    df_scenario = community[community['Scenario'] == scenario]

    fig, axes = plt.subplots(*plt_shape, figsize=(16, 9), sharex=True)

    for i, col in enumerate(cols_to_plot):
        r, c = np.unravel_index(i, plt_shape)
        ax = axes[r, c]
        for forest in df_scenario['ForestID'].unique():
            df_forest = df_scenario[df_scenario['ForestID'] == forest]

            # Group by timeStep and compute mean & std across SpeciesPool
            stats = df_forest.groupby('Year')[col].agg(
                ['mean', 'std']).reset_index()

            ax.plot(stats['Year'], stats['mean'], label=f'Forest {forest}')
            ax.fill_between(stats['Year'],
                            stats['mean'] - stats['std'],
                            stats['mean'] + stats['std'],
                            alpha=0.2)

        ax.set_ylabel(col)

    ax.set_xlabel('Year')
    ax.legend()
    axes[-1, -1].axis('off')
    axes[-1, -2].axis('off')
    plt.tight_layout()
    plt.show()
    plt.savefig(DirectoryPlots / f"demography_forests_{scenario}.pdf")

# --- Plot on species pool level

for scenario in scenarios_str:
    df_scenario = community[community['Scenario'] == scenario]

    fig, axes = plt.subplots(*plt_shape, figsize=(16, 9), sharex=True)

    for i, col in enumerate(cols_to_plot):
        r, c = np.unravel_index(i, plt_shape)
        ax = axes[r, c]
        for sp in np.sort(df_scenario['SpeciesPool'].unique()):
            df_forest = df_scenario[df_scenario['SpeciesPool'] == sp]

            # Group by Year and compute mean & std across SpeciesPool
            stats = df_forest.groupby('Year')[col].agg(
                ['mean', 'std']).reset_index()

            ax.plot(stats['Year'], stats['mean'], label=f'Species pool {sp}')
            ax.fill_between(stats['Year'],
                            stats['mean'] - stats['std'],
                            stats['mean'] + stats['std'],
                            alpha=0.2)

        ax.set_ylabel(col)

    ax.set_xlabel('Year')
    axes[-1, -1].axis('off')
    axes[-1, -2].axis('off')

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower right', bbox_to_anchor=(0.6, 0.02))

    plt.tight_layout()
    plt.show()
    plt.savefig(DirectoryPlots / f"demography_species_pools_{scenario}.pdf")