import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use("MacOSX")  # Use the appropriate backend for your environment

poolid = 1
sampling_method = "rand"

env_path = 'R/tests/data/output_a2/UniqueEnvVarsComb_t131-t198.csv'
niche_path = f'R/tests/data/output_a2_100spec_{sampling_method}_v2/SpeciesPool{poolid}.csv'

out_path = f'../../figs/a2_plots_test/niches_vs_env_{sampling_method}_pool{poolid}_v3.png'

env = pd.read_csv(env_path)
niches = pd.read_csv(niche_path)

niches.sort_values(by="OptimumLight", inplace=True)
niches["SpeciesID"] = niches.SpeciesID.astype(str)

env_stats = env.describe()
env_ranges = env_stats.loc[["min", "mean", "max"], :]
env_ranges['Light'] = env_ranges['Light'] * 900

# Assume `niches` is your DataFrame from Out[18]
# Example: niches = pd.read_csv('your_niches_data.csv')

# Set up figure
fig, axes = plt.subplots(4, figsize=(12, 10), sharex=True)
axes = axes.flatten()

variables = ['Light', 'Hum', 'Temp', 'Wind']

for i, var in enumerate(variables):

    print(f"{var}, - {env_ranges[var]['min']}, {env_ranges[var]['max']}")

    ax = axes[i]

    # Extract min, opt, max columns for the variable
    min_col = f'Min{var}'
    opt_col = f'Optimum{var}'
    max_col = f'Max{var}'

    # Melt the dataframe for seaborn
    df_plot = niches[['SpeciesID', min_col, opt_col, max_col]].copy()
    df_plot = df_plot.rename(columns={
        min_col: 'Min',
        opt_col: 'Opt',
        max_col: 'Max'
    })

    # Add vertical lines for Min-Max ranges
    for species_id in df_plot['SpeciesID']:
        min_val = df_plot.loc[df_plot['SpeciesID'] == species_id, 'Min'].values[
            0]
        max_val = df_plot.loc[df_plot['SpeciesID'] == species_id, 'Max'].values[
            0]
        ax.vlines(x=species_id, ymin=min_val, ymax=max_val, color='lightblue',
                  linewidth=1)
    sns.scatterplot(data=df_plot, x='SpeciesID', y='Opt',
                    label = 'Optimum' if i == 0 else "", marker=".", s=70,
                    color="coral", ax=ax, legend=(i == 0))

    # Add horizontal lines for env min and max
    ax.axhline(env_ranges[var]['min'], color='black', linestyle='--',
               label='Env Min' if i == 0 else "")
    ax.axhline(env_ranges[var]['max'], color='black', linestyle='--',
               label='Env Max' if i == 0 else "")
    ax.axhline(env_ranges[var]['mean'], color='black', linestyle=':',
               label='Env Mean' if i == 0 else "")

    lower = min(env_ranges[var]['min'], min(df_plot["Min"]))
    upper = max(env_ranges[var]['max'], max(df_plot["Max"]))
    if var == 'Light':
        ax.set_ylim([lower - 1, upper + 20])
    else:
        ax.set_ylim([lower - 1, upper + 1])

    ax.set_xlabel('SpeciesID')
    ax.set_ylabel(var)
    ax.tick_params(axis='x', labelrotation=90)

# Shared legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right')

plt.tight_layout(rect=[0, 0, 0.95, 1])
plt.show()

plt.savefig(out_path, dpi=700)

plt.close("all")