import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("MacOSX")
from pathlib import Path

# ------ 1. Load and process CMIP6 data ------

start = "2024-01-01"
end = "2024-12-30"
out_dir = Path("../../data/mc_input/climate/cmip6_ceda")

ds_tas = xr.open_dataset(out_dir / f"baf_tas_day_ssp245_{start}_{end}.nc")
ds_hurs = xr.open_dataset(out_dir / f"baf_hurs_day_ssp245_{start}_{end}.nc")

# ------ 2. Process the downloaded data ------

# Load the data and create dataframe with model and values
df_hurs = ds_hurs.to_dataframe().reset_index()
df_tas = ds_tas.to_dataframe().reset_index()

# Group by model and date and calculate mean temperature
df_tas['time'] = pd.to_datetime(df_tas['time'].astype(str))

# Add avaerage of all models to df

# Average temperature across all models
df_model_avg = df_tas.groupby(['time', 'lat', 'lon'], as_index=False)['tas'].mean()
df_spatial_avg = df_model_avg.groupby('time', as_index=False)['tas'].mean()
df_spatial_avg["model"] = "CMIP6 Average"

df_tas = pd.concat(
    [df_tas,
     df_spatial_avg]
)

# ------ 3. Load ERA5 ------

in_dir = Path("/Users/johanna/Uni/masterarbeit/data/mc_input/regua")
climdata_reg = pd.read_csv(in_dir / "era5_climdata_2024_v2.csv")

# Get daily means
climdata_reg['obs_time'] = pd.to_datetime(climdata_reg['obs_time'], format='mixed')
era5_daily = climdata_reg.set_index('obs_time').resample('D').mean(numeric_only=True).reset_index()
era5_daily.rename(columns={"temp": "tas",
                           "obs_time": "time"}, inplace=True)
era5_daily["model"] = "ERA5"
df_tas_all = pd.concat([df_tas, era5_daily[["time", "model", "tas"]]], ignore_index=True)

# Plot the time series of the mean temperature for each model

plt.figure(figsize=(12, 6))
for model in df_tas_all['model'].unique():
    model_data = df_tas_all[df_tas_all['model'] == model]
    plt.plot(model_data['time'], model_data['tas'], label=model, alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig("../../figs/mc_input/compare_climate_regua_tas_cmip6_era5.png", dpi=700)

print(df_tas_all.groupby("model")['tas'].mean())


# -- Same for hurs

era5_daily.rename(columns={"relhum": "hurs"}, inplace=True)
df_hurs_all = pd.concat([df_hurs, era5_daily[["time", "model", "hurs"]]], ignore_index=True)

# Plot the time series of the mean temperature for each model

plt.figure(figsize=(12, 6))
for model in df_hurs_all['model'].unique():
    model_data = df_hurs_all[df_hurs_all['model'] == model]
    plt.plot(model_data['time'], model_data['hurs'], label=model, alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig("../../figs/mc_input/compare_climate_regua_hurs_cmip6_era5.png", dpi=700)

print(df_hurs_all.groupby("model")['hurs'].mean())
