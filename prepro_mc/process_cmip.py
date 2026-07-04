import pandas as pd
from glob import glob
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

data_directory = "/Users/johanna/Uni/masterarbeit/data/mc_input/climate/cmip6/piControl/tas"
nc_files = glob(f"{data_directory}/baf*.nc")

combined_data_lst = []

for path in nc_files:

    model = path.split('/')[-1].split("_")[3]

    ds = xr.open_dataset(path, engine="h5netcdf")

    tas_data = ds['tas'].sel(lat=ds.lat.values[0], lon=ds.lon.values[0])  # shape: (time,)
    tas_values = tas_data.values  # numpy array of shape (time,)
    time_values = ds['time'].values  # numpy datetime64 array

    # Create DataFrame
    df = pd.DataFrame({
        'time': time_values,
        'tas': tas_values,
        'model': model,
        'lat': ds.lat.values[0],
        'lon': ds.lon.values[0]
    })

    combined_data_lst.append(df)

combined_data = pd.concat(combined_data_lst, ignore_index=True)

# Extract the year from the time column
combined_data['year'] = combined_data['time'].apply(lambda date: date.year)

# Check number of years per model
print(combined_data.groupby(['model', "year"]).first().reset_index().groupby('model').size())

# Check number of lat/lon pairs per model
print(combined_data.groupby(['model', 'lat', 'lon']).size().reset_index(name='count'))

# Get unique years per model and map to 1-500
model_year_map = (
    combined_data[['model', 'year']]
    .drop_duplicates()
    .sort_values(['model', 'year'])
    .groupby('model')
    .cumcount() + 1
)

# Assign this normalized year to a new column
combined_data['normalized_year'] = (
    combined_data.set_index(['model', 'year'])
      .index
      .map(dict(zip(combined_data[['model', 'year']]
                    .drop_duplicates()
                    .apply(tuple, axis=1),
                   model_year_map)))
)

# Convert tas from Kelvin to Celsius
combined_data['temp'] = combined_data['tas'] - 273.15

# Plot time series with hue based on model
plt.figure(figsize=(12, 6))
sns.lineplot(data=combined_data, x='normalized_year', y='temp', hue='model')
