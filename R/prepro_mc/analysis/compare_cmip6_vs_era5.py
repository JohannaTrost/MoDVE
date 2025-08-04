import xarray as xr
from pyesgf.search import SearchConnection
import pandas as pd
import os
import requests
from tqdm import tqdm
import numpy as np
from pathlib import Path

# ------ 1. Download and process CMIP6 data ------

model_selection = ['CESM2-WACCM', 'GFDL-ESM4', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'INM-CM4-8']


def cfnoleap_to_datetime(da_or_ds):
    time_index = da_or_ds.indexes['time']

    try:
        # Try to convert from CFTimeIndex
        datetimeindex = time_index.to_datetimeindex()
    except AttributeError:
        # Already a DatetimeIndex
        datetimeindex = time_index

    # Convert to dataset if it's a DataArray
    if hasattr(da_or_ds, 'to_dataset'):
        ds = da_or_ds.to_dataset()
    else:
        ds = da_or_ds.copy()

    ds['time_dt'] = ('time', datetimeindex)
    ds = ds.swap_dims({'time': 'time_dt'})
    assert len(da_or_ds.time) == len(ds.time_dt)
    return ds


# Adapted from: https://stackoverflow.com/a/37573701
def download(url, filename):
    """
    Download a file hosted at <url> and write to <filename>
    """
    print("Downloading ", filename)
    r = requests.get(url, stream=True)
    total_size, block_size = int(r.headers.get('content-length', 0)), 1024
    with open(filename, 'wb') as f:
        for data in tqdm(r.iter_content(block_size),
                         total=total_size//block_size,
                         unit='KiB', unit_scale=True):
            f.write(data)

    if total_size != 0 and os.path.getsize(filename) != total_size:
        print("Downloaded size does not match expected size!\n",
              "FYI, the status code was ", r.status_code)


# Search for ssp245 data
conn = SearchConnection('http://esgf-data.dkrz.de/esg-search', distrib=True)
query = conn.new_context(project='CMIP6', experiment_id='ssp245',
                         variable_id="tas", frequency="day",
                         nominal_resolution="100 km",
                         source_id=model_selection, latest=True)

results = query.search()

print(len(results))

files = []
for i in range(len(results)):
    try:
        hit = results[i].file_context().search()
    except:
        hit = results[i].file_context().search()
    files += list(map(lambda f: {'filename': f.filename,
                                 'download_url': f.download_url,
                                 'opendap_url': f.opendap_url}, hit))
files = pd.DataFrame.from_dict(files)
len(files)

for fname in files['filename'].sort_values(): # print in alphabetical order
    print(fname)

# filter the DataFrame to drop duplicate filenames
files = files.drop_duplicates('filename')

# create a directory (inside the current working directory) to save the data to
data_directory = "/Users/johanna/Uni/masterarbeit/data/mc_input/climate/cmip6/ssp245/tas"
# only make the directory if it doesn't already exist
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

# Extract models and years
files['model'] = files['filename'].str.split('_').str[2]
files['first_year'] = files['filename'].apply(lambda x: int(x.split('_')[6].split('-')[0][:4]))
files['last_year'] = files['filename'].apply(lambda x: int(x.split('_')[6].split('-')[1].split('.')[0][:4]))

filtered_files = files[(files['first_year'] <= 2024) & (files['last_year'] >= 2024)].reset_index(drop=True)

# download the data, one file at a time
for i in range(1, len(filtered_files)):
    url = filtered_files['download_url'].loc[i]
    filename = filtered_files['filename'].loc[i]
    path_to_write = os.path.join(data_directory, filename)
    regional_filename = os.path.join(data_directory, "baf_" + filename)

    print(f"Downloading {filename},", i)

    # only download if the filtered_files doesn't already exist.
    if not os.path.exists(path_to_write) and not os.path.exists(regional_filename):
        try:
            download(url, path_to_write)
        except Exception as e:
            print(f"Error downloading {filename}: {e}")
            continue

    # Create new filename
    if not os.path.exists(regional_filename):
        # Load and subset with xarray
        try:
            ds = xr.open_dataset(path_to_write, engine="h5netcdf")
        except Exception as e:
            print(f"{filename} corrupted, retry download: {e}")

            try:
                download(url, path_to_write)
            except Exception as e:
                print(f"Error downloading {filename}: {e}")

            try:
                ds = xr.open_dataset(path_to_write, engine="h5netcdf")
            except Exception as e:
                print(f"{filename} corrupted, skip and remove file: {e}")
                os.remove(path_to_write)
                continue

        # Inspect lon range
        if ds.lon.min().item() >= 0:
            ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180))
            ds = ds.sortby('lon')  # Optional: to reorder the longitude axis

        # Subset the dataset to our area in atlantic brazilian forest
        regional_data = ds.sel(lat=slice(-42.9, -42.2), lon=slice(-22.8, -22.2))

        # Try point location in case the region is not available
        if regional_data.sizes['lat'] == 0 or regional_data.sizes['lon'] == 0:
            regional_data = ds.sel(lat=(-42.9 + -42.2 ) / 2,
                                   lon=(-22.8 + -22.2 ) / 2,
                                   method="nearest")

        # In case we have only a single grid celll left convert to a dataframe
        if np.prod(list(regional_data.sizes.values())) > 0:
            # Delete the old file
            os.remove(path_to_write)
            # Write regional data to netcdf file
            regional_data.to_netcdf(regional_filename)
        else:
            print(f"Skipping {path_to_write} as it does not contain the expected region.")
    else:
        print(f"Regional file {regional_filename} already exists, skipping.")

# ------ 2. Process the downloaded data ------

# Load the data and create dataframe with model and values
cmip6_files = list(Path(data_directory).glob("baf*.nc"))
data = []
for file in cmip6_files:
    ds = xr.open_dataset(file)
    # Convert to DataFrame
    ds = cfnoleap_to_datetime(ds)
    df = ds.to_dataframe().reset_index()
    if 'bnds' in df.columns:
        df = df[df['bnds'] == 0]
    elif 'nbnd' in df.columns:
        df = df[df['nbnd'] == 0]
    else:
        print(df.columns)
    df['model'] = file.stem.split('_')[3]  # Extract model name from filename

    df = df[['time', 'lat', 'lon', 'tas', 'model']]

    data.append(df)

cmip6 = pd.concat(data, ignore_index=True)

# Print number of models
print(f"Number of models: {cmip6['model'].nunique()}")

# Filter for year 2024
cmip6_2024 = cmip6[cmip6['time'].apply(lambda x: x.year == 2024)].reset_index(drop=True)

# Group by model and date and calculate mean temperature
cmip6_2024['time'] = pd.to_datetime(cmip6_2024['time'].astype(str))
cmip6_2024 = cmip6_2024.groupby(['model', cmip6_2024['time']])['tas'].mean().reset_index()
print(cmip6_2024['model'].value_counts())

# Remove leap day (Feb 29) from all models to standardize to 365 days
cmip6_2024 = cmip6_2024[~((cmip6_2024['time'].dt.month == 2) & (cmip6_2024['time'].dt.day == 29))]

# Convert from Kelvin to Celsius
cmip6_2024['temp'] = cmip6_2024['tas'] - 273.15

# Plot the time series of the mean temperature for each model
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("MacOSX")
plt.figure(figsize=(12, 6))
for model in cmip6_2024['model'].unique():
    model_data = cmip6_2024[cmip6_2024['model'] == model]
    plt.plot(model_data['time'], model_data['temp'], label=model)

# ------ 3. Load ERA5 ------

in_dir = Path("/Users/johanna/Uni/masterarbeit/data/mc_input/pirineus")
climdata_reg = pd.read_csv(in_dir / "era5_climdata_2024.csv")

# Filter data from 20.9.2024 to 23.9.2024
climdata_reg['obs_time'] = pd.to_datetime(climdata_reg['obs_time'], format='mixed')
start_date = pd.Timestamp('2024-09-20 00:00:00')
end_date = pd.Timestamp('2024-09-23 23:59:59')
filtered_data = climdata_reg[(climdata_reg['obs_time'] >= start_date) & (climdata_reg['obs_time'] <= end_date)]

# - Compute daily mean values
filtered_data.set_index('obs_time', inplace=True)

# Resample by day and calculate mean
daily_means = filtered_data.resample('D').mean(numeric_only=True)
daily_means = daily_means.reset_index()

# Display result
print(daily_means)
