import xarray as xr
from pyesgf.search import SearchConnection
import pandas as pd
import os
import requests
from tqdm import tqdm
import numpy as np

# Create ESGF search connection
conn = SearchConnection('http://esgf-data.dkrz.de/esg-search', distrib=True)

model_data = []

for exp in ['piControl', 'historical', 'ssp245']:
    for var in ['tas', 'hurs']:
        for res in ['100 km', '50 km']:
            print(f"Searching: var={var}, res={res}")
            ctx = conn.new_context(
                project='CMIP6',
                experiment_id='piControl',
                variable_id=var,
                frequency='day',
                nominal_resolution=res,
                latest=True
            )

            try:
                results = ctx.search()
                for result in results:
                    # dataset_id looks like: 'CMIP6.CMIP.MPI-M.MPI-ESM1-2-HR.piControl.day.gr.v20200205'
                    parts = result.dataset_id.split('.')
                    if len(parts) >= 4:
                        model = parts[3]  # source_id is the 4th element
                        model_data.append({
                            'model': model,
                            'var': var,
                            'res': res,
                            'experiment': exp,
                        })
            except Exception as e:
                print(f"Error with var={var}, res={res}: {e}")

# Convert to DataFrame
df_models = pd.DataFrame(model_data).drop_duplicates()

print(df_models)

# Step 1: Group by (var, experiment) and get sets of models for each group
grouped_models = df_models.groupby(['var', 'experiment'])['model'].apply(set)

# Step 2: Compute intersection across all (var, experiment) combinations
common_models = set.intersection(*grouped_models)

# Result: models common to all (var, experiment) groups
print("Models present in all (var, experiment) combinations:")
print(common_models)

# Check intersection of models with the below list
model_names = [  # doi: 10.1088/1748-9326/abd7ad
    "BCC-CSM2-MR",
    "CanESM5",
    "CESM2-WACCM",
    "CNRM-CM6-1",
    "CNRM-ESM2-1",
    "FGOALS-g3",
    "GFDL-ESM4",
    "INM-CM4-8",
    "INM-CM5-0",
    "IPSL-CM6A-LR",
    "MIROC6",
    "MPI-ESM1-2-HR",
    "MPI-ESM1-2-LR",
    "MRI-ESM2-0",
    "UKESM1-0-LL"
]
model_selection = list(common_models.intersection(set(model_names)))
model_selection = ['CESM2-WACCM', 'GFDL-ESM4', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'INM-CM4-8']

print(df_models[df_models["model"].isin(model_selection)])

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


# Search for piControl data
conn = SearchConnection('http://esgf-data.dkrz.de/esg-search', distrib=True)
query = conn.new_context(project='CMIP6', experiment_id='piControl',
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
data_directory = "/Users/johanna/Uni/masterarbeit/data/mc_input/climate/cmip6/piControl/tas_v2"
# only make the directory if it doesn't already exist
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

# Extract models and years
files['model'] = files['filename'].str.split('_').str[2]
files['first_year'] = files['filename'].apply(lambda x: int(x.split('_')[6].split('-')[0][:4]))
files['last_year'] = files['filename'].apply(lambda x: int(x.split('_')[6].split('-')[1].split('.')[0][:4]))

# Only keep those with more than 10 years of data
print(max(files['last_year'] - files['first_year'])) # 49
filtered_files = files[files['last_year'] - files['first_year'] >= 10]

filtered_files.reset_index(drop=True, inplace=True)

print(filtered_files)

# download the data, one file at a time
for i in range(len(filtered_files)):
    url = filtered_files['download_url'].loc[i]
    filename = filtered_files['filename'].loc[i]
    path_to_write = os.path.join(data_directory, filename)
    regional_filename = os.path.join(data_directory, "baf_" + filename)

    print(i)

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


# if there is an error downloading a file, you'll need to either delete the file,
# or change this block of code to call download(url, path_to_write) without
# first checking if path_to_write exists,




