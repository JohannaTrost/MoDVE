import numpy as np
import xarray as xr
import requests
import pandas as pd
from pathlib import Path

class TokenSession(requests.Session):
    def __init__(self, token):
        super().__init__()
        self.headers.update({"Authorization": f"Bearer {token}"})


def open_opendap_xarray(url, token):
    # Create a requests session with the token
    session = TokenSession(token)

    # Open with Pydap (which supports OPeNDAP slicing)
    ds = xr.open_dataset(url, engine="pydap",
                         backend_kwargs={"session": session})
    return ds

token = "eyJhbGciOiJSUzI1NiIsInR5cCIgOiAiSldUIiwia2lkIiA6ICI4ZjhmaUpyaUtDY3hmaHhzdU5vazVEekdJdFZ4amhhTWNJa05ZX2U4MnhJIn0.eyJleHAiOjE3NTUyNDczNTksImlhdCI6MTc1NDk4ODE1OSwianRpIjoiMTY5OGMyM2QtYTNmOC00NjFjLTgzMzItMjNkMWViMGFmMTEzIiwiaXNzIjoiaHR0cHM6Ly9hY2NvdW50cy5jZWRhLmFjLnVrL3JlYWxtcy9jZWRhIiwic3ViIjoiMmQxZTEzNGYtN2FjNi00YWQ4LTg1YjQtNmMyMmNkNzg0Y2RmIiwidHlwIjoiQmVhcmVyIiwiYXpwIjoic2VydmljZXMtcG9ydGFsLWNlZGEtYWMtdWsiLCJzZXNzaW9uX3N0YXRlIjoiYmIxOTg4YzctYzdhZi00YTM5LWIxOTktZGY3MWZjMDhhY2ZkIiwiYWNyIjoiMSIsInNjb3BlIjoiZW1haWwgb3BlbmlkIHByb2ZpbGUgZ3JvdXBfbWVtYmVyc2hpcCIsInNpZCI6ImJiMTk4OGM3LWM3YWYtNGEzOS1iMTk5LWRmNzFmYzA4YWNmZCIsImVtYWlsX3ZlcmlmaWVkIjp0cnVlLCJuYW1lIjoiSm9oYW5uYSBUcm9zdCIsInByZWZlcnJlZF91c2VybmFtZSI6Imp0MjgxIiwibG9jYWxlIjoiZW4iLCJnaXZlbl9uYW1lIjoiSm9oYW5uYSIsImZhbWlseV9uYW1lIjoiVHJvc3QiLCJlbWFpbCI6ImpvaGFubmEudHJvc3QuMTk5N0BnbWFpbC5jb20ifQ.ZlVOLVtfiOGFVcTTlcQifpGJclJ7bDFMcAA3yBZavxwTjZD08fPvvLFpmqTwBYKdvQKpHUO2Z4fZZ5lNa5jsruk4Tv9ML1ObHRRHfqPy0_9IFwlA90Neip-xPN2uYaDh3v0jm-YqKUOUdsXkT2f2dxO6_1C4wTZ29VhsQ9BQZgeeqL7SKeLzVlHYzJv-4WHyff-eIerhULksO3SV_sWgwPvHYuOyFdagD7IWrpSkaC0ZFIU2-AJ6T4TDEYywrjGYGiS7OXkOYdYBt3jN8iPnVxtxG0aJ2tOLa3tC3CBtkbV7DEW3zZ3SSpEyBtT5HsV3AUn9n8zzuh325iW_716JCw"

out_dir = Path("../../data/mc_input/climate/cmip6_ceda")
download_range = range(2093, 2101)
#download_range = range(2042, 2050)

# Configure the area of interest and dates
area = {"lat": [-42.9, -42.2], "lon": [-22.8, -22.2]}

# Precipitation: GFDL-ESM4, CMCC_CM2-SR5, BCC-CSM2-MR, GISS-E2-1-G, INM-CM4-8 (LPB region)
# Temperatire: FIO-ESM-2-0, CanESM5 and HadGEM3-GC31-MM -> none of them available in CEDA data (https://doi.org/10.5285/c107618f1db34801bb88a1e927b82317, https://doi.org/10.1038/s41597-023-02528-x)

models = {
    "tas": ["CMCC-ESM2", "CMCC-CM2-SR5", "MRI-ESM2-0"],
    "hurs": ["IPSL-CM6A-LR", "INM-CM4-8", "GFDL-ESM4"],
    'pr': ["INM-CM4-8", "GFDL-ESM4", "CMCC-CM2-SR5", "BCC-CSM2-MR"], # IPSL-CM6A-LR has 31,382 dates -> don't know which calendar to use
    'sfcWind': ["GFDL-ESM4", "BCC-CSM2-MR", "INM-CM4-8"],
    'ps': ["GFDL-ESM4", "BCC-CSM2-MR", "IPSL-CM6A-LR", "CMCC-ESM2", "CMCC-CM2-SR5",
           "MRI-ESM2-0"]
}
    # "HadGEM3-GC31-MM", -> Not available in CEDA
    # "CMCC-ESM2", # Best for temp
    # "CMCC-CM2-SR5", # Best for temp, and precip in LPB region
    # "MRI-ESM2-0", # Best for temp
    # BCC-ESM1 -> not available in CEDA
    # ACCESS-ESM1-5, -> not available in CEDA
    # "IPSL-CM6A-LR", # Best for precipitation
    # IPSL-CM6A-LR-INCA, -> not available in CEDA ()
    # "INM-CM4-8", # Best for precipitation in LPB
    #"GFDL-ESM4" # Best for precipitation in LPB

# calendar categories gregorian, no_leap_days, 360_days
calendars = {
    "tas": ["no_leap_days", "no_leap_days", "gregorian"],
    "hurs": ["gregorian", "no_leap_days", "no_leap_days"],
    "pr": ["no_leap_days", "no_leap_days", "no_leap_days",
        "no_leap_days"],
    "sfcWind": ["no_leap_days", "no_leap_days", "no_leap_days"],
    "ps": ["no_leap_days", "no_leap_days", "gregorian", "no_leap_days",
           "no_leap_days", "gregorian"]
    # Also chose models with overall good performance in prep and temp over all of Brazil
}
micropoint_vars = {"tas": "temp", "hurs": "relhum", "pr": "precip",
                   "sfcWind": "windspeed", "ps": "pres"}

# - Extract time indices
date_range = pd.date_range(start="2015-01-01", end="2100-12-31", freq="D")

dates_gregorian = date_range.strftime("%Y-%m-%d").to_numpy()
filtered_range = date_range[
        ~((date_range.month == 2) & (date_range.day == 29))]
dates_no_leap_days = filtered_range.strftime("%Y-%m-%d").to_numpy()
years = range(2015, 2101)  # Inclusive
months = range(1, 13)
days = range(1, 31)  # 30 days per month
dates_360days = np.asarray(
    [f"{year:04d}-{month:02d}-{day:02d}"
     for year in years
     for month in months
     for day in days]
)

for year in download_range:
    print(year)
    start = f"{year}-01-01"
    end = f"{year}-12-31"

    date_idx_range = {
        "gregorian": [np.where(dates_gregorian == start)[0][0], np.where(dates_gregorian == end)[0][0]], # -> 31411
        "no_leap_days": [np.where(dates_no_leap_days == start)[0][0], np.where(dates_no_leap_days == end)[0][0]] # -> 31390
    }

    # Dates in date format
    start_date = pd.to_datetime(start)
    end_date = pd.to_datetime(end)
    datetime_range = pd.date_range(start=start_date, end=end_date, freq='D')

    # Extract correct indices for latitude and longitude
    lats = np.linspace(-90, 90, 580)
    lons = np.linspace(-180, 180, 1421)
    lat_idx_range = [np.argmin(np.abs(lat - lats)).astype(int) for lat in area["lat"]]
    lon_idx_range = [np.argmin(np.abs(lon - lons)).astype(int) for lon in area["lon"]]
    # Regional coordinates
    area_coords = {
        "lat": lats[lat_idx_range[0]: lat_idx_range[1] + 1],
        "lon": lons[lon_idx_range[0]: lon_idx_range[1] + 1]
    }

    cmip_ensemble = []

    for var in models.keys():
        print(var)

        output_path = out_dir / f"baf_{var}_day_ssp245_{year}.nc"
        if output_path.exists():
            models_ds = xr.open_dataset(output_path)
            print(f"{output_path} already exists.")
        else:
            ds_lst = []
            for i, model in enumerate(models[var]):

                calendar = calendars[var][i]
                subset_str = (f"[{date_idx_range[calendar][0]}:1:{date_idx_range[calendar][1]}]"
                              f"[{lat_idx_range[0]}:1:{lat_idx_range[1]}]"
                              f"[{lon_idx_range[0]}:1:{lon_idx_range[1]}]")
                url3 = f"dap2://dap.ceda.ac.uk/thredds/dodsC/badc/evoflood/data/Downscaled_CMIP6_Climate_Data/{var}/SSP245/Global_{var}_Downscaled_{model}_2015-2100_ssp245_compressed.nc?{var}{subset_str}"

                print(subset_str)

                session = TokenSession(token)
                ds = open_opendap_xarray(url3, token)

                # Linear interpolation over time for missing days to align calendars
                if calendar == "no_leap_days":
                    all_dates = pd.to_datetime(dates_no_leap_days)
                    target_dates = all_dates[
                        (all_dates >= start_date) & (all_dates <= end_date)]

                    assert len(target_dates) == ds["time"].shape[0], "Different number of dates than in ds."

                    # Assign these dates to the dataset's time coordinate
                    ds = ds.assign_coords(time=target_dates)
                    # Reindex and interpolate missing days
                    var_filled = (
                        ds[var]
                        .reindex(time=datetime_range)
                        .interpolate_na(dim='time', method='linear')
                        #.ffill(dim='time')  # Fill forward to cover the end TODO check if necessary
                        #.bfill(dim='time')  # Optional: fill backwards at the start
                    )
                    # Assign interpolated values back to a dataset (if needed)
                    ds = var_filled.to_dataset(name=var)
                else:
                    ds = ds.assign_coords(time=datetime_range)

                # Assign proper decimal degrees for longitude
                ds = ds.assign_coords(lat=area_coords["lat"])
                ds = ds.assign_coords(lon=area_coords["lon"])

                ds_lst.append(ds[var].expand_dims(model=[model]))

            # Concatenate along the new 'model' dimension and turn into a Dataset
            models_ds = xr.concat(ds_lst, dim='model').to_dataset(name=var)

            # Save the dataset to a netCDF file
            models_ds.to_netcdf(output_path, mode='w', format='NETCDF4')
            print(f"Saved {var} data to {output_path}")

        # Save aggregated data
        var_ensemble = (models_ds
                        .mean(dim="model")
                        .rename({var: micropoint_vars[var]}))

        cmip_ensemble.append(var_ensemble)

    cmip_ensemble_ds = xr.merge(cmip_ensemble)
    output_path = out_dir / f"baf_ensemble_day_ssp245_{year}.nc"
    cmip_ensemble_ds.to_netcdf(output_path, mode='w', format='NETCDF4')