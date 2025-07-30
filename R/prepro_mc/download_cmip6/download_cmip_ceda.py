import numpy as np
import xarray as xr
import requests
import pandas as pd

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

token = "eyJhbGciOiJSUzI1NiIsInR5cCIgOiAiSldUIiwia2lkIiA6ICI4ZjhmaUpyaUtDY3hmaHhzdU5vazVEekdJdFZ4amhhTWNJa05ZX2U4MnhJIn0.eyJleHAiOjE3NTQxMzk4MTAsImlhdCI6MTc1Mzg4MDYxMCwianRpIjoiM2YwYzU4ZjAtM2VhMy00MzRiLTg2NDItMGRjMTEzZjQ4MGFhIiwiaXNzIjoiaHR0cHM6Ly9hY2NvdW50cy5jZWRhLmFjLnVrL3JlYWxtcy9jZWRhIiwic3ViIjoiMmQxZTEzNGYtN2FjNi00YWQ4LTg1YjQtNmMyMmNkNzg0Y2RmIiwidHlwIjoiQmVhcmVyIiwiYXpwIjoic2VydmljZXMtcG9ydGFsLWNlZGEtYWMtdWsiLCJzZXNzaW9uX3N0YXRlIjoiMzU0MDhlY2UtNmRlYy00ZWY4LWI5ODUtNDJlYWU0NTJiMDYyIiwiYWNyIjoiMSIsInNjb3BlIjoiZW1haWwgb3BlbmlkIHByb2ZpbGUgZ3JvdXBfbWVtYmVyc2hpcCIsInNpZCI6IjM1NDA4ZWNlLTZkZWMtNGVmOC1iOTg1LTQyZWFlNDUyYjA2MiIsImVtYWlsX3ZlcmlmaWVkIjp0cnVlLCJuYW1lIjoiSm9oYW5uYSBUcm9zdCIsInByZWZlcnJlZF91c2VybmFtZSI6Imp0MjgxIiwibG9jYWxlIjoiZW4iLCJnaXZlbl9uYW1lIjoiSm9oYW5uYSIsImZhbWlseV9uYW1lIjoiVHJvc3QiLCJlbWFpbCI6ImpvaGFubmEudHJvc3QuMTk5N0BnbWFpbC5jb20ifQ.E0S9wMOQ0MAmXp3U8e9FYHz4KgudUa7peGBs3yLysQDi5iWia5G6bvRT8heeFaO0rrmL2D_-WnKSHzVQEP2nusyZjF2D18b1koN1BB1zAlEGtq7cl68T9jQ89gDIZm9C5OJqHPs3vpldkiW441Y9Ai0bgJxvVdjWJ8S-ZoDojDESg7DQBWj3XER8UWZrlAdZGhUGrRyNSXidVVbcysk-XYF1XtTIbsFD4pNd-Z8Z3CSxNMhcn-Akl5jbOHUnbHIzOBwk1QLTj5nJBBasB5gTaO8xK1W0ZDv_0pY9hz9lGIGa6WKiP6uDFb5mp2zJzgqj3iFKJywggwVCUgvydpOW7Q"

# Configure the area of interest and dates
area = {"lat": [-42.9, -42.2], "lon": [-22.8, -22.2]}
start = "2024-01-01"
end = "2024-12-30"

# Precipitation: GFDL-ESM4, CMCC_CM2-SR5, BCC-CSM2-MR, GISS- E2-1-G, INM-CM4-8 (LPB region)
# Temperatire: FIO-ESM-2-0, CanESM5 and HadGEM3-GC31-MM -> none of them available in CEDA data (https://doi.org/10.5285/c107618f1db34801bb88a1e927b82317, https://doi.org/10.1038/s41597-023-02528-x)

models = {
    "temperature": ["CMCC-ESM2", "CMCC-CM2-SR5", "MRI-ESM2-0"],
    "precipitation": ["IPSL-CM6A-LR", "INM-CM4-8", "GFDL-ESM4"],
}
    # "HadGEM3-GC31-MM", -> Not available in CEDA
    # "CMCC-ESM2", # Best for temp
    # "CMCC-CM2-SR5", # Best for temp, and precip in LPB region
    # "MRI-ESM2-0", # Best for temp
    # BCC-ESM1 -> not available in CEDA
    # ACCESS-ESM1-5, -> not available in CEDA
    # "IPSL-CM6A-LR", # Best for precipitation
    # IPSL-CM6A-LR-INCA, -> not available in CEDA
    # "INM-CM4-8  ", # Best for precipitation in LPB
    #"GFDL-ESM4" # Best for precipitation in LPB

gregorian_lst = [False, True, ]

# Extract correct indices for latitude and longitude
lats = np.linspace(-90, 90, 580)
lons = np.linspace(-180, 180, 1421)

lat_idx_range = [np.argmin(np.abs(lat - lats)).astype(int) for lat in area["lat"]]
lon_idx_range = [np.argmin(np.abs(lon - lons)).astype(int) for lon in area["lon"]]

# - Extract time indices
if gregorian:
    date_range = pd.date_range(start="2015-01-01", end="2100-12-31", freq="D")
    dates = date_range.strftime("%Y-%m-%d").to_numpy()

years = range(2015, 2101)  # Inclusive
months = range(1, 13)
days = range(1, 31)  # 30 days per month
dates = np.asarray(
    [f"{year:04d}-{month:02d}-{day:02d}"
     for year in years
     for month in months
     for day in days]
)
date_idx_range = [np.where(dates == start)[0][0], np.where(dates == end)[0][0]]

subset_str = f"[{date_idx_range[0]}:1:{date_idx_range[1]}][{lat_idx_range[0]}:1:{lat_idx_range[1]}][{lon_idx_range[0]}:1:{lon_idx_range[1]}]"
url3 = f"dap2://dap.ceda.ac.uk/thredds/dodsC/badc/evoflood/data/Downscaled_CMIP6_Climate_Data/hurs/SSP245/Global_hurs_Downscaled_{model}_2015-2100_ssp245_compressed.nc?hurs{subset_str}"

print(subset_str)

session = TokenSession(token)
ds = open_opendap_xarray(url3, token)



full_url = "https://dap.ceda.ac.uk/thredds/dodsC/badc/evoflood/data/Downscaled_CMIP6_Climate_Data/hurs/SSP245/Global_hurs_Downscaled_HadGEM3-GC31-LL_2015-2100_ssp245_compressed.nc?hurs"
session = TokenSession(token)
ds = open_opendap_xarray(full_url, token)