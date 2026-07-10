#!/usr/bin/env python3
"""
CMIP6 Climate Data Downloader CLI

Downloads and processes CMIP6 climate data from CEDA for specified years and variables.
"""

import argparse
import sys
import numpy as np
import xarray as xr
import requests
import pandas as pd
from pathlib import Path


class TokenSession(requests.Session):
    """Custom requests session with bearer token authentication."""

    def __init__(self, token):
        super().__init__()
        self.headers.update({"Authorization": f"Bearer {token}"})


def open_opendap_xarray(url, token):
    """Open OPeNDAP dataset with xarray using token authentication."""
    session = TokenSession(token)
    ds = xr.open_dataset(url, engine="pydap",
                         backend_kwargs={"session": session})
    return ds


def create_date_arrays(first, last):
    """Create date arrays for different calendar types."""
    date_range = pd.date_range(start=f"{first}-01-01", end=f"{last}-12-31", freq="D")

    dates_gregorian = date_range.strftime("%Y-%m-%d").to_numpy()
    filtered_range = date_range[
        ~((date_range.month == 2) & (date_range.day == 29))]
    dates_no_leap_days = filtered_range.strftime("%Y-%m-%d").to_numpy()

    years = range(first, last)
    months = range(1, 13)
    days = range(1, 31)  # 30 days per month
    dates_360days = np.asarray([
        f"{year:04d}-{month:02d}-{day:02d}"
        for year in years
        for month in months
        for day in days
    ])

    return dates_gregorian, dates_no_leap_days, dates_360days


def get_model_configurations():
    """Return model and calendar configurations."""
    models = {
        "tas": ["CMCC-ESM2", "CMCC-CM2-SR5", "MRI-ESM2-0"],
        "hurs": ["IPSL-CM6A-LR", "INM-CM4-8", "GFDL-ESM4"],
        "pr": ["INM-CM4-8", "GFDL-ESM4", "CMCC-CM2-SR5", "BCC-CSM2-MR"],
        "sfcWind": ["GFDL-ESM4", "BCC-CSM2-MR", "INM-CM4-8"],
        "ps": ["GFDL-ESM4", "BCC-CSM2-MR", "IPSL-CM6A-LR", "CMCC-ESM2",
               "CMCC-CM2-SR5", "MRI-ESM2-0"]
    }

    calendars = {
        "tas": ["no_leap_days", "no_leap_days", "gregorian"],
        "hurs": ["gregorian", "no_leap_days", "no_leap_days"],
        "pr": ["no_leap_days", "no_leap_days", "no_leap_days", "no_leap_days"],
        "sfcWind": ["no_leap_days", "no_leap_days", "no_leap_days"],
        "ps": ["no_leap_days", "no_leap_days", "gregorian", "no_leap_days",
               "no_leap_days", "gregorian"]
    }

    micropoint_vars = {
        "tas": "temp",
        "hurs": "relhum",
        "pr": "precip",
        "sfcWind": "windspeed",
        "ps": "pres"
    }

    return models, calendars, micropoint_vars


def process_year(year, out_dir, token, models, calendars, micropoint_vars,
        dates_gregorian, dates_no_leap_days, area, scenario_path, scenario, default_first_year, default_last_year):
    """Process climate data for a single year."""
    print(f"Processing year: {year}")

    start = f"{year}-01-01"
    end = f"{year}-12-31"

    date_idx_range = {
        "gregorian": [
            np.where(dates_gregorian == start)[0][0],
            np.where(dates_gregorian == end)[0][0]
        ],
        "no_leap_days": [
            np.where(dates_no_leap_days == start)[0][0],
            np.where(dates_no_leap_days == end)[0][0]
        ]
    }

    # Create date range for the year
    start_date = pd.to_datetime(start)
    end_date = pd.to_datetime(end)
    datetime_range = pd.date_range(start=start_date, end=end_date, freq='D')

    # Extract correct indices for latitude and longitude
    lats = np.linspace(-90, 90, 580)
    lons = np.linspace(-180, 180, 1421)
    lat_idx_range = [np.argmin(np.abs(lat - lats)).astype(int) for lat in
                     area["lat"]]
    lon_idx_range = [np.argmin(np.abs(lon - lons)).astype(int) for lon in
                     area["lon"]]

    # Regional coordinates
    area_coords = {
        "lat": lats[lat_idx_range[0]: lat_idx_range[1] + 1],
        "lon": lons[lon_idx_range[0]: lon_idx_range[1] + 1]
    }

    cmip_ensemble = []

    for var in models.keys():
        print(f"  Processing variable: {var}")

        output_path = out_dir / f"baf_{var}_day_{scenario}_{year}.nc"

        if output_path.exists():
            print(f"    {output_path} already exists, loading...")
            models_ds = xr.open_dataset(output_path)
        else:
            ds_lst = []

            for i, model in enumerate(models[var]):
                calendar = calendars[var][i]
                subset_str = (
                    f"[{date_idx_range[calendar][0]}:1:{date_idx_range[calendar][1]}]"
                    f"[{lat_idx_range[0]}:1:{lat_idx_range[1]}]"
                    f"[{lon_idx_range[0]}:1:{lon_idx_range[1]}]"
                )

                url = (
                    f"dap2://dap.ceda.ac.uk/thredds/dodsC/badc/evoflood/data/"
                    f"Downscaled_CMIP6_Climate_Data/{var}/{scenario_path}/"
                    f"Global_{var}_Downscaled_{model}_{default_first_year}-{default_last_year}_"
                    f"{'' if scenario_path == 'Historical' else 'ssp245_'}compressed.nc?"
                    f"{var}{subset_str}"
                )

                print(f"    Downloading {model} with subset: {subset_str}")

                ds = open_opendap_xarray(url, token)

                # Handle different calendar types
                if calendar == "no_leap_days":
                    all_dates = pd.to_datetime(dates_no_leap_days)
                    target_dates = all_dates[
                        (all_dates >= start_date) & (all_dates <= end_date)
                        ]

                    assert len(target_dates) == ds["time"].shape[0], \
                        "Different number of dates than in dataset."

                    ds = ds.assign_coords(time=target_dates)
                    var_filled = (
                        ds[var]
                        .reindex(time=datetime_range)
                        .interpolate_na(dim='time', method='linear')
                    )
                    ds = var_filled.to_dataset(name=var)
                else:
                    ds = ds.assign_coords(time=datetime_range)

                # Assign proper coordinates
                ds = ds.assign_coords(lat=area_coords["lat"])
                ds = ds.assign_coords(lon=area_coords["lon"])

                ds_lst.append(ds[var].expand_dims(model=[model]))

            # Concatenate models and save
            models_ds = xr.concat(ds_lst, dim='model').to_dataset(name=var)
            models_ds.to_netcdf(output_path, mode='w', format='NETCDF4')
            print(f"    Saved {var} data to {output_path}")

        # Create ensemble mean
        var_ensemble = (
            models_ds
            .mean(dim="model")
            .rename({var: micropoint_vars[var]})
        )
        cmip_ensemble.append(var_ensemble)

    # Save ensemble dataset
    cmip_ensemble_ds = xr.merge(cmip_ensemble)

    ensemble_output_path = out_dir / f"baf_ensemble_day_{scenario}_{year}.nc"
    cmip_ensemble_ds.to_netcdf(ensemble_output_path, mode='w', format='NETCDF4')
    print(f"  Saved ensemble data to {ensemble_output_path}")


def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(
        description="Download and process CMIP6 climate data from CEDA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --first-year 1981 --last-year 1985 --out-dir /path/to/output
  %(prog)s -f 2000 -l 2010 -o ./climate_data --token-file my_token.txt
        """
    )

    parser.add_argument(
        "-f", "--first-year",
        type=int,
        required=True,
        help="First year to process (e.g., 1981)"
    )

    parser.add_argument(
        "-l", "--last-year",
        type=int,
        required=True,
        help="Last year to process (e.g., 2015)"
    )

    parser.add_argument(
        "-o", "--out-dir",
        type=str,
        required=True,
        help="Output directory path"
    )

    parser.add_argument(
        "-s", "--scenario",
        type=str,
        choices=["historical", "ssp245"],
        required=True,
        help="Scenario to download: 'historical' or 'ssp245'"
    )

    parser.add_argument(
        "-t", "--token",
        type=str,
        help="CEDA authentication token"
    )

    parser.add_argument(
        "--token-file",
        type=str,
        help="File containing the CEDA authentication token"
    )

    parser.add_argument(
        "--area-lat-min",
        type=float,
        default=-42.9,
        help="Minimum latitude for area of interest (default: -42.9)"
    )

    parser.add_argument(
        "--area-lat-max",
        type=float,
        default=-42.2,
        help="Maximum latitude for area of interest (default: -42.2)"
    )

    parser.add_argument(
        "--area-lon-min",
        type=float,
        default=-22.8,
        help="Minimum longitude for area of interest (default: -22.8)"
    )

    parser.add_argument(
        "--area-lon-max",
        type=float,
        default=-22.2,
        help="Maximum longitude for area of interest (default: -22.2)"
    )

    args = parser.parse_args()

    # Validate arguments
    if args.scenario == "historical":
        default_first_year, default_last_year = 1981, 2014
    else:  # ssp245
        default_first_year, default_last_year = 2015, 2100

    if args.first_year < default_first_year or args.last_year > default_last_year:
        print(f"Error: Years must be between {default_first_year} and {default_last_year} for {args.scenario} scenario",
              file=sys.stderr)
        sys.exit(1)

    scenario = args.scenario
    scenario_path = "Historical" if scenario == "historical" else "SSP245"

    # Handle token
    if args.token_file:
        try:
            with open(args.token_file, 'r') as f:
                token = f.read().strip()
        except FileNotFoundError:
            print(f"Error: Token file '{args.token_file}' not found",
                  file=sys.stderr)
            sys.exit(1)
    elif args.token:
        token = args.token
    else:
        # Throw error
        raise ValueError("Missing token file or token argument")

    # Create output directory
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Define area of interest
    area = {
        "lat": [args.area_lat_min, args.area_lat_max],
        "lon": [args.area_lon_min, args.area_lon_max]
    }

    print(f"Processing years {args.first_year} to {args.last_year}")
    print(f"Output directory: {out_dir}")
    print(f"Area of interest: lat {area['lat']}, lon {area['lon']}")

    # Get configurations
    models, calendars, micropoint_vars = get_model_configurations()
    dates_gregorian, dates_no_leap_days, dates_360days = create_date_arrays(default_first_year, default_last_year)

    # Process each year
    for year in range(args.first_year, args.last_year + 1):
        try:
            process_year(
                year, out_dir, token, models, calendars, micropoint_vars, dates_gregorian, dates_no_leap_days, area,
                scenario_path, scenario, default_first_year, default_last_year
            )
        except Exception as e:
            print(f"Error processing year {year}: {e}", file=sys.stderr)
            continue

    print("Processing complete!")


if __name__ == "__main__":
    main()