#!/usr/bin/env Rscript

# ------ Download and process meteorological data - Command Line Version

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage: Rscript era5_download.R <output_directory>\n")
  cat("Example: Rscript era5_download.R /path/to/output/directory\n")
  quit(status = 1)
}

dirout <- args[1]

# Validate output directory
if (!dir.exists(dirout)) {
  cat("Error: Output directory does not exist:", dirout, "\n")
  cat("Please create the directory first or provide a valid path.\n")
  quit(status = 1)
}

cat("Output directory:", dirout, "\n")

first_year <- as.integer(args[2])
last_year <- as.integer(args[3])

if (first_year < 1940 || last_year > 2024 || first_year > last_year) {
  cat("Error: Invalid year range. Please provide years between 1940 and 2024.\n")
  quit(status = 1)
}

# obs_time – UTC POSIXlt object of observation times for each climate variable, 2017-01-01 00:00:00
# temp – temperatures (deg C)
# relhum - relative humidity (percentage)
# pres - atmospheric pressure (kPa)
# swdown - total downward shortwave radiation received by a horizontal surface (W/m^2)
# difrad - diffuse radiation (W/m^2)
# lwdown - total downward longward radiation (W/m^2)
# windspeed - wind speed at reference height (m/s)
# winddir - wind direction in degrees
# precip - hourly precipitation (mm).

# -- Install and load libraries

cat("Loading required libraries...\n")

# Check and install required packages
required_packages <- c("sf", "terra", "mcera5", "ncdf4", "tools", "ggplot2", "patchwork")
# needs install from github "microclimf", "microclimdata"

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

library(microclimf)
library(microclimdata)

# -- Helper functions

rename_nc_subset <- function(file_path, old_name, new_name) {

  tmp_path <- tempfile(fileext = ".nc")

  # Open the original NetCDF file
  nc_old <- nc_open(file_path)

  # Gather dimensions
  dims <- nc_old$dim
  dim_defs <- lapply(dims, function(d) ncdim_def(d$name, d$units, d$vals))

  # Gather variables
  vars <- nc_old$var
  var_defs <- list()
  for (v in names(vars)) {
    var <- vars[[v]]
    var_name <- if (v == old_name) new_name else v
    var_defs[[var_name]] <- ncvar_def(
      name = var_name,
      units = var$units,
      dim = var$dim,
      missval = var$missval,
      longname = var$longname,
      prec = var$prec
    )
  }

  # Create new NetCDF file with updated variable name
  nc_new <- nc_create(tmp_path, vars = var_defs)

  # Copy over data for each variable
  for (v in names(vars)) {
    old_name_v <- v
    new_name_v <- if (v == old_name) new_name else v
    data <- ncvar_get(nc_old, old_name_v)
    ncvar_put(nc_new, new_name_v, data)
  }

  # Close files
  nc_close(nc_old)
  nc_close(nc_new)

  # Replace the original file with the modified one
  file.copy(tmp_path, file_path, overwrite = TRUE)
  unlink(tmp_path)  # Clean up
}

# -- Setup credentials

cat("Setting up credentials...\n")

# Access datasets for download
data(credentials)

uid <- "05d8ef03-8aa5-4508-af48-4c740eeff0d3"
access_token <- "ac0fedd8-e162-4742-949a-a3f40e26a68b"

credentials[3, "username"] <- uid
credentials[3, "password"] <- access_token

# -- Climate Data Setup

cat("Setting up spatial extent and projection...\n")

# Define extent in decimal degrees
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
r <- rast(e)
crs(r) <- "EPSG:4326"

# Reproject to meters
r <- project(r, "EPSG:31983", res = 10)

# -- Process years

years <- seq(first_year, last_year + 1)
total_years <- length(years)

cat("Downloading", total_years, "years of data...\n")

for (i in seq_along(years)) {
  year <- years[i]

  cat("\n--- Downloading year", year, "(", i, "of", total_years, ") ---\n")

  # Time sequence
  all_tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-01-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
  ))
  if (year == 1981) {
      tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-11-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
      ))
      cat("Time range: Nov 1 - Dec 31,", year, "\n")
  } else if (year == 2001) {
      tme <- as.POSIXlt(seq(
          as.POSIXct(paste0(year, "-10-01 00:00"), tz = "UTC"),
          as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
          by = "1 hour"
      ))
  } else {
    tme <- all_tme
      cat("Time range: Jan 1 - Dec 31,", year, "\n")
  }

  # Define file prefix
  file_prefix <- paste0(gsub("-", "_", substr(Sys.time(), 1, 10)), "_")

  cat("Downloading ERA5 data...\n")

  # Download ERA5 data
  req <- era5_download(r, tme, credentials, file_prefix, dirout)
}

cat("\n=== Downloading complete! ===\n")
cat("Output directory:", dirout, "\n")
cat("Downloaded", total_years, "years of ERA5 climate data.\n")