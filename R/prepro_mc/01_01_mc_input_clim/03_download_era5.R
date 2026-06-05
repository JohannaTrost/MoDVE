#!/usr/bin/env Rscript

# ------ Download and process meteorological data - Command Line Version

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage: Rscript 03_download_era5.R <output_directory> <year>\n")
  cat("Example: Rscript 03_download_era5.R /path/to/output/directory 2042\n")
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

year <- args[2]  # Year to download, e.g., "2024"

# ERA5 data month and year that are to be downloaded
era5_months <- list(
  "2025" = 1:12
)

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
required_packages <- c("sf", "terra", "mcera5", "ncdf4", "tools", "ggplot2", "patchwork",
                       "foreach", "doParallel", "doRNG", "parallel")

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

# -- Process year

cat("Downloading", year, "data...\n")

# Time sequence
all_tme <- as.POSIXlt(seq(
      as.POSIXct(paste0(year, "-01-01 00:00"), tz = "UTC"),
      as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
      by = "1 hour"
))
if (year == 2001) {
    tme <- as.POSIXlt(seq(as.POSIXct("2001-01-01 00:00", tz="UTC"), as.POSIXct("2001-01-31 23:00", tz="UTC"), by="1 hour"))
} else {
  tme <- all_tme
    cat("Time range: Jan 1 - Dec 31,", year, "\n")
}

# Define file prefix
file_prefix <- paste0(gsub("-", "_", substr(Sys.time(), 1, 10)), "_")

cat("Downloading ERA5 data...\n")

# Download ERA5 data
req <- era5_download(r, tme, credentials, file_prefix, dirout)


# Check files month by month
year_str <- as.character(year)  # current year as string
months_to_check <- era5_months[[year_str]]

for (month in months_to_check) {
  month_str <- sprintf("%02d", month)
  nc_file_path <- file.path(dirout, paste0(file_prefix, year_str, "_", month, ".nc"))

  # Check if dates match
  nc_file <- rast(nc_file_path)
  depth_values <- depth(nc_file)
  dates <- as.POSIXct(depth_values, origin = "1970-01-01", tz = "UTC")
  actual_month <- unique(format(dates, "%m"))

  cat("Expected month:", month_str, "\n")
  cat("Actual month:", actual_month, "\n")
  cat(length(dates), "timestamps\n")
  cat(length(dates) / 24, "days\n")
  cat("----------\n")

  if (length(actual_month) == 1 && actual_month != month_str) {
    cat("⚠️ Mismatch:", month_str, "→", actual_month, "\n\n")
  } else if (length(actual_month) > 1) {
    cat("⚠️ Multiple months detected in", basename(nc_file_path), ":", paste(actual_month, collapse = ", "), "\n\n")
  }
}

cat("\n=== Downloading complete! ===\n")
cat("Output directory:", dirout, "\n")
cat("Downloaded", year, " ERA5 climate data.\n")