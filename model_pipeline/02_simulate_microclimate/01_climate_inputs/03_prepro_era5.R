#!/usr/bin/env Rscript

# ------ Download and process meteorological data - Command Line Version
# Usage:
# model_pipeline/02_simulate_microclimate/01_climate_inputs/03_prepro_era5.R path/to/raw/era5 path/to/processed/output year1 year2 ...

library(sf)
library(terra)
library(mcera5)
library(ncdf4)
library(tools)
library(ggplot2)
library(patchwork)
library(stringr)
library(microclimf)
library(microclimdata)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

dirin <- args[1]
dirout <- args[2]

# Validate output directory
if (!dir.exists(dirin)) {
  cat("Error: Input directory does not exist:", dirin, "\n")
  cat("Please create the directory first or provide a valid path.\n")
  quit(status = 1)
}

cat("Input directory:", dirin, "\n")

dir.create(dirout)

# Read all available remaining args and convert them to integers for year range
years <- c()
for (arg in args[3:length(args)]) {
  if (grepl("^\\d{4}$", arg)) {  # Check if the argument is a valid year
    years <- c(years, as.integer(arg))
  } else {
    cat("Warning: Ignoring invalid argument:", arg, "\n")
  }
}
cat("Years provided:", paste(years, collapse = ", "), "\n")

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

# Define output path
if (!dir.exists(dirin)) {
  dir.create(dirin, recursive = TRUE)
  cat("Created directory:", dirin, "\n")
}

# -- Process years

total_years <- length(years)

cat("Processing", total_years, "years of data...\n")

for (i in seq_along(years)) {
  year <- years[i]

  cat("\n--- Processing year", year, "(", i, "of", total_years, ") ---\n")

  # 1. List all files matching the year
  files <- list.files(
    path = dirin,
    pattern = paste0("(^|__)", year, "_[0-9]{1,2}\\.nc$"),
    full.names = TRUE
  )

  # 2. Extract month number
  months <- as.integer(str_extract(files, paste0("(?<=", year, "_)\\d{1,2}(?=\\.nc$)")))

  # 3. Order months
  sorted_indices <- order(months)
  files <- files[sorted_indices]
  months <- months[sorted_indices]

  # 4. Create new names (year + month only)
  if (length(months) == 0) {
    months <- seq(1, 12)
  }
  new_files <- file.path(
    dirin,
    sprintf("%d_%02d.nc", year, months)
  )

  # Rename (move) the files
  incorrect_names <- !file.exists(new_files)
  file.rename(files[incorrect_names], new_files[incorrect_names])

  # Check if months are duplicate
  unique_months <- !duplicated(months)
  new_files <- new_files[unique_months]
  months <- months[unique_months]

  if (length(new_files) != 12) {
    print(new_files)
    print(months)
    cat("Error: Expected 12 files for year", year, "but found months", months)
    next
  }

  tme <- as.POSIXlt(seq(
      as.POSIXct(paste0(year, "-01-01 00:00"), tz = "UTC"),
      as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
      by = "1 hour"
  ))

  # Recover req object for downloaded files to move on
  req <- build_era5_request(
    xmin = ext(r)$xmin,
    xmax = ext(r)$xmax,
    ymin = ext(r)$ymin,
    ymax = ext(r)$ymax,
    start_time = tme[1],
    end_time = tme[length(tme)],
    outfile_name = year
  )

  cat("Correcting NetCDF files...\n")

  # Correct req object and nc file
  for (month in seq(1, 12)) {
    month_str <- sprintf("%02d", month)
    req[[month]]$target <- paste0(year, "_", month_str, ".nc")

    # Correct subset name in nc files
    nc_path <- file.path(dirin, req[[month]]$target)
    if (file.exists(nc_path)) {
      tryCatch({
        rename_nc_subset(nc_path, "avg_sdlwrf", "msdwlwrf")
        cat("  Processed month", month, "\n")
      }, error = function(e) {
        cat("  Error processing month", month, ":", e$message, "\n")
      })

      # Check if dates match
      nc_file <- rast(nc_path)
      depth_values <- depth(nc_file)
      dates <- as.POSIXct(depth_values, origin = "1970-01-01", tz = "UTC")
      actual_month <- unique(format(dates, "%m"))
      cat("Actucal month:", actual_month, "/n")
      cat(length(dates), "\n")
      cat(length(dates) / 24, "\n")
      cat("----------")

      if (length(actual_month) == 1 && actual_month != month_str) {
        # construct corrected filename
        new_file <- file.path(dirin, paste0(year, "_", actual_month, ".nc"))

        # rename file
        file.rename(nc_path, new_file)

        # short message
        cat("⚠️ File renamed:", basename(nc_path), "→", basename(new_file), "\n\n")
      } else if (length(actual_month) > 1) {
        cat("⚠️ Multiple months detected in", basename(nc_path), ":", paste(actual_month, collapse = ", "), "\n\n")
      }

    } else {
      cat("  Warning: File not found for month", month, "\n")
    }
  }

  cat("Processing ERA5 data...\n")

  # Process data
  tryCatch({
    era5climdata <- era5_process(req, paste(dirin, "/", sep = ""), r, tme)

    # Save processed data
    output_file <- file.path(dirout, paste0("era5_climdata_", year, ".RDS"))
    saveRDS(era5climdata, output_file)
    cat("Saved processed data to:", output_file, "\n")

  }, error = function(e) {
    cat("Error processing data for year", year, ":", e$message, "\n")
  })

cat("\n=== Processing complete! ===\n")
cat("Output directory:", dirout, "\n")
cat("Processed", total_years, "years of ERA5 climate data.\n")