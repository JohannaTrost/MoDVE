#!/usr/bin/env Rscript

# ------ Download and process meteorological data - Command Line Version

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

# Read all available remaining args
# and convert them to integers for year range
# 2006, 2008, 2017, 2023
# 2092, 2078, 2042, 2089
# -> era5 here: 2024, 2022, 2021, 2020
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
library(stringr)
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

  # Time sequence
  tme <- as.POSIXlt(seq(
        as.POSIXct(paste0(year, "-01-01 00:00"), tz = "UTC"),
        as.POSIXct(paste0(year, "-12-31 23:00"), tz = "UTC"),
        by = "1 hour"
  ))

  # Define file prefix
  file_prefix <- paste0(gsub("-", "_", substr(Sys.time(), 1, 10)), "_")
  file <- list.files(
    path = dirin,
    pattern = paste0("__", year, "_12\\.nc$"),
    full.names = TRUE
  )
  file_prefix <- str_extract(file, "\\d{4}_\\d{2}_\\d{2}_")

  # Recover req object for downloaded files to move on
  req <- build_era5_request(
    xmin = ext(r)$xmin,
    xmax = ext(r)$xmax,
    ymin = ext(r)$ymin,
    ymax = ext(r)$ymax,
    start_time = tme[1],
    end_time = tme[length(tme)],
    outfile_name = file_prefix
  )

  cat("Correcting NetCDF files...\n")

  # Correct req object and nc file
  for (month in seq(1, 12)) {
    req[[month]]$target <- paste0(gsub("-", "_", file_prefix),
                                  paste0("_", year, "_"), month, ".nc")

    # Correct subset name in nc files
    nc_path <- file.path(dirin, req[[month]]$target)
    if (file.exists(nc_path)) {
      tryCatch({
        rename_nc_subset(nc_path, "avg_sdlwrf", "msdwlwrf")
        cat("  Processed month", month, "\n")
      }, error = function(e) {
        cat("  Error processing month", month, ":", e$message, "\n")
      })
    } else {
      cat("  Warning: File not found for month", month, "\n")
    }
  }

  cat("Processing ERA5 data...\n")

  # Process data
  tryCatch({
    era5climdata <- era5_process(req, paste(dirin, "/", sep = ""), r, tme)

    # Save processed data
    output_file <- file.path(dirout, paste0("/era5_climdata_", year, ".RDS"))
    saveRDS(era5climdata, output_file)
    cat("Saved processed data to:", output_file, "\n")

  }, error = function(e) {
    cat("Error processing data for year", year, ":", e$message, "\n")
  })

  # Optional: Create plots (commented out for command line efficiency)
  # Uncomment the following section if you want to generate plots

  # cat("Creating plots...\n")
  #
  # roi <- ext(-42.9, -42.2, -22.8, -22.2)
  #
  # # Example palettes
  # pal_temp <- colorRampPalette(c("white", "red"))(255)
  # pal_relhum <- colorRampPalette(c("white", "blue"))(255)
  # pal_pres <- terrain.colors(255)
  # pal_swdown <- heat.colors(255)
  # pal_difrad <- heat.colors(255)
  # pal_lwdown <- heat.colors(255)
  # pal_windspeed <- colorRampPalette(c("white", "darkgreen"))(255)
  # pal_winddir <- colorRampPalette(c("white", "purple"))(255)
  # pal_precip <- colorRampPalette(c("white", "blue", "darkblue"))(255)
  #
  # # Variable name to palette mapping
  # var_pal <- list(
  #   temp = pal_temp,
  #   relhum = pal_relhum,
  #   pres = pal_pres,
  #   swdown = pal_swdown,
  #   difrad = pal_difrad,
  #   lwdown = pal_lwdown,
  #   windspeed = pal_windspeed,
  #   winddir = pal_winddir,
  #   precip = pal_precip
  # )
  #
  # # Plotting function
  # plot_clim_var_roi <- function(r, title, pal, roi) {
  #   r_unwr <- unwrap(r)
  #   r_crop <- crop(r_unwr, roi)
  #   r_mean <- mean(r_crop, na.rm = TRUE)
  #   r_df <- as.data.frame(r_mean, xy = TRUE, na.rm = TRUE)
  #   colnames(r_df)[3] <- "value"
  #
  #   ggplot(r_df, aes(x = x, y = y, fill = value)) +
  #     geom_raster() +
  #     scale_fill_gradientn(colours = pal, name = title) +
  #     coord_equal() +
  #     theme_minimal() +
  #     labs(title = title) +
  #     theme(axis.text = element_blank(), axis.ticks = element_blank())
  # }
  #
  # # Loop over all variables
  # tryCatch({
  #   plots <- lapply(names(era5climdata), function(varname) {
  #     plot_clim_var_roi(era5climdata[[varname]], varname, var_pal[[varname]], roi)
  #   })
  #
  #   # Save combined plot
  #   combined_plot <- wrap_plots(plots, ncol = 3)
  #   plot_file <- file.path(climate_dir, paste0("era5_plots_", year, ".png"))
  #   ggsave(plot_file, combined_plot, width = 15, height = 10)
  #   cat("Saved plots to:", plot_file, "\n")
  # }, error = function(e) {
  #   cat("Error creating plots for year", year, ":", e$message, "\n")
  # })
}

cat("\n=== Processing complete! ===\n")
cat("Output directory:", dirout, "\n")
cat("Processed", total_years, "years of ERA5 climate data.\n")