#require(devtools)
# install_github("ilyamaclean/microclimf")
# remotes::install_local("/Users/johanna/Uni/masterarbeit/code/forks/micropoint",
#                        force = TRUE)
#!/usr/bin/env Rscript

# Microclimate Simulation CLI with TOML Configuration
# Usage: Rscript microclimate_cli.R --config config.toml [--verbose]

suppressPackageStartupMessages({
  library(optparse)
  library(RcppTOML)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(terra)
  library(micropoint)
  library(readr)
  library(viridis)
  library(microclimf)
  library(dplyr)
  library(lubridate)
})

get_xy_pairs <- function(xdim, ydim, outdir) {
  missing_x <- c()
  missing_y <- c()

  # 1. Check if folder is empty
  if (length(list.files(outdir)) == 0) {
      message("Folder is empty. Marking all pairs as missing.")
      for (x in seq(xdim)) {
          for (y in seq(ydim)) {
              missing_x <- c(missing_x, x)
              missing_y <- c(missing_y, y)
          }
      }
  } else {
      # 2. Loop through and check files
      for (x in seq(xdim)) {
          for (y in seq(ydim)) {
              out_path <- file.path(outdir, paste0("mc_x", x, "_y", y, ".rds"))

              if (!file.exists(out_path)) {
                  # Missing file
                  missing_x <- c(missing_x, x)
                  missing_y <- c(missing_y, y)
              } else {
                  # 3. Check contents
                mc <- readRDS(out_path)
                if (is.null(mc$data)) {
                    missing_x <- c(missing_x, x)
                    missing_y <- c(missing_y, y)
                }
            }
        }
    }

  }
  message("No. missing cells: ", length(missing_x))

  return(list(x = missing_x, y = missing_y))
}

# Vectorized function to calculate mean diurnal temperature range
mean_diurnal_temp_vectorized <- function(tair_matrix, timestamps) {
  # Convert to dates for grouping
  dates <- as.Date(timestamps)
  unique_dates <- unique(dates)

  # Pre-allocate results
  n_dates <- length(unique_dates)
  n_heights <- nrow(tair_matrix)
  daily_ranges <- matrix(NA, nrow = n_heights, ncol = n_dates)

  # Calculate daily min/max for all heights simultaneously
  for (i in seq_along(unique_dates)) {
    day_indices <- which(dates == unique_dates[i])
    if (length(day_indices) >= 12) {  # Filter days with sufficient data
      day_data <- tair_matrix[, day_indices, drop = FALSE]
      daily_max <- apply(day_data, 1, max, na.rm = TRUE)
      daily_min <- apply(day_data, 1, min, na.rm = TRUE)
      daily_ranges[, i] <- daily_max - daily_min
    }
  }

  # Return mean diurnal range for each height
  apply(daily_ranges, 1, mean, na.rm = TRUE)
}

# Vectorized monthly statistics function
month_stats_vectorized <- function(data_matrix, timestamps) {
  # Get months from timestamps
  months <- month(timestamps)
  unique_months <- 1:12

  n_heights <- nrow(data_matrix)
  n_months <- length(unique_months)

  # Pre-allocate monthly means matrix
  monthly_means <- matrix(NA, nrow = n_heights, ncol = n_months)
  monthly_max <- matrix(NA, nrow = n_heights, ncol = n_months)
  monthly_min <- matrix(NA, nrow = n_heights, ncol = n_months)

  # Calculate monthly statistics for all heights simultaneously
  for (m in unique_months) {
    month_indices <- which(months == m)
    if (length(month_indices) > 0) {
      month_data <- data_matrix[, month_indices, drop = FALSE]
      monthly_means[, m] <- apply(month_data, 1, mean, na.rm = TRUE)
      monthly_max[, m] <- apply(month_data, 1, max, na.rm = TRUE)
      monthly_min[, m] <- apply(month_data, 1, min, na.rm = TRUE)
    }
  }

  # Calculate derived metrics for all heights
  seasonality <- apply(monthly_means, 1, sd, na.rm = TRUE) * 100

  # Find max and min values across months for each height
  max_of_warmest <- apply(monthly_max, 1, max, na.rm = TRUE)
  min_of_coldest <- apply(monthly_min, 1, min, na.rm = TRUE)
  annual_range <- max_of_warmest - min_of_coldest

  return(list(
    seasonality = seasonality,
    max_value = max_of_warmest,
    min_value = min_of_coldest,
    annual_range = annual_range
  ))
}

# Function to convert hourly to daily medians more efficiently
hourly_to_daily_medians <- function(data_matrix, timestamps) {
  dates <- as.Date(timestamps)
  unique_dates <- unique(dates)

  n_heights <- nrow(data_matrix)
  n_days <- length(unique_dates)
  daily_medians <- matrix(NA, nrow = n_heights, ncol = n_days)

  for (i in seq_along(unique_dates)) {
    day_indices <- which(dates == unique_dates[i])
    if (length(day_indices) > 0) {
      day_data <- data_matrix[, day_indices, drop = FALSE]
      daily_medians[, i] <- apply(day_data, 1, median, na.rm = TRUE)
    }
  }

  return(list(data = daily_medians, timestamps = unique_dates))
}


aggregate_mc <- function(mc, timestamps, n_temp_metrics) {

  start_time_cell <- Sys.time()

  actual_max_hgt <- dim(mc$tair)[1]

  # Pre-allocate result matrix for this cell
  cell_matrix <- array(NA, dim = c(actual_max_hgt, n_temp_metrics))

  # Verify data dimensions
  if (dim(mc$tair)[2] != length(timestamps)) {
    stop(paste("Timestamp mismatch in file:", mc_file))
  }

  # Process all heights simultaneously where possible

  # 1. Mean annual temperature (vectorized)
  cell_matrix[, 1] <- apply(mc$tair, 1, mean, na.rm = TRUE)

  # 2. Mean diurnal temperature range (vectorized)
  cell_matrix[, 2] <- mean_diurnal_temp_vectorized(mc$tair, timestamps)

  # Convert to daily medians for remaining calculations
  daily_tair <- hourly_to_daily_medians(mc$tair, timestamps)
  daily_relhum <- hourly_to_daily_medians(mc$relhum, timestamps)
  daily_windspeed <- hourly_to_daily_medians(mc$windspeed, timestamps)

  # 3-6. Temperature statistics (vectorized)
  tair_stats <- month_stats_vectorized(daily_tair$data, daily_tair$timestamps)
  cell_matrix[, 3] <- tair_stats$annual_range
  cell_matrix[, 4] <- tair_stats$max_value
  cell_matrix[, 5] <- tair_stats$min_value
  cell_matrix[, 6] <- (cell_matrix[, 2] / tair_stats$annual_range) * 100  # Isothermality

  # 7-10. Humidity statistics (vectorized)
  relhum_stats <- month_stats_vectorized(daily_relhum$data, daily_relhum$timestamps)
  cell_matrix[, 7] <- apply(daily_relhum$data, 1, mean, na.rm = TRUE)
  cell_matrix[, 8] <- relhum_stats$annual_range
  cell_matrix[, 9] <- relhum_stats$max_value
  cell_matrix[, 10] <- relhum_stats$min_value

  # 11-14. Wind speed statistics (vectorized)
  ws_stats <- month_stats_vectorized(daily_windspeed$data, daily_windspeed$timestamps)
  cell_matrix[, 11] <- apply(daily_windspeed$data, 1, mean, na.rm = TRUE)
  cell_matrix[, 12] <- ws_stats$annual_range
  cell_matrix[, 13] <- ws_stats$max_value
  cell_matrix[, 14] <- ws_stats$min_value

  end_time_cell <- Sys.time()
  processing_time <- end_time_cell - start_time_cell
  #cat("Finished processing cell:", x, y, "in", round(processing_time, 2), "seconds\n")

  return(list(
    data = cell_matrix,
    processing_time = processing_time
  ))
}

immediateMessage <- function(..., domain = NULL, appendLF = TRUE) {
  msg <- .makeMessage(..., domain = domain, appendLF = appendLF)
  call <- sys.call()
  m <- simpleMessage(msg, call)

  cls <- class(m)
  cls <- setdiff(cls, "condition")
  cls <- c(cls, "immediateCondition", "condition")
  class(m) <- cls

  message(m)
  invisible(m)
}


get_monthly_mc_stats <- function(min_arr, max_arr, month_labs_mn, 
                                 month_labs_mx){
  months <- unique(month_labs_mn)
  # Initialize 4D output array: 50 x 50 x 12 x 3
  result_arr <- array(NA, 
                      dim = c(dim(min_arr)[1:2], length(months), 3))
  
  for (i in seq_along(months)) {
    m <- months[i]
    
    # Get indices for this month
    idx_mn <- which(month_labs_mn == m)
    idx_mx <- which(month_labs_mx == m)
    
    # Compute means for this month
    min_mean <- apply(min_arr[,,idx_mn, drop = FALSE], c(1,2), mean, na.rm = TRUE)
    max_mean <- apply(max_arr[,,idx_mx, drop = FALSE], c(1,2), mean, na.rm = TRUE)
    avg_mean <- (min_mean + max_mean) / 2
    
    result_arr[,,i,1] <- min_mean
    result_arr[,,i,2] <- max_mean
    result_arr[,,i,3] <- avg_mean
  }
  
  return(result_arr)
}


extract_params <- function(raster_list, lon, lat, crs = "EPSG:4326") {
  # Create a SpatVector point in WGS84 (decimal degrees)
  point <- terra::vect(cbind(lon, lat), crs = crs)
  
  # Transform the point to the CRS of the rasters 
  target_crs <- crs(unwrap(raster_list[[1]]))
  point_proj <- project(point, target_crs)
  
  # Extract the value at the projected point for each raster
  result <- lapply(raster_list, function(r) {
    r <- unwrap(r)
    terra::extract(r, point_proj)[[2]]
  })

  return(result)
}


indices2coords <- function(x, y, raster, crs_out = "EPSG:4326") {
  # Get the cell number from row and column indices
  cell <- terra::cellFromRowCol(raster, row = x, col = y)
  
  # Get the center coordinates of the cell in the raster's CRS
  coords <- terra::xyFromCell(raster, cell)
  
  # Create point in specified crs
  pts <- terra::vect(coords, type = "points", crs = terra::crs(raster))
  pts_proj <- terra::project(pts, crs_out)
  
  # Extract coordinates
  coords_proj <- terra::geom(pts_proj)[, c("x", "y")]
  names(coords_proj) <- c("x", "y")
  
  return(coords_proj)
}

# ----------------------------------------------------------
# Set up parallel processing
# ----------------------------------------------------------

# Function to process a single (x, y) cell
process_cell <- function(x, y, year, n_temp_metrics = 14, microhab_path,
                         vegp_path, soilc_path, climdata_path, out_dir) {

  start_time <- Sys.time()

  out_path <- paste0(out_dir, "/mc_x", x, "_y", y, ".rds")
  if (file.exists(out_path)) {
    agg_res <- readRDS(out_path)
    if (!is.null(agg_res$data)) {
      immediateMessage(paste("Output file already exists for cell:", x, y))
      return(list(
        x = x,
        y = y,
        data = agg_res$data,
        processing_time = agg_res$processing_time
      ))
    }
  }

  library(microclimf) # Have to load inside the function otherwise there is "future" error

  # Load data for one year
  vegp_reg <- readRDS(vegp_path)
  soilc_reg <- readRDS(soilc_path)
  climdata_reg <- read_csv(climdata_path) %>% filter(year(obs_time) == year)
  pai <- readRDS(microhab_path)[,,,4]

  # Create hourly timestamps for 2024
  timestamps <- seq(ymd_h(paste0(year, "-01-01 00")), ymd_h(paste0(year, "-12-31 23")), by = "hour")

  # Veg heights
  max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
  heights <- seq(0.5, max_veg_height + 1)

  immediateMessage(paste("Processing cell:", x, y))

  coords <- indices2coords(x, y, terra::unwrap(vegp_reg$pai))[c("x", "y")]
  lon <- coords[[1]]
  lat <- coords[[2]]

  vegparams <- extract_params(vegp_reg, lon, lat)
  grndparams <- extract_params(soilc_reg, lon, lat)
  #paii <- pai[x, y, 1:max(vegparams$h, 0.5)]
  max_hgt <- max(values(terra::unwrap(vegp_reg$h)))
  paii <- pai[x, y, 1:max_hgt]

    res <- list()

  for (i in seq_along(heights)) {
    h <- heights[i]
    immediateMessage(paste(x, y, "height:", h))

    mout <- micropoint::runpointmodel(climdata_reg, reqhgt = h, vegparams,
                                      paii, grndparams, lat = lat, long = lon)

    immediateMessage("Point model run completed.")

    # Which variables to save
    if (i == 1 & x == 1 & y == 1) {
      exp_vars <- c("tair", "relhum", "windspeed", "obs_time")
    } else {
      exp_vars <- c("tair", "relhum", "windspeed")
    }

    for (var in exp_vars) {
      message(paste(x, y, var, h))
      if (i == 1) {
        res[[var]] <- array(NA, dim = c(length(heights), length(mout[[var]])))
      }
      res[[var]][i, ] <- mout[[var]]
    }
  }

  immediateMessage("Point model run completed.")

  # aggregate yearly microclimate data
  agg_res <- aggregate_mc(res, timestamps, n_temp_metrics) # Expects heights x time and not time x heights

  end_time <- Sys.time()
  processing_time <- end_time - start_time
  immediateMessage(paste("Finished cell: x:", x, ", y:", y, "in", round(processing_time, 2), "seconds"))

  # Save results to file
  saveRDS(agg_res, file = out_path)

  return(list(
    x = x,
    y = y,
    data = agg_res$data,
    processing_time = processing_time
  ))
}

# Define command line options
option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.toml",
              help = "Path to TOML configuration file [default %default]", metavar = "FILE"),
  make_option(c("-k", "--chunk"), type = "integer", default = "10",
              help = "Number of chunk that determines which cells to simulate [default %default]"),
  make_option(c("-y", "--year"), type = "integer", default = "2024",
              help = "Year for which microclimate will be simulated (must be within the climdata)  [default %default]"),
  make_option(c("-t", "--timestep"), type = "integer", default = "1",
              help = "Time step to take from forest simulation [default %default]"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Enable verbose output")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                          description = "Microclimate simulation processing script using TOML configuration")

opt <- parse_args(opt_parser)

# Check if config file exists
if (!file.exists(opt$config)) {
  cat("Error: Configuration file not found:", opt$config, "\n")
  quit(status = 1)
}

# Load configuration
tryCatch({
  config <- parseTOML(opt$config)
}, error = function(e) {
  cat("Error parsing TOML configuration file:", opt$config, "\n")
  cat("Error message:", e$message, "\n")
  quit(status = 1)
})

# Validate required sections
required_sections <- c("simulation", "processing", "paths")
missing_sections <- setdiff(required_sections, names(config))
if (length(missing_sections) > 0) {
  cat("Error: Missing required sections in configuration file:\n")
  for (section in missing_sections) {
    cat("  -", section, "\n")
  }
  quit(status = 1)
}

# Extract configuration values with defaults
# Simulation parameters
x_dim <- ifelse(is.null(config$simulation$x_dim), 50, config$simulation$x_dim)
y_dim <- ifelse(is.null(config$simulation$y_dim), 50, config$simulation$y_dim)
n_temp_metrics <- ifelse(is.null(config$simulation$n_temp_metrics), 14, config$simulation$n_temp_metrics)
year <- opt$year
ts <- opt$timestep

# Processing parameters
chunk <- ifelse(is.null(opt$chunk), 1, opt$chunk)
chunk_size <- ifelse(is.null(config$processing$chunk_size), 2, config$processing$chunk_size)
cores_config <- config$processing$cores

# Path parameters
region <- ifelse(is.null(config$paths$region), "pirineus", config$paths$region)
veg_indir_base <- ifelse(is.null(config$paths$veg_indir),
                        "/Users/johanna/Uni/masterarbeit/data/modve_output",
                        config$paths$veg_indir)
in_dir_base <- ifelse(is.null(config$paths$in_dir),
                     "/Users/johanna/Uni/masterarbeit/data/mc_input",
                     config$paths$in_dir)
climdata_path <- ifelse(is.null(config$paths$clim_path),
                    "/Users/johanna/Uni/masterarbeit/data/mc_input/pirineus/scenarios/climdata_era5_cmip6_1906-2024_ssp245_119ts_v1.csv",
                    config$paths$clim_path)
outdir_base <- ifelse(is.null(config$paths$outdir),
                     "/Users/johanna/Uni/masterarbeit/data/mc_output/v5",
                     config$paths$outdir)

# Output parameters (can be overridden by command line)
verbose <- opt$verbose || (!is.null(config$output$verbose) && config$output$verbose)

# ---------------------------------------------------------- Configure

# Define full directory paths
veg_indir <- file.path(veg_indir_base, region, "a1_1")
in_dir <- file.path(in_dir_base, region)
outdir <- file.path(outdir_base, region, year)

if (verbose) {
  cat("Configuration loaded from:", opt$config, "\n")
  cat("Settings:\n")
  cat("  Region:", region, "\n")
  cat("  Year:", year, "\n")
  cat("  Dimensions:", x_dim, "x", y_dim, "\n")
  cat("  Time step:", ts, "\n")
  cat("  Temperature metrics:", n_temp_metrics, "\n")
  cat("  Chunk:", chunk, "of size", chunk_size, "\n")
  cat("  Vegetation input dir:", veg_indir, "\n")
  cat("  Main input dir:", in_dir, "\n")
  cat("  Output directory:", outdir, "\n")
}

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  if (verbose) cat("Created output directory:", outdir, "\n")
}

# Define file paths
microhab_path <- file.path(veg_indir, paste0("MicrohabitatMatrix", ts, ".rds"))
vegp_path <- file.path(in_dir, paste0("vegp_mof3d_ptm_", ts, "_v4.RDS"))
soilc_path <- file.path(in_dir, "soilc_v2.RDS")

# Validate input files exist
required_files <- c(microhab_path, vegp_path, soilc_path, climdata_path)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  cat("Error: Missing required input files:\n")
  for (file in missing_files) {
    cat("  -", file, "\n")
  }
  quit(status = 1)
}

if (verbose) {
  cat("All required input files found.\n")
}

# Save heights of MC simulations
vegp_reg <- readRDS(vegp_path)
max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)
max_hgt <- length(heights)

heights_file <- file.path(outdir, "mc_heights.rds")
if (!file.exists(heights_file)) {
  saveRDS(heights, heights_file)
  if (verbose) cat("Saved heights to:", heights_file, "\n")
}

# Get all cells for which we need to run the microclimate model
cells <- get_xy_pairs(x_dim, y_dim, outdir)

# Get the ith chunk of cells
first_cells_idx <- (chunk - 1) * chunk_size + 1
last_cells_idx <- min(chunk * chunk_size, length(cells$x))

cat("Simulating microclimate for cells", first_cells_idx, "to", last_cells_idx, "\n")

# Determine number of cores
if (cores_config == "auto") {
  nTasks <- Sys.getenv("SLURM_NTASKS")
  if (nTasks != "") {
    numCores <- strtoi(nTasks)
  } else {
    numCores <- parallel::detectCores() - 1
  }
} else {
  numCores <- cores_config
}

registerDoParallel(numCores)
cat("Using", numCores, "cores for parallel processing.\n")

# Main processing loop
output <- foreach(pair_idx = seq(first_cells_idx, last_cells_idx),
                 .export = c("process_cell", "mean_diurnal_temp_vectorized", "month_stats_vectorized",
                            "hourly_to_daily_medians", "aggregate_mc", "immediateMessage",
                            "get_monthly_mc_stats", "extract_params", "indices2coords")) %dorng% {
  x <- cells$x[pair_idx]
  y <- cells$y[pair_idx]

  c <- process_cell(x, y, year, n_temp_metrics, microhab_path,
                   vegp_path, soilc_path, climdata_path, outdir)

  out_path <- paste0(outdir, "/mc_x", x, "_y", y, ".rds")
  mc <- readRDS(out_path)

  if (is.null(mc$data)) {
    warning(paste("No data for cell:", x, y))
  } else if (verbose) {
    message(paste("Processed cell:", x, y, "with data."))
  }

  return(NULL)
}

# Check for missing cells after processing
missing_x <- c()
missing_y <- c()
total_nas <- 0

for (i in seq(first_cells_idx, last_cells_idx)) {
  x <- cells$x[i]
  y <- cells$y[i]
  out_path <- paste0(outdir, "/mc_x", x, "_y", y, ".rds")

  if (file.exists(out_path)) {
    mc <- readRDS(out_path)
    if (is.null(mc$data)) {
      missing_x <- c(missing_x, x)
      missing_y <- c(missing_y, y)
    } else {
      # Count NAs in data
      nas_count <- sum(is.na(mc$data))
      total_nas <- total_nas + nas_count
      if (nas_count > 0 && verbose) {
        message(paste("Cell:", x, y, "has", nas_count, "missing values in data."))
      }
    }
  } else {
    warning(paste("Output file not found for cell:", x, y))
    missing_x <- c(missing_x, x)
    missing_y <- c(missing_y, y)
  }
}

# Summary report
cat("\n--- Processing Summary ---\n")
cat("Configuration file:", opt$config, "\n")
cat("Processed cells:", first_cells_idx, "to", last_cells_idx, "\n")
cat("Missing cells:", length(missing_x), "\n")
cat("Total NA values across all cells:", total_nas, "\n")

if (length(missing_x) > 0) {
  cat("Missing cell coordinates:\n")
  for (i in seq_along(missing_x)) {
    cat("  x =", missing_x[i], ", y =", missing_y[i], "\n")
  }
  quit(status = 1)
} else {
  cat("All cells processed successfully!\n")
  quit(status = 0)
}