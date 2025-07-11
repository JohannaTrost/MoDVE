#require(devtools)
# install_github("ilyamaclean/microclimf")
remotes::install_local("/Users/johanna/Uni/masterarbeit/code/forks/micropoint",
                       force = TRUE)

library(micropoint)
library(terra)
library(readr)
library(viridis)
library(microclimf)
library(dplyr)
library(lubridate)
library("foreach")
library("doParallel")
library("doRNG")


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
    extract(r, point_proj)[[2]]
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

  out_path <- paste0(out_dir, "/v3_mc_x", x, "_y", y, ".rds")
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
  climdata_reg <- read_csv(climdata_path)
  pai <- readRDS(microhab_path)[,,,5]

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
    print(h)

    mout <- micropoint::runpointmodel(climdata_reg, reqhgt = h, vegparams,
                                      paii, grndparams, lat = lat, long = lon)

    immediateMessage("Point model run completed.")

    # Which variables to save
    if (i == 1 & x == 1 & y == 1) {
      exp_vars <- c("tair", "tcanopy", "relhum", "windspeed", "obs_time")
    } else {
      exp_vars <- c("tair", "tcanopy", "relhum", "windspeed")
    }

    for (var in exp_vars) {
      message(paste(x, y, var, h))
      if (i == 1) {
        res[[var]] <- array(NA, dim = c(length(heights), length(mout[[var]])))
      }
      res[[var]][i, ] <- mout[[var]]
    }
  }

  # aggregate yearly microclimate data
  agg_res <- aggregate_mc(res, timestamps, n_temp_metrics)

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

# vegp 1m resolution
# dtm 1m resolution
# soilc 1m resolution
# climdata is equal for the whole grid

# ---------------------------------------------------------- Configure

x_dim <- 50
y_dim <- 50
n_temp_metrics <- 14
year <- 2024

# Define output directory
outdir <- "/Users/johanna/Uni/masterarbeit/data/mc_output"

# Save heights of MC simulations
in_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
microhab_path <- "/Users/johanna/Uni/masterarbeit/output/modev_zach_25_01_07/MicrohabitatMatrix99.rds"
vegp_path <- paste(in_dir, "vegp_mof3d_ptm_v3.RDS", sep = "/")
soilc_path <- paste(in_dir, "soilc_v2.RDS", sep = "/")
climdata_path <- paste(in_dir, "era5_climdata_2024_v2.csv", sep = "/")

# Get forest heights for which microclimate will be simulated
vegp_reg <- readRDS(vegp_path)
max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)
max_hgt <- length(heights)
saveRDS(heights, paste0(outdir, "/v2_mc_heights.rds"))

missing_x <- c()
missing_y <- c()
for (x in seq(50)) {
    for (y in seq(50)) {
        out_path <- paste0(outdir, "/v3_mc_x", x, "_y", y, ".rds")
        mc <- readRDS(out_path)
        if(is.null(mc$data)) {
            missing_x <- c(missing_x, x)
            missing_y <- c(missing_y, y)
        }
    }
}
message("No. missing cells: ", length(missing_x))

nTasks <- Sys.getenv("SLURM_NTASKS")
if (nTasks != "") {
    numCores <- strtoi(nTasks)
} else {
    numCores <- detectCores() - 1
}

registerDoParallel(numCores)

output <- foreach (pair_idx=seq_len(length(missing_x)),
                   .export=c("process_cell", "mean_diurnal_temp_vectorized", "month_stats_vectorized",
                             "hourly_to_daily_medians", "aggregate_mc", "immediateMessage",
                             "get_monthly_mc_stats", "extract_params", "indices2coords"
                   )) %dorng% {
    x <- missing_x[pair_idx]
    y <- missing_y[pair_idx]

    c <- process_cell(x, y, year, n_temp_metrics, microhab_path,
                 vegp_path, soilc_path, climdata_path, outdir)

    out_path <- paste0(outdir, "/v3_mc_x", x, "_y", y, ".rds")
    mc <- readRDS(out_path)
    if(is.null(mc$data)) {
      warning(paste("No data for cell:", x, y))
    } else {
      message(paste("Processed cell:", x, y, "with data."))
    }
    return(NULL)
}

# Check for missing cells after processing
missing_x <- c()
missing_y <- c()
for (x in seq(50)) {
    for (y in seq(50)) {
        out_path <- paste0(outdir, "/v3_mc_x", x, "_y", y, ".rds")
        mc <- readRDS(out_path)
        if(is.null(mc$data)) {
            missing_x <- c(missing_x, x)
            missing_y <- c(missing_y, y)
        }
        else {
            # Count nas in data
            nas_count <- sum(is.na(mc$data))
            if (nas_count > 0) {
                message(paste("Cell:", x, y, "has", nas_count, "missing values in data."))
            }
        }
    }
}
message("No. missing cells: ", length(missing_x))


# DEEBUG new method

# Load data for one year
vegp_reg <- readRDS(vegp_path)
soilc_reg <- readRDS(soilc_path)
climdata_reg <- read_csv(climdata_path)
pai <- readRDS(microhab_path)[,,,5]

# Check for missing cells after processing
nImgs <- 0
while (nImgs < 4) {
  x <- sample(1:50, 1)
  y <- sample(1:50, 1)
  out_path <- paste0(outdir, "/v3_mc_x", x, "_y", y, ".rds")
  mc <- readRDS(out_path)
  if(!is.null(mc$data)) {

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

      # try plotprofile
      ppout <- plotprofile(climdata_reg, hr = 1029, "relhum",
                            vegparams, paii = paii[1:vegparams$h], grndparams, lat = lat, long= lon)
      #
      #pmout <- runprofilemodel(climdata_reg, vegparams,
      #                           paii = paii[1:vegparams$h], grndparams, lat = lat, long= lon)

      pmout <- micropoint::runpointmodel(climdata_reg, reqhgt = seq(0.5, length(paii)), vegparams,
                                         paii, grndparams, lat = lat, long = lon)

      # Extract x values from all three line data sources
      x1 <- ppout$var
      x2 <- mc$data[1:length(ppout$z), 7] # 1 airt , 7 relhum, 11 windspeed
      #x3 <- apply(pmout$profile[,,2], 2, mean)
      x3 <- apply(pmout$relhum, 2, mean)

      # Determine the global x-axis range
      x_range <- range(c(x2, x3), na.rm = TRUE)
      y_range <- range(c(ppout$z, seq(length(pmout$height))), na.rm = TRUE)

      pdf(paste0("../../figs/mc_output/test_profile_relhum_cppfct4_comp_", x, "_", y, ".pdf"))

      # Plot the first line with axis labels
      plot(ppout$z ~ x2, type = "l", xlim = x_range, ylim = y_range,
           xlab = "Relative Humidity (%)", ylab = "Height (m)")

      # Add the second line
      lines(pmout$heights ~ x3, col = "red")

      # Add a legend
      legend("topright",
             legend = c("Original 'runmodel' for each height", "New cpp 'runmodelProfile'"),
             col = c("black", "red"), lty = 1, bty = "n")

      dev.off()

      nImgs <- nImgs + 1
  }
}

# Benchmarking the speed of the new function

library(microbenchmark)

# Simulate one iteration's data:
x <- sample(1:50, 1)
y <- sample(1:50, 1)

coords <- indices2coords(x, y, terra::unwrap(vegp_reg$pai))[c("x", "y")]
lon <- coords[[1]]
lat <- coords[[2]]

vegparams <- extract_params(vegp_reg, lon, lat)
grndparams <- extract_params(soilc_reg, lon, lat)
max_hgt <- max(values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
paii <- pai[x, y, 1:max_hgt]
heights <- seq(0.5, length(paii))

# Benchmark both approaches
benchmark_results <- microbenchmark(
  new_approach = {
    pmout <- micropoint::runpointmodel(
      climdata_reg, reqhgt = heights, vegparams,
      paii, grndparams, lat = lat, long = lon
    )
  },
  old_approach = {
    for (h in heights) {
      mout <- micropoint::runpointmodel(
        climdata_reg, reqhgt = h, vegparams,
        paii, grndparams, lat = lat, long = lon
      )
    }
  },
  times = 10L # Repeat 10 times for statistical accuracy
)

print(benchmark_results)

library(ggplot2)
pdf("../../figs/benchmarking/mc_benchmarking.pdf")
autoplot(benchmark_results)
dev.off()
