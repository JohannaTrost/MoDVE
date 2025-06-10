library(dplyr)
library(tidyr)
library(lubridate)
library(parallel)

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

# Optimized processing function for a single cell
process_single_cell <- function(x, y, timestamps, x_dim, y_dim, max_hgt, n_temp_metrics) {

  start_time_cell <- Sys.time()

  # Load microclimate data
  mc_file <- paste0("/Users/johanna/Uni/masterarbeit/data/mc_output/v2_mc_x", x, "_y", y, ".rds")

  if (!file.exists(mc_file)) {
    warning(paste("File not found:", mc_file))
    return(NULL)
  }

  mc <- readRDS(mc_file)
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
    x = x,
    y = y,
    data = cell_matrix,
    processing_time = processing_time
  ))
}

# Main processing function with optional parallelization
process_microclimate_data <- function(x_dim = 50, y_dim = 50, use_parallel = TRUE, n_cores = NULL) {

  # Create hourly timestamps for 2024
  timestamps <- seq(ymd_h("2024-01-01 00"), ymd_h("2024-12-31 23"), by = "hour")

  # Initialize the microclimate matrix
  max_hgt <- 80
  n_temp_metrics <- 14
  mc_matrix <- array(rep(NA, x_dim * y_dim * max_hgt * n_temp_metrics),
                     dim = c(x_dim, y_dim, max_hgt, n_temp_metrics))

  # Create grid of coordinates
  coords <- expand.grid(x = 1:x_dim, y = 1:y_dim)

  if (use_parallel) {
    # Set up parallel processing
    if (is.null(n_cores)) {
      n_cores <- min(detectCores() - 1, nrow(coords))
    }

    cat("Processing", nrow(coords), "cells using", n_cores, "cores\n")

    # Process cells in parallel
    results <- mclapply(1:nrow(coords), function(i) {
      process_single_cell(coords$x[i], coords$y[i], timestamps,
                         x_dim, y_dim, max_hgt, n_temp_metrics)
    }, mc.cores = n_cores)

  } else {
    # Sequential processing with progress tracking
    cat("Processing", nrow(coords), "cells sequentially\n")

    results <- vector("list", nrow(coords))
    for (i in 1:nrow(coords)) {
      results[[i]] <- process_single_cell(coords$x[i], coords$y[i], timestamps,
                                         x_dim, y_dim, max_hgt, n_temp_metrics)
      if (i %% 100 == 0) {
        # Print progress every 10 cells
        cat("Processed", i, "of", nrow(coords), "cells\n")
      }
    }
  }

  # Populate the main matrix with results
  total_time <- 0
  successful_cells <- 0

  for (result in results) {
    if (!is.null(result)) {
      x <- result$x
      y <- result$y
      actual_heights <- nrow(result$data)

      # Fill the matrix (up to the actual number of heights)
      mc_matrix[x, y, 1:actual_heights, ] <- result$data

      total_time <- total_time + result$processing_time
      successful_cells <- successful_cells + 1

      if (successful_cells %% 100 == 0) {
        avg_time <- total_time / successful_cells
        cat("Processed", successful_cells, "cells. Average time per cell:",
            round(avg_time, 3), "seconds\n")
      }
    }
  }

  cat("Processing complete! Total successful cells:", successful_cells, "\n")
  cat("Total processing time:", round(total_time, 2), "seconds\n")
  cat("Average time per cell:", round(total_time / successful_cells, 3), "seconds\n")

  return(mc_matrix)
}

# Usage examples:
# Sequential processing:
#mc_matrix <- process_microclimate_data(x_dim = 50, y_dim = 50, use_parallel = FALSE)

# Parallel processing (recommended):
mc_matrix <- process_microclimate_data(x_dim = 50, y_dim = 50, use_parallel = TRUE)

# For testing with smaller subset:
#mc_matrix <- process_microclimate_data(x_dim = 5, y_dim = 5, use_parallel = TRUE)

# Save the result
saveRDS(mc_matrix, "/Users/johanna/Uni/masterarbeit/data/mc_output/v2_mc_matrix.rds")

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Define variable names and units for better labeling
var_names <- c(
  "Mean Annual Temperature (°C)",
  "Mean Diurnal Temperature Range (°C)",
  "Annual Temperature Range (°C)",
  "Max Temperature of Warmest Month (°C)",
  "Min Temperature of Coldest Month (°C)",
  "Isothermality (%)",
  "Mean Annual Relative Humidity (%)",
  "Annual Humidity Range (%)",
  "Max Humidity (%)",
  "Min Humidity (%)",
  "Mean Annual Wind Speed (m/s)",
  "Annual Wind Speed Range (m/s)",
  "Max Wind Speed (m/s)",
  "Min Wind Speed (m/s)"
)

# Convert height index to actual height in meters
# Height 1 = 0.5m, then increases by 1m: 0.5, 1.5, 2.5, ..., 58.5
heights <- seq(0.5, 58.5, by = 1)

# Create a function to plot individual vertical profiles
plot_vertical_profile <- function(var_index, x_coord, y_coord) {
  # Extract data for this variable
  values <- mc_profile[x_coord, y_coord, 1:length(heights), var_index]

  # Create data frame
  profile_data <- data.frame(
    height = heights,
    value = values
  )

  # Create the plot
  ggplot(profile_data, aes(x = value, y = height)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 1.5, alpha = 0.7) +
    labs(
      title = var_names[var_index],
      x = "Value",
      y = "Height (m)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    ) +
    scale_y_continuous(breaks = seq(0, 60, by = 10))
}

# Create all 14 plots (assuming x and y coordinates - you'll need to specify these)
# For demonstration, I'll use x = 25, y = 25 (middle of your 50x50 grid)
x_coord <- 1  # Change this to your desired x coordinate
y_coord <- 1  # Change this to your desired y coordinate

mc_profile <- array(rep(NA, 1 * 1 * 80 * n_temp_metrics),
                     dim = c(1, 1, 80, n_temp_metrics))
mc_profile[1, 1, , ] <- apply(mc_matrix, c(3, 4), mean, na.rm = TRUE)

# Generate all plots
plot_list <- list()
for (i in 1:14) {
  plot_list[[i]] <- plot_vertical_profile(i, x_coord, y_coord)
}

# Option 1: Display all plots in a grid (might be crowded)
grid_plot <- grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)

# Option 2: Create separate plots for different variable groups
# Temperature variables (1-6)
temp_plots <- grid.arrange(grobs = plot_list[1:6], ncol = 3, nrow = 2,
                          top = "Temperature Variables - Vertical Profiles")

# Humidity variables (7-10)
humidity_plots <- grid.arrange(grobs = plot_list[7:10], ncol = 2, nrow = 2,
                              top = "Humidity Variables - Vertical Profiles")

# Wind speed variables (11-14)
wind_plots <- grid.arrange(grobs = plot_list[11:14], ncol = 2, nrow = 2,
                          top = "Wind Speed Variables - Vertical Profiles")

# Option 3: Create a combined data frame for faceted plotting
create_combined_data <- function(x_coord, y_coord) {
  combined_data <- data.frame()

  for (i in 1:14) {
    values <- mc_profile[x_coord, y_coord, 1:max_hgt, i]
    temp_data <- data.frame(
      height = heights,
      value = values,
      variable = var_names[i],
      var_group = case_when(
        i <= 6 ~ "Temperature",
        i <= 10 ~ "Humidity",
        TRUE ~ "Wind Speed"
      )
    )
    combined_data <- rbind(combined_data, temp_data)
  }
  return(combined_data)
}

# Create faceted plot
combined_data <- create_combined_data(x_coord, y_coord)

faceted_plot <- ggplot(combined_data, aes(x = value, y = height)) +
  geom_line(color = "blue", size = 0.8) +
  geom_point(color = "red", size = 1, alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_x", ncol = 4) +
  labs(
    title = paste("Vertical Profiles of Microclimate Variables (Grid Point:", x_coord, ",", y_coord, ")"),
    x = "Value",
    y = "Height (m)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  ) +
  scale_y_continuous(breaks = seq(0, 60, by = 20))

print(faceted_plot)

pdf("../../figs/mc_output/mc_plot_avg_metrics.pdf")
print(faceted_plot)
dev.off()