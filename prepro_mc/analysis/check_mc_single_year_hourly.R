library(microclimf)
library(terra)
library(readr)
library(viridis)
library(ggplot2)
library(lubridate)
library(reshape2)
library(dplyr)

# Function to load tair and average over time
load_avg_tair <- function(x, y, base_path, var = "tair") {
  # Map var to index tair -> 1, relhum -> 7, windspeed -> 11
  idx <- switch(var,
              tair = 1,
              relhum = 7,
              windspeed = 11)
  file <- sprintf("%s/mc_x%d_y%d.rds", base_path, x, y)
  data <- readRDS(file)
  return(data$data[, idx])
}

# Settings
base_path <- "/Users/johanna/Uni/masterarbeit/data/mc_output/v5/regua"
#z_levels <- 59
x_range <- c(1)
y_range <- 1:13

# Sample random slices
set.seed(42)
random_x <- sample(x_range, 1)  # 3 random X values for YZ slices
random_y <- sample(y_range, 3)  # 3 random Y values for XZ slices

# -------------- Check MC in XZ and YZ slices --------------

for (var in c("tair", "relhum", "windspeed")) {
  # Gather data for XZ slices
  xz_data <- do.call(rbind, lapply(random_y, function(y) {
    do.call(rbind, lapply(x_range, function(x) {
      tair_avg <- load_avg_tair(x, y, base_path, var)
      z_levels <- length(tair_avg)
      data.frame(x = x, z = 1:z_levels, value = tair_avg, slice = paste("y =", y), variable = var)
    }))
  }))

  # Gather data for YZ slices
  yz_data <- do.call(rbind, lapply(random_x, function(x) {
    do.call(rbind, lapply(y_range, function(y) {
      tair_avg <- load_avg_tair(x, y, base_path, var)
      z_levels <- length(tair_avg)
      data.frame(y = y, z = 1:z_levels, value = tair_avg, slice = paste("x =", x), variable = var)
    }))
  }))

  # Plot XZ slices
  ggplot(xz_data, aes(x = x, y = z, fill = value)) +
    geom_tile() +
    scale_fill_viridis(option = "viridis", name = paste("Avg", var)) +
    facet_wrap(~ slice + variable, ncol = 1) +
    labs(title = paste("XZ Slices of Average", var), x = "X", y = "Height (Z)") +
    theme_minimal()
  ggsave(paste0("../../figs/mc_output/", "v5_xz_slices_", var, ".png"))

  # Plot YZ slices
  ggplot(yz_data, aes(x = y, y = z, fill = value)) +
    geom_tile() +
    scale_fill_viridis(option = "magma", name = paste("Avg", var)) +
    facet_wrap(~ slice + variable, ncol = 1) +
    labs(title = paste("YZ Slices of Average", var), x = "Y", y = "Height (Z)") +
    theme_minimal()
  ggsave(paste0("../../figs/mc_output/", "v5_yz_slices_", var, ".png"))
}

# -------------- Check time series of avg MC for this one year --------------

# Load all data and convet from hourly to mean daily values (50 x 50 x 59 x 365) for each variable
dates_2024 <- seq.Date(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")

x_range <- 1:50
y_range <- 1:50
z_levels <- dim(readRDS("/Users/johanna/Uni/masterarbeit/data/mc_output/v5/regua/mc_x1_y1.rds")$data)[1]
vars <- c("tair", "relhum")

# Initialize storage
sim <- array(NA, dim = c(50, 50, z_levels, 2))

# Load all data
for (x in x_range) {
  for (y in y_range) {
    file <- sprintf("%s/mc_x%d_y%d.rds", base_path, x, y)
    if (!file.exists(file)) {
      print(paste0("File does not exist, skipping...", file))
      next  # Skip if file does not exist
    }
    data <- readRDS(file)
    sim[x, y, , ] <- data$data[, c(1, 7)]
  }
}

# Compute average profiles
avg_tair   <- apply(sim[,,,1], 3, mean, na.rm = TRUE)
avg_relhum <- apply(sim[,,,2], 3, mean, na.rm = TRUE)

# Get dimensions
nx <- dim(sim)[1]
ny <- dim(sim)[2]
nz <- dim(sim)[3]

# Heights (replace with actual height vector if available)
heights <- seq(0.5, nz - 0.5, by = 1)

# RELATIVE HUMIDITY PLOT
pdf("../../figs/mc_output/regua_relhum_24_v5.pdf")
plot(avg_relhum, heights, type = "n",
     xlab = "Average Relative Humidity (%)",
     ylab = "Height (m)")

# Loop over all (x, y) points
for (ix in 1:nx) {
  for (iy in 1:ny) {
    lines(sim[ix, iy, , 2], heights,
          col = rgb(0.5, 0.5, 0.5, alpha = 0.05), lwd = 0.5)
  }
}

# Average on top
lines(avg_relhum, heights, col = "blue", lwd = 2)
dev.off()

# TEMPERATURE PLOT
pdf("../../figs/mc_output/regua_tair_24_v5.pdf")
plot(avg_tair, heights, type = "n",
     xlab = "Average Temperature (°C)",
     ylab = "Height (m)")

for (ix in 1:nx) {
  for (iy in 1:ny) {
    lines(sim[ix, iy, , 1], heights,
          col = rgb(0.5, 0.5, 0.5, alpha = 0.05), lwd = 0.5)
  }
}

lines(avg_tair, heights, col = "blue", lwd = 2)
dev.off()