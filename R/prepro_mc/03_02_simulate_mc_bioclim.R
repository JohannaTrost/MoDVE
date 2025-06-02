require(devtools)
# install_github("ilyamaclean/microclimf")

library(microclimf)
library(terra)
library(readr)
library(viridis)

# vegp 1m resolution
# dtm 1m resolution
# soilc 1m resolution
# climdata is equal for the whole grid 


# Load data for one year 
in_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
vegp_reg <- readRDS(paste(in_dir, "vegp.RDS", sep = "/"))
dtm_reg <- rast(paste(in_dir, "dtm.tif", sep = "/"))
soilc_reg <- readRDS(paste(in_dir, "soilc.RDS", sep = "/"))
climdata_reg <- read_csv(paste(in_dir, "era5_climdata_2024.csv", sep = "/"))

max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$hgt)), na.rm = TRUE)
min_veg_height <- min(terra::values(terra::unwrap(vegp_reg$hgt)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)
n_heights <- length(heights)

# Initialize an empty list to store output variables
avg_mc_all <- list()

# First run to initialize output structure with dimensions
height <- heights[1]
mout <- runbioclim(climdata_reg, reqhgt = height, vegp_reg, soilc_reg, dtm_reg)

# Determine dimensions (assumes all variables are same size)
dims <- dim(mout$bio1)  # e.g., [50, 50, 288]

# Preallocate 4D arrays for each variable
for (var in names(mout)) {
  avg_mc_all[[var]] <- array(NA, dim = c(dims[1], dims[2], dims[3], n_heights))
  avg_mc_all[[var]][,,,1] <- as.array(mout[[var]])
}

# Loop through remaining heights
for (i in 2:n_heights) {
  height <- heights[i]
  
  # Run model
  mout <- runbioclim(climdata_reg, reqhgt = height, vegp_reg, soilc_reg, dtm_reg)

  # Store in 4D arrays
  for (var in names(mout)) {
    avg_mc_all[[var]][,,,i] <- as.array(mout[[var]])
  }
  
  cat("Stored height index", i, "for height", height, "\n")
}

# -------
# Plot grid of 3 heights
# Plot air temperatures on hottest hour in micropoint (2017-06-20 13:00:00 UTC)
mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "orange",  "red"))(255)

for (i in seq(1, 288, 50)) {
  
  lower <- min(avg_mc_all$Tz[,,i, c(1, 12, 25)])
  upper <- max(avg_mc_all$Tz[,,i, c(1, 12, 25)])
  
  plot(rast(avg_mc_all$Tz[,,i, 1]), col = mypal, 
       range = c(lower, upper),
       main = "Mean Tz 0.5m")
  plot(rast(avg_mc_all$Tz[,,i, 12]), col = mypal, 
       range = c(lower, upper),
       main = "Mean Tz 11.5m")
  plot(rast(avg_mc_all$Tz[,,i, 25]), col = mypal, 
       range = c(lower, upper),
       main = "Mean Tz 24.5m")
}

# Plot avg profiles
for (var in names(avg_mc_all)) {
  avg_profile <- apply(avg_mc_all[[var]], 4, mean)
  plot(avg_profile, heights, type="b", 
       main = var)
}



# -------
# Show example of a slice 

cols <- viridis(10)

# Height 
hgt_raster <- terra::unwrap(vegp_reg$hgt)
height_means <- rowMeans(hgt_raster)

# Loop over time indices
for (i in seq(1, 288, 20)) {
  
  all_means <- list()
  all_heights <- list()
  
  # First loop: collect all profiles for this time index
  for (y_index in seq(1, 50, 5)) {
    xz_slice <- avg_mc_all$Tz[, y_index, i, ]
    mean_temp_by_height <- colMeans(xz_slice, na.rm = TRUE)
    xz_slice_rast <- rast(t(xz_slice))
    height_vals <- yFromRow(xz_slice_rast, nrow(xz_slice_rast):1)
    
    all_means[[length(all_means) + 1]] <- mean_temp_by_height
    all_heights[[length(all_heights) + 1]] <- height_vals
  }
  
  # Get global plot limits
  all_means_combined <- unlist(all_means)
  all_heights_combined <- unlist(all_heights)
  xlim_vals <- range(all_means_combined, na.rm = TRUE)
  ylim_vals <- range(all_heights_combined, na.rm = TRUE)
  
  # Start plot
  plot(NULL, xlim = xlim_vals, ylim = ylim_vals,
       xlab = "Average Temperature (°C)", ylab = "Height (m)",
       main = paste("Average Tz Profiles\n Day =", i))
  
  for (k in seq_along(all_means)) {
    print(k)
    col <- cols[k]
    lines(all_means[[k]], all_heights[[k]], col = col, lwd = 1.5, alpha=0.5)
    points(all_means[[k]], all_heights[[k]], col = col, pch = 16, cex = 0.5, , alpha=0.5)
    abline(h=row_means[k*5], col=col, lty=2, lwd = 1.5, alpha=0.5)
  }
}

