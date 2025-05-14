require(devtools)
install_github("ilyamaclean/microclimf")

library(microclimf)
library(terra)

# vegp 1m resolution
# dtmcaerth 1m resolution
# soilc 1m resolution
# climdata is equal for the whole grid 


# Load data for one year 
# (...)

max_veg_height <- max(terra::values(terra::unwrap(vegp$hgt)), na.rm = TRUE)
min_veg_height <- min(terra::values(terra::unwrap(vegp$hgt)), na.rm = TRUE)
heights <- seq(min_veg_height + 0.5, max_veg_height + 1)
n_heights <- length(heights)

# Initialize an empty list to store output variables
avg_mc_all <- list()

# First run to initialize output structure with dimensions
height <- heights[1]
micropoint <- runpointmodel(climdata, reqhgt = height, dtmcaerth, vegp, soilc)
micropoint_mx <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
micropoint_mn <- subsetpointmodel(micropoint, tstep = "month", what = "tmin")
mout_mx <- runmicro(micropoint_mx, reqhgt = height, vegp, soilc, dtmcaerth)
mout_mn <- runmicro(micropoint_mn, reqhgt = height, vegp, soilc, dtmcaerth)

# Average the two outputs
avg_mc <- list()
for (var in names(mout_mx)) {
  avg_mc[[var]] <- (mout_mx[[var]] + mout_mn[[var]]) / 2
}

# Determine dimensions (assumes all variables are same size)
dims <- dim(avg_mc$Tz)  # e.g., [50, 50, 288]

# Preallocate 4D arrays for each variable
for (var in names(avg_mc)) {
  avg_mc_all[[var]] <- array(NA, dim = c(dims[1], dims[2], dims[3], n_heights))
  avg_mc_all[[var]][,,,1] <- avg_mc[[var]]
}

# Loop through remaining heights
for (i in 2:n_heights) {
  height <- heights[i]
  
  # Run point and grid models
  micropoint <- runpointmodel(climdata, reqhgt = height, dtmcaerth, vegp, soilc)
  micropoint_mx <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
  micropoint_mn <- subsetpointmodel(micropoint, tstep = "month", what = "tmin")
  mout_mx <- runmicro(micropoint_mx, reqhgt = height, vegp, soilc, dtmcaerth)
  mout_mn <- runmicro(micropoint_mn, reqhgt = height, vegp, soilc, dtmcaerth)
  
  # Average variables
  avg_mc <- list()
  for (var in names(mout_mx)) {
    avg_mc[[var]] <- (mout_mx[[var]] + mout_mn[[var]]) / 2
  }
  
  # Store in 4D arrays
  for (var in names(avg_mc)) {
    avg_mc_all[[var]][,,,i] <- avg_mc[[var]]
  }
  
  cat("Stored height index", i, "for height", height, "\n")
}

# -------
# Plot grid of 3 heights
# Plot air temperatures on hottest hour in micropoint (2017-06-20 13:00:00 UTC)
mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow", "orange",  "red"))(255)
plot(rast(avg_mc_all$Tz[,,134, 1]), col = mypal, range = c(20, 48),
     main = "Mean Tz 1m")
plot(rast(avg_mc_all$Tz[,,134, 2]), col = mypal, range = c(20, 48),
     main = "Mean Tz 2m")
plot(rast(avg_mc_all$Tz[,,134, 3]), col = mypal, range = c(20, 48),
     main = "Mean Tz 3m")

# -------
# Show example of a slice 

# Choose y index and time index
y_index <- 25
time_index <- 134

# Extract [x, z] matrix at fixed y and time
xz_slice <- avg_mc_all$Tz[, y_index, time_index, ]

# Transpose for raster format: [rows = height, columns = x]
xz_slice_rast <- rast(t(xz_slice))

# Plot
plot(xz_slice_rast, col = mypal, range = c(20, 48),
     main = "Tz cross-section (x vs height)",
     xlab = "x coordinate", ylab = "Height (m)")



