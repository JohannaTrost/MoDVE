require(devtools)
# install_github("ilyamaclean/microclimf")

library(microclimf)
library(terra)
library(readr)
library(viridis)
library(suncalc)
library(dplyr)
library(abind)


agg_model_var <- function(data, time_vec, time_df) {
  
  # Reshape mout$Tz from 3D to 2D for easier processing
  mat <- matrix(data, nrow = 50*50, ncol = 8784)
  
  # Create vectors of day index (1 to 366), and logicals for day/night
  day_idx <- as.numeric(factor(as.Date(time_vec)))
  is_day <- time_df$period == "day"
  is_night <- time_df$period == "night"
  
  # Initialize result arrays
  mean_daily <- array(NA, dim = c(50*50, 366))
  mean_daytime <- array(NA, dim = c(50*50, 366))
  mean_nighttime <- array(NA, dim = c(50*50, 366))
  
  # Loop over each day
  for (d in 1:366) {
    idx_day <- which(day_idx == d & is_day)
    idx_night <- which(day_idx == d & is_night)
    idx_all <- which(day_idx == d)
    
    if (length(idx_all) > 0) {
      mean_daily[, d] <- rowMeans(mat[, idx_all], na.rm = TRUE)
    }
    if (length(idx_day) > 0) {
      mean_daytime[, d] <- rowMeans(mat[, idx_day], na.rm = TRUE)
    }
    if (length(idx_night) > 0) {
      mean_nighttime[, d] <- rowMeans(mat[, idx_night], na.rm = TRUE)
    }
  }
  
  # Reshape
  mean_daily_arr <- array(mean_daily, dim = c(50, 50, 366))
  mean_daytime_arr <- array(mean_daytime, dim = c(50, 50, 366))
  mean_nighttime_arr <- array(mean_nighttime, dim = c(50, 50, 366))
  
  return(list(
    mean_daily = mean_daily_arr,
    mean_daytime = mean_daytime_arr,
    mean_nighttime = mean_nighttime_arr
  ))
}


# vegp 1m resolution
# dtm 1m resolution
# soilc 1m resolution
# climdata is equal for the whole grid 

out = "/Users/johanna/Uni/masterarbeit/code/output/microclimf/all_days_sim_regua_2024"

# Load data for one year 
in_dir <- "/Users/johanna/Uni/masterarbeit/code/data/mc_input/regua"
vegp_reg <- readRDS(paste(in_dir, "vegp.RDS", sep = "/"))
dtm_reg <- rast(paste(in_dir, "dtm.tif", sep = "/"))
soilc_reg <- readRDS(paste(in_dir, "soilc.RDS", sep = "/"))
climdata_reg <- read_csv(paste(in_dir, "era5_climdata_2024.csv", sep = "/"))

max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$hgt)), na.rm = TRUE)
min_veg_height <- min(terra::values(terra::unwrap(vegp_reg$hgt)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)
n_heights <- length(heights)

# First run to initialize output structure with dimensions
height <- heights[1]
micropoint <- runpointmodel(climdata_reg, reqhgt = height, dtm_reg, vegp_reg, soilc_reg)

# --- Get nighttime and daytime
time_vec <- as.POSIXct(micropoint$weather$obs_time)

# Create a sequence of days for the year
days <- unique(as.Date(time_vec))

# Location
lat <- -22.426880
lon <- -42.765096

# Get sunrise and sunset for each day
sun_times <- getSunlightTimes(date = days, lat = lat, lon = lon, keep = c("sunrise", "sunset"), tz = "UTC")

# Create a data.frame for hourly times
time_df <- data.frame(datetime = time_vec)
time_df$date <- as.Date(time_df$datetime)

# Join with sun_times to get daily sunrise/sunset
time_df <- left_join(time_df, sun_times, by = c("date" = "date"))

# Classify as daytime or nighttime
time_df$period <- ifelse(
  time_df$datetime >= time_df$sunrise & time_df$datetime < time_df$sunset,
  "day", "night"
)
# ---

mean_daily <- list()
mean_daytime <- list()
mean_nighttime <- list()

mout <- runmicro(micropoint, reqhgt = height, vegp_reg, soilc_reg, dtm_reg,
                 method = "Cpp")
saveRDS(mout, paste(out, "/", height, "m_mout.rds", sep = ""))

# Preallocate 4D arrays for each variable
for (var in c("Tz")) {
  gc()
  res <- agg_model_var(mout[[var]], time_vec, time_df)
  dims <- dim(res$mean_daily)
  
  # Save res
  mean_daily[[var]] <- array(NA, dim = c(dims[1], dims[2], dims[3], n_heights))
  mean_daily[[var]][,,,1] <- res$mean_daily
  mean_daytime[[var]] <- array(NA, dim = c(dims[1], dims[2], dims[3], n_heights))
  mean_daytime[[var]][,,,1] <- res$mean_daytime
  mean_nighttime[[var]] <- array(NA, dim = c(dims[1], dims[2], dims[3], n_heights))
  mean_nighttime[[var]][,,,1] <- res$mean_nighttime
}

# Loop through remaining heights
for (i in 2:n_heights) {
  height <- heights[i]
  
  # Run point and grid models
  micropoint <- runpointmodel(climdata_reg, reqhgt = height, dtm_reg, vegp_reg, soilc_reg)
  mout <- runmicro(micropoint, reqhgt = height, vegp_reg, soilc_reg, dtm_reg,
                   method = "Cpp")
  saveRDS(mout, paste(out, "/", height, "m_mout.rds", sep = ""))
  
  # Store in 4D arrays
  for (var in c("Tz")) {
    gc()
    res <- agg_model_var(mout[[var]], time_vec, time_df)
    
    # Save res
    mean_daily[[var]][,,,i] <- res$mean_daily
    mean_daytime[[var]][,,,i] <- res$mean_daytime
    mean_nighttime[[var]][,,,i] <- res$mean_nighttime
  }
  
  cat("Stored height index", i, "for height", height, "\n")
}

# Save simulation outputs
saveRDS(mean_daily[["Tz"]], paste(out, "mean_daily_tz.rds", sep = "/"))
saveRDS(mean_daytime[["Tz"]], paste(out, "mean_daytime_tz.rds", sep = "/"))
saveRDS(mean_nighttime[["Tz"]], paste(out, "mean_nighttime_tz.rds", sep = "/"))

