
library(micropoint)
library(terra)
library(readr)
library(viridis)
library(microclimf)
library(lubridate)
library(furrr)
library(future)
library(microclimf)
library(dplyr)
library(ggplot2)

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
  point <- vect(cbind(lon, lat), crs = crs)

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

coords2indices <- function(lon, lat, raster, coord_crs = "EPSG:4326") {
  # Create a point in WGS84
  point <- terra::vect(cbind(lon, lat), crs = coord_crs)
  point_proj <- terra::project(point, terra::crs(raster)) # match raster crs

  # Get the cell number from the transformed point
  cell <- terra::cellFromXY(raster,
                            t(as.matrix(terra::geom(point_proj)[, c("x", "y")])))

  # Convert the cell number to row and column
  rc <- terra::rowColFromCell(raster, cell)
  names(rc) <- c("x", "y")
  rc <- rc[c("x", "y")]

  return(rc)
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


# vegp 1m resolution
# dtm 1m resolution
# soilc 1m resolution
# climdata is equal for the whole grid

# Function to process a single (x, y) cell
x <- 25
y <- 25

# Load data for one year
in_dir <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
vegp_reg <- readRDS(paste(in_dir, "vegp_mof3d_ptm_v2.RDS", sep = "/"))
soilc_reg <- readRDS(paste(in_dir, "soilc_v2.RDS", sep = "/"))
climdata_reg <- read_csv(paste(in_dir, "era5_climdata_2024_v2.csv", sep = "/"))

# Get PAI
microhab_file <- "/Users/johanna/Uni/masterarbeit/code/output/modev_zach_25_01_07/MicrohabitatMatrix98.rds"
pai <- readRDS(microhab_file)[,,,5]
paii <- pai[x, y, 1:max(vegparams$h, 0.5)]

# Get coordiantes
coords <- indices2coords(x, y, terra::unwrap(vegp_reg$pai))[c("x", "y")]
lon <- coords[[1]]
lat <- coords[[2]]

# Get params for the point
vegparams <- extract_params(vegp_reg, lon, lat)
grndparams <- extract_params(soilc_reg, lon, lat)

# Veg heights
max_veg_height <- max(terra::values(terra::unwrap(vegp_reg$h)), na.rm = TRUE)
heights <- seq(0.5, max_veg_height + 1)
i <- 2
h <- heights[i]

mout <- micropoint::runpointmodel(climdata_reg, reqhgt = h, vegparams,
                                  paii, grndparams, lat = lat, long = lon)

# put into dataframe
mc_sim <- data.frame(
  obs_time = mout$obs_time,
  tair = mout$tair,
  tleaf = mout$tleaf,
  relhum = mout$relhum
)

# Load empirical data
emp_path <- "/Users/johanna/Uni/masterarbeit/data/empirical/Datalogger 400m elevation REGUA understory Trilha Verde"
emp_dir <- "3600m, 446mASL"
emp_file <- "21983350 2024-11-06 21_41_33 Brazil Standard Time (Data Brazil Standard Time).xlsx"

emp_data <- readxl::read_excel(file.path(emp_path, emp_dir, emp_file))

# Format the empirical data and join with simulation data
emp_sim_data <- emp_data %>%
  select("Date-Time (Brazil Standard Time)", "Temperature (°C)", "RH (%)") %>%
  rename("obs_time" = "Date-Time (Brazil Standard Time)",
         "tair_emp" = "Temperature (°C)",
         "relhum_emp" = "RH (%)") %>%
    mutate(
    # Parse obs_time as POSIXct in Brazil Standard Time (UTC-3)
    obs_time = force_tz(as.POSIXct(obs_time), tzone = "America/Sao_Paulo"),

    # Convert to UTC
    obs_time_utc = with_tz(obs_time, tzone = "UTC"),

    # Floor to the hour
    obs_hour_utc = floor_date(obs_time_utc, unit = "hour")
  ) %>%
  # Group by floored hourly time and average the values
  group_by(obs_hour_utc) %>%
  summarise(
    tair_emp = mean(tair_emp, na.rm = TRUE),
    relhum_emp = mean(relhum_emp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(obs_time = obs_hour_utc) %>%
  # Join with the microclimate simulation data
  left_join(., mc_sim, by = "obs_time")

# Assuming emp_sim_data is already loaded and is a tibble

# First plot: Air temperature (empirical vs. simulated)
plot_airt <- ggplot(emp_sim_data, aes(x = obs_time)) +
  geom_line(aes(y = tair_emp, color = "Empirical")) +
  geom_line(aes(y = tair, color = "Simulated")) +
  labs(
    title = "Air Temperature: Empirical vs. Simulated",
    x = "Time",
    y = "Air Temperature (°C)",
    color = "Measurement"
  ) +
  theme_minimal()

# Second plot: Relative Humidity (empirical vs. simulated)
plot_relhum <- ggplot(emp_sim_data, aes(x = obs_time)) +
  geom_line(aes(y = relhum_emp, color = "Empirical")) +
  geom_line(aes(y = relhum, color = "Simulated")) +
  labs(
    title = "Relative Humidity: Empirical vs. Simulated",
    x = "Time",
    y = "Relative Humidity (%)",
    color = "Measurement"
  )
  theme_minimal()

# Print plots to a pdf file
pdf("../../figs/mc_output/airt_emp_vs_sim_mc_regua_3600m_446mASL.pdf")
print(plot_airt)
dev.off()

pdf("../../figs/mc_output/relhum_emp_vs_sim_mc_regua_3600m_446mASL.pdf")
print(plot_relhum)
dev.off()