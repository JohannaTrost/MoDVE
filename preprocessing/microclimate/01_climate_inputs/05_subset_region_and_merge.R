# ------
# Subset the processed data to a 50m x 50m grid and format

library(sf)
library(terra)
library(dplyr)
library(microclimf)
library(readr)
library(ncdf4)
library(tidyr)
library(circular)
library(lubridate)
library(tidyr)

# GLOBAL variables (to configure)
REGION <- "regua"  # set to "regua" or "pirineus"

# Directories
data_dir <- "../modve_data"

# Add a diurnal cycle to temperature and relative humidity in CMIP6 from era5 anomaly
add_era5_diurnal_cycle <- function (climdata_regua, cmip6_hourly, cmip6_year) {

  # 1. Compute daily anolmaly for Era5 data
  climdata_regua_anomaly <- climdata_regua %>%
    mutate(
      day = as.Date(obs_time),
    ) %>%
    group_by(day) %>%
    mutate(
      temp_anom = temp - mean(temp, na.rm = TRUE),
      relhum_anom = relhum - mean(relhum, na.rm = TRUE),
      windspeed_anom = windspeed - mean(windspeed, na.rm = TRUE),
    ) %>%
    ungroup()

  # 2. Equalize number of days in CMIP6 and ERA5
  # Count rows to detect leap vs non-leap
  n_clim  <- nrow(climdata_regua_anomaly)
  n_cmip6 <- nrow(cmip6_hourly)

  # Aligning
  if (n_cmip6 == 8760 && n_clim == 8784) {
    # Drop Feb 29
    clim_aligned <- climdata_regua_anomaly %>%
      filter(!(month(obs_time) == 2 & day(obs_time) == 29)) %>%
      arrange(obs_time)
  } else if (n_cmip6 == 8784 && n_clim == 8760) {
    # Insert Feb 29
    clim_aligned <- insert_feb29(climdata_regua_anomaly, year = cmip6_year)
  } else {
    if (n_clim != n_cmip6){
      cat("CMIP6 != ERA5, ERA5: ", n_clim, "CMIP6: ", n_cmip6)
    }
    clim_aligned <- climdata_regua_anomaly
  }
  return(clim_aligned)
}

# Convert daily CMIP6 data to hourly splitting precipitation
cmip62hourly <- function (cmip6_long) {
  cmip6_long %>%
    # Create 24 hours for each day
    slice(rep(1:n(), each = 24)) %>%
    # Add hour to each day
    mutate(
      hour = rep(0:23, times = nrow(cmip6_long)),
      obs_time = obs_time + hours(hour)
    ) %>%
    # Rescale variables appropriately
    mutate(
      # DIVIDE by 24 (cumulative → rate):
      precip = precip / 24,      # Daily total → hourly rate

      # REPEAT (instantaneous values):
      pres = pres,               # Pressure stays same
      relhum = relhum,           # Relative humidity stays same
      temp = temp,               # Temperature stays same
      windspeed = windspeed      # Wind speed stays same
    ) %>%
    select(-hour)
}

# Load cmip6 data for location and generate long format df
extract_cmip6 <- function (data_dir, cmip6_year, lat, lon) {

  if (cmip6_year < 2015) {
    experimt <- "historical"
  } else {
    experimt <- "ssp245"
  }

  # Load CMIP6 data
  cmip_baf <- rast(
    file.path(data_dir, "mc_input", "climate", "cmip6_ceda",
              paste0("baf_ensemble_day_", experimt, "_", cmip6_year, ".nc"))
  )

  # Get the coordinates of that cell (i.e., the raster grid point)
  nearest_xy <- xyFromCell(cmip_baf, cellFromXY(cmip_baf, matrix(c(lat, lon), ncol = 2)))

  # Extract data at that cell
  extracted_cmip6 <- terra::extract(cmip_baf, nearest_xy)

  # Convert to long df
  cmip6_long <- extracted_cmip6 %>%
    pivot_longer(
      everything(),
      names_to = c("variable", "day"),
      names_pattern = "(.+)_(\\d+)",
      values_to = "value"
    ) %>%
    mutate(day = as.integer(day)) %>%
    arrange(variable, day) %>%
    group_by(variable) %>%
    mutate(obs_time = seq(as.Date(paste0(cmip6_year, "-01-01")), by = "1 day", length.out = n())) %>%
    ungroup() %>%
    select(obs_time, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)

  return(cmip6_long)
}

# Load era5 data for location and load into df
extract_era5 <- function(data_dir, era5_year, lon, lat) {
  climdata_baf <- readRDS(
    file.path(data_dir, "mc_input", "climate", "era5_processed", paste0("era5_climdata_", era5_year, ".RDS"))
  )

  # Unwrap all PackedSpatRasters in the list
  climdata_baf_unpacked <- lapply(climdata_baf, terra::unwrap)

  # Define REGUA Point
  xy <- matrix(c(lon, lat), ncol = 2)
  p <- vect(xy, crs = "EPSG:4326")

  # 3. Extract time values from any one raster (assumes all rasters are time-aligned)
  layer_names <- names(climdata_baf_unpacked$precip)
  timestamps <- sub(".*=", "", layer_names)
  timestamps <- as.POSIXct(as.numeric(timestamps), origin = "1970-01-01", tz = "UTC")

  # 4. Extract values from each variable and build named list
  extracted_values <- lapply(climdata_baf_unpacked, function(r) {
    vals <- terra::extract(r, p)[1, -1, drop = FALSE]  # drop ID column
    as.numeric(vals)
  })

  # 5. Combine into final data frame
  climdata_regua <- data.frame(
    obs_time = timestamps,
    extracted_values,
    check.names = FALSE  # keep original variable names
  )

  # Clean data
  climdata_regua$swdown[climdata_regua$swdown < 0] <- 0
  climdata_regua$difrad[climdata_regua$difrad < 0] <- 0

  return(climdata_regua)
}

# Function to insert Feb 29 by extrapolation
insert_feb29 <- function(df, year) {
  # Build Feb 29 sequence
  feb29_times <- seq(ymd_hms(paste0(year, "-02-29 00:00:00")),
                     ymd_hms(paste0(year, "-02-29 23:00:00")), by = "1 hour")

  # Use Feb 28 and Mar 1 as neighbors
  feb28 <- df %>% filter(month(obs_time) == 2, day(obs_time) == 28)
  mar01 <- df %>% filter(month(obs_time) == 3, day(obs_time) == 1)

  # Simple average interpolation (between Feb 28 and Mar 1)
  feb29 <- feb28
  feb29$obs_time <- feb29_times
  num_cols <- sapply(feb29, is.numeric)
  feb29[num_cols] <- (feb28[num_cols] + mar01[num_cols]) / 2

  # Bind into full dataset
  bind_rows(df, feb29) %>% arrange(obs_time)
}


terra_crop <- function(data, ext, res = NULL) {
  data_unpacked <- lapply(data, terra::unwrap)

  if (!is.null(res)) {
    # Expand extent for buffer
    expanded_vect <- buffer(ext, width = 10)

    # Crop all rasters first to the same buffered extent
    data_unpacked <- lapply(data_unpacked, function(rast) {
      terra::crop(rast, expanded_vect)
    })

    # Define a common template from the buffered extent
    template <- terra::rast(extent = terra::ext(expanded_vect),
                            resolution = res,
                            crs = terra::crs(data_unpacked[[1]]))  # assume all same CRS

    # Resample all rasters to the common template
    data_unpacked <- lapply(data_unpacked, function(rast) {
      rast_names <- names(rast)
      is_cat <- any(rast_names %in% c("soiltype", "clump"))
      method <- if (is_cat) "near" else "bilinear"
      terra::resample(rast, template, method = method)
    })
  }

  # Final crop to the exact extent
  data_square <- lapply(data_unpacked, function(rast) {
    terra::crop(rast, ext)
  })

  data_square_packed <- lapply(data_square, terra::wrap)
  return(data_square_packed)
}

dms_to_decimal <- function(dms) {
  parts <- as.numeric(unlist(regmatches(dms, gregexpr("[0-9.]+", dms))))
  sign <- ifelse(grepl("[SW]", dms), -1, 1)
  return(sign * (parts[1] + parts[2] / 60 + parts[3] / 3600))
}

if (REGION == "regua") {
  # -- REGUA

  out <- file.path(data_dir, "mc_input", "regua")
  dir.create(file.path(out))

  # Define the center point - 22°23′44.75″S, 42°44′15.78″W
  lat <- dms_to_decimal("22°23′44.75″S")
  lon <- dms_to_decimal("42°44′15.78″W")
  center_coords <- st_sfc(st_point(c(lon, lat)), crs = 4326)
  center_proj <- st_transform(center_coords, crs = 31983)
  center_coords_m <- st_coordinates(center_proj)

} else if (REGION == "pirineus") {
  # -- Pirineus

  out <- file.path(data_dir, "mc_input", "pirineus")
  dir.create(file.path(out))

  # Define the center point
  lat <- dms_to_decimal("22°26′46.74″S")
  lon <- dms_to_decimal("42°30′06.16″W")
  center_coords <- st_sfc(st_point(c(lon, lat)), crs = 4326)
  center_proj <- st_transform(center_coords, crs = 31983)
  center_coords_m <- st_coordinates(center_proj)
}

# Create a 50m x 50m square (i.e., a 50m buffer in all directions)
extent_box <- st_buffer(center_proj, dist = 25, endCapStyle = "SQUARE")
extent_vect <- vect(extent_box)
# Ceil the extent if decimals
e <- ext(extent_vect)
ce <- ext(round(e[1]), round(e[2]), round(e[3]), round(e[4]))
# Create new polygon from ceiled extent
ceiled_box <- as.polygons(ce)
crs(ceiled_box) <- crs(extent_vect)  # Set CRS to match original

# - Prep climate data
# List of years for "aligning" ERA5 and CMIP6
era5_2025_2100 <- rep(seq(2020, 2024), 16)[((16*5)-75):(16*5)] # for 2025-2100 repeat 2020-2024 16 times
era5_yrs <- c(seq(1980, 2024), era5_2025_2100)
cmip6_yrs <- seq(1980, 2100)

for (i in seq_along(era5_yrs)) { # Download each year

  era5_year <- era5_yrs[i]
  cmip6_year <- cmip6_yrs[i]

  cat(paste0("Processing climate data for year ", cmip6_year, "\n"))

  # Prep ERA5 data
  climdata_regua <- extract_era5(data_dir, era5_year, lon, lat)

  # Prep CMIP6 data
  cmip6_long <- extract_cmip6(data_dir, cmip6_year, lon, lat)

  # - Rescale CMIP6 data from daily to hourly

  # Create hourly timestamps for each day
  cmip6_hourly <- cmip62hourly(cmip6_long)

  # Add a diurnal cycle to temperature and relative humidity
  clim_aliged <- add_era5_diurnal_cycle(climdata_regua, cmip6_hourly, cmip6_year)

  # Merge CMIP6 with remaining ERA5 data
  clim_aligned <- clim_aligned %>%
    mutate(obs_time = cmip6_hourly$obs_time) %>%
    select(-c(precip, temp, relhum, pres, windspeed)) %>%
    inner_join(., cmip6_hourly, by = "obs_time") %>%
      mutate(
          # Add diurnal cycle to temperature and relative humidity
          temp = temp + temp_anom,
          relhum = relhum + relhum_anom,
          windspeed = windspeed + windspeed_anom
      ) %>%
    select(-temp_anom, -relhum_anom, -day)

  # Save data
  write.csv(clim_aligned,
          paste(out, paste0("climdata_era5_cmip6_", cmip6_year, ".csv"), sep = "/"),
          row.names=FALSE)
  cat(paste0("Saved climate data for year ", cmip6_year, " in",
             paste(out, paste0("climdata_era5_cmip6_", cmip6_year, ".csv"), sep = "/"), "\n"))
}