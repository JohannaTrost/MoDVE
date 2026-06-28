library(terra)
library(microclimdata)

# Authenticat
data(credentials) # See

# Function for downloading RGB and CIR
modis_rgb_cir_download <- function(r, tme, pathout, credentials) {
  # Time sanity check
  tsed <- tme[length(tme)]
  tmod <- as.POSIXlt(0, origin = "2000-02-18", tz = "UTC")
  if (tsed < tmod) stop("No MODIS data available prior to 2000-02-18")

  # Prepare AOI from input raster
  e <- terra::ext(r)
  r2 <- terra::rast(e)
  terra::crs(r2) <- terra::crs(r)
  r2 <- terra::project(r2, "EPSG:4326")  # MODIS needs lat/lon
  e2 <- terra::ext(r2)

  # Format time range
  st <- substr(as.character(tme[1]), 1, 10)
  ed <- substr(as.character(tme[length(tme)]), 1, 10)

  # Extract login credentials
  username <- credentials$username[1]
  password <- credentials$password[1]

  # Search for MODIS MOD09GA (surface reflectance) data
  mf <- luna::getNASA("MOD09GA", start = st, end = ed, aoi = e2,
                      version = "061", download = FALSE)

  if (length(mf) > 0) {
    # Download if available
    luna::getNASA("MOD09GA", start = st, end = ed, aoi = e2, version = "061",
                  download = TRUE, path = pathout,
                  username = username, password = password,
                  server = "LPDAAC_ECS")
  } else {
    stop("No MOD09GA data found for the specified area and time period.")
  }
}

normalize <- function(x) {
  x <- terra::clamp(x, lower = 0, upper = 10000)
  return((x / 10000) * 255)
}

# Merge and process RGB and CIR
rgb_cir_process <- function(r, pathin) {
  lst <- list.files(pathin, pattern = "\\.hdf$", full.names = TRUE)
  if (length(lst) == 0) stop("No files to process!")

  # MODIS MOD09GA band indices for Surface Reflectance:
  # Band 1 = Red (index 1)
  # Band 2 = NIR  (index 2)
  # Band 3 = Blue (index 3)
  # Band 4 = Green (index 4)

  pb <- utils::txtProgressBar(min = 0, max = length(lst), style = 3)

  rgb_stack <- list()
  cir_stack <- list()

  for (i in seq_along(lst)) {
    fi <- lst[i]

    # Read all 4 bands (surface reflectance)
    mod_bands <- terra::rast(fi)[[1:4]]  # Bands 1–4

    # AOI transformation to MODIS CRS
    bbx <- terra::rast(terra::ext(r))
    terra::crs(bbx) <- terra::crs(r)
    bbx <- terra::project(bbx, terra::crs(mod_bands))

    # Crop to area of interest
    mod_bands <- terra::crop(mod_bands, bbx)

    # - Build composites

    # Use Band 4 (Green) as reference
    green <- mod_bands[[4]]

    # Resample all bands to match 'green'
    red <- terra::resample(mod_bands[[1]], green, method = "bilinear")
    nir <- terra::resample(mod_bands[[2]], green, method = "bilinear")
    blue <- terra::resample(mod_bands[[3]], green, method = "bilinear")

    rgb <- terra::rast(list(red, green, blue))
    cir <- terra::rast(list(nir, red, green))

    # Stack for temporal processing
    rgb_stack[[i]] <- rgb
    cir_stack[[i]] <- cir

    utils::setTxtProgressBar(pb, i)
  }

  # - Create mean composite over time for each

  # Stack all RGB rasters over time
  rgb_stack_rast <- terra::rast(rgb_stack)

  # Assuming input is 3 bands per date, reshape accordingly
  n_dates <- length(rgb_stack)
  r_means <- terra::app(rgb_stack_rast[[seq(1, n_dates * 3, 3)]], mean, na.rm = TRUE)  # Red
  g_means <- terra::app(rgb_stack_rast[[seq(2, n_dates * 3, 3)]], mean, na.rm = TRUE)  # Green
  b_means <- terra::app(rgb_stack_rast[[seq(3, n_dates * 3, 3)]], mean, na.rm = TRUE)  # Blue

  rgb_mean <- terra::rast(list(r_means, g_means, b_means))

  # Same for CIR
  cir_stack_rast <- terra::rast(cir_stack)
  n_dates <- length(cir_stack)
  n_means <- terra::app(cir_stack_rast[[seq(1, n_dates * 3, 3)]], mean, na.rm = TRUE)  # NIR
  g_means <- terra::app(cir_stack_rast[[seq(2, n_dates * 3, 3)]], mean, na.rm = TRUE)  # Green
  b_means <- terra::app(cir_stack_rast[[seq(3, n_dates * 3, 3)]], mean, na.rm = TRUE)  # Blue

  cir_mean <- terra::rast(list(n_means, g_means, b_means))

  # Reproject to match input raster
  rgb_mean <- terra::project(rgb_mean, terra::crs(r))
  cir_mean <- terra::project(cir_mean, terra::crs(r))

  # Crop to match original extent
  rgb_mean <- terra::crop(rgb_mean, terra::ext(r))
  cir_mean <- terra::crop(cir_mean, terra::ext(r))

  rgb_mean_vis <- terra::rast(lapply(terra::as.list(rgb_mean), normalize))

  # Combine into a 6-band raster: RGB (1–3), CIR (4–6)
  out <- c(rgb_mean_vis, cir_mean)

  # Assign names only if layer count matches expected
  if (terra::nlyr(out) == 6) {
    names(out) <- c("RGB_R", "RGB_G", "RGB_B", "CIR_NIR", "CIR_R", "CIR_G")
  } else {
    warning(paste("Expected 6 bands, found", terra::nlyr(out)))
  }

  return(out)
}


# Directories
figs_dir <- file.path("..", "modve_figs", "mc_input")
data_dir <- file.path("..", "modve_data")
out_dir <- file.path(data_dir, "mc_input", "albedo")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Define extent in decimal degrees
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
r <- rast(e)
crs(r) <- "EPSG:4326"

# Reproject to meters 
r <- project(r, "EPSG:31983", res = 500)

# --- Get albedo

tme <- as.POSIXlt(c(0:30) * 3600 * 24, origin = "2024-02-01 00:00", tz = "UTC")
albedo_download(r, tme, out_dir, credentials)

# Questionable albedo
albmodis <- albedo_process(r, out_dir)

# Save
writeRaster(albmodis, file.path(out_dir, "albedo_modis_2024-01-01.tif"), filetype = "GTiff", overwrite = TRUE)

# --- Get RGB and CIR

pathout <- file.path(out_dir, "rgb_cir_2024-02-01")

# Download data for RGB and CIR 
modis_rgb_cir_download(r, tme, pathout, credentials)

# Process RGB and CIR 
modis_rgb <- rgb_cir_process(r = r, pathin = pathout)

# Save RGB and CIR
if(!dir.exists(file.path(out_dir, "rgb_cir_2024-02-01"))) {
  dir.create(file.path(out_dir, "rgb_cir_2024-02-01"))
}
writeRaster(modis_rgb[[1:3]], file.path(out_dir, "rgb_cir_2024-02-01", "rgb.tif"), filetype = "GTiff", overwrite = TRUE)
writeRaster(modis_rgb[[4:6]], file.path(out_dir, "rgb_cir_2024-02-01", "cir.tif"), filetype = "GTiff", overwrite = TRUE)

# -- Plot the resulting albedo, RGB and CIR data

albphoto <- albedo_fromaerial(modis_rgb[[1:3]], modis_rgb[[4:6]])

# Plot 1: Albedo from aerial
png(file.path(figs_dir, "01_albedo_from_aerial.png"))
plot(albphoto, col=gray.colors(255), range = c(0, 1))
dev.off()

# Plot 2: Albedo from aerial
albadjusted <- albedo_adjust(albphoto, albmodis)
png(file.path(figs_dir, "02_albedo_from_aerial_adjusted.png"))
plot(albadjusted, col=gray.colors(255), range = c(0, 1))
dev.off()

# Plot 3: Albedo from modis
png(file.path(figs_dir, "03_albedo.png"))
plot(albmodis, col = gray.colors(255), range = c(0, 1))
dev.off()

# Plot 4: RGB
png(file.path(figs_dir, "04_rgb.png"))
plotRGB(modis_rgb[[1:3]], stretch = "lin")
dev.off()

# Plot 5: CIR
png(file.path(figs_dir, "04_cir.png"))
plotRGB(modis_rgb[[4:6]], stretch = "lin")
dev.off()



