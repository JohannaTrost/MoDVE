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


# CONFIGURE Region
REGION <- "regua"
# Center point of forest stand
lat <- dms_to_decimal("22°23′44.75″S")
lon <- dms_to_decimal("42°44′15.78″W")

# For Pirineus uncomment below
# REGION <- "pirineus"
# Center point of forest stand
# lat <- dms_to_decimal("22°26′46.74″S")
# lon <- dms_to_decimal("42°30′06.16″W")

# Function to crop Spatratser object to given region
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

# Directories
figs_dir <- file.path("..", "modve_figs", "mc_input")
data_dir <- file.path("..", "modve_data")
veg_dir <- file.path(data_dir, "mc_input", "vegetation")
soil_dir <- file.path(data_dir, "mc_input", "soil")
out <- file.path(data_dir, REGION)

if(!dir.exists(out)) {
  dir.create(file.path(out))
}

# Load data
vegp_baf <- readRDS(veg_dir, "vegp.RDS")
soilc_baf <- readRDS(soil_dir, "soil.RDS")

# Project coordinates
center_coords <- st_sfc(st_point(c(lon, lat)), crs = 4326)
center_proj <- st_transform(center_coords, crs = 31983)
center_coords_m <- st_coordinates(center_proj)

# Create a 50m x 50m square (i.e., a 50m buffer in all directions)
extent_box <- st_buffer(center_proj, dist = 25, endCapStyle = "SQUARE")
extent_vect <- vect(extent_box)

# Ceil the extent if decimals
e <- ext(extent_vect)
ce <- ext(round(e[1]), round(e[2]), round(e[3]), round(e[4]))

# Create new polygon from ceiled extent
ceiled_box <- as.polygons(ce)
crs(ceiled_box) <- crs(extent_vect)  # Set CRS to match original

# --------------------------- Prep Vegetation --------------------------- #

# Correct layer names (turns multi-layer SpatRaster into a named list of single-layer SpatRaster objects)
layer_names <- names(vegp_baf)
unwrapped_list <- setNames(
  lapply(layer_names, function(nm) {
    r <- terra::unwrap(vegp_baf[[nm]])
    names(r) <- nm
    terra::wrap(r)
  }),
  layer_names
)
vegp_square <- terra_crop(unwrapped_list, ceiled_box, res = 1)

# --------------------------- Prep Soil data --------------------------- #

soil_square <- terra_crop(soilc_baf, ceiled_box, res = 1)

# Potentially fill soil variables if there are NAs
for (name in names(soil_square)) {

  raster <- terra::unwrap(soil_square[[name]])
  
  # Check if there are any NAs
  if (anyNA(values(raster))) {
    # Get unique non-NA values
    unique_vals <- unique(na.omit(values(raster)))
    
    if (length(unique_vals) == 1) {
      # Fill NAs with the single unique value
      values(raster)[is.na(values(raster))] <- unique_vals
      soil_square[[name]] <- raster  # Replace with modified raster
      cat(paste0("Filled NAs in '", name, "' with value: ", unique_vals, "\n"))
    } else {
      cat(paste0("'", name, "': Raster contains NAs and multiple values.\n"))
    }
  }
}

# Fix NAs in ground radiation
soil_square$gref <- terra::ifel(is.na(soil_square$gref), unique(soil_square$gref)[[1]], soil_square$gref)

# Wrap unwrapped soil variables
soil_square$gref <- terra::wrap(soil_square$gref)
soil_square$slope <- terra::wrap(soil_square$slope)
soil_square$aspect <- terra::wrap(soil_square$aspect)
soil_square$em <- terra::wrap(soil_square$em)

# Correct names for Micropoint
soil_square <- soil_square[names(micropoint::groundparams)]

# Save data 
saveRDS(vegp_square, file.path(out, "vegp.RDS"))
saveRDS(soil_square, file.path(out, "soil.RDS"))

# -- Visually check data

# Vegetation
for (name in names(vegp_square)) {
  vegp_var <- terra::unwrap(vegp_square[[name]])
  plot(vegp_var, main = name)
}
print(res(vegp_var))

# Soil
for (name in names(soil_square)) {
  soil_var <- terra::unwrap(soil_square[[name]])
  plot(soil_var, main = name)
}
print(res(soil_var))

