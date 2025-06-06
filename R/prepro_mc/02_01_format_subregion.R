# ------
# Subset the processed data to a 50m x 50m grid and format

library(sf)
library(terra)
library(dplyr)
library(microclimf)
library(readr)


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

# Load data
vegp_baf <- readRDS("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/vegp.RDS")
dtm_baf <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/soil/dtm.tif")
soilc_baf <- readRDS("/Users/johanna/Uni/masterarbeit/data/mc_input/soil/soilc_ptm.RDS")
climdata_baf <- readRDS("/Users/johanna/Uni/masterarbeit/data/mc_input/climate/era5_climdata_2024.RDS")

# -- REGUA

# out <- "/Users/johanna/Uni/masterarbeit/data/mc_input/regua"
# dir.create(file.path(out))
#
# # Define the center point - 22°23′44.75″S, 42°44′15.78″W
# lat <- dms_to_decimal("22°23′44.75″S")
# lon <- dms_to_decimal("42°44′15.78″W")
# center_coords <- st_sfc(st_point(c(lon, lat)), crs = 4326)
# center_proj <- st_transform(center_coords, crs = 31983)
# center_coords_m <- st_coordinates(center_proj)

# -- Pirineus

out <- "/Users/johanna/Uni/masterarbeit/data/mc_input/pirineus"
dir.create(file.path(out))

# Define the center point
lat <- dms_to_decimal("22°26′46.74″S")
lon <- dms_to_decimal("42°30′06.16″W")
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

# - Prep climate data 

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

# -  Prep Vegetation

# Correct layer names 
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
plot(terra::unwrap(vegp_square$hgt))

# - Prep Soil data
soil_square <- terra_crop(soilc_baf, ceiled_box, res = 1)

plot(terra::unwrap(soil_square$soiltype))

# Potentially fill up soil variables if there are NAs
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

# Wrap unwrapped soil variables
soil_square$gref <- terra::wrap(soil_square$gref)
soil_square$slope <- terra::wrap(soil_square$slope)
soil_square$aspect <- terra::wrap(soil_square$aspect)
soil_square$em <- terra::wrap(soil_square$em)

# - Prep DTM
em <- deepcopy(terra::unwrap(soil_square$em))
dtm_baf_res <- resample(dtm_baf, em) 
dtm_square <- dtm_baf_res$lyr1

plot(dtm_square)

# Save data 
saveRDS(vegp_square, paste(out, "vegp_v2.RDS", sep = "/"))
writeRaster(dtm_square, paste(out, "dtm_v2.tif", sep = "/"),
            filetype = "GTiff", overwrite = TRUE)
saveRDS(soil_square, paste(out, "soilc_v2.RDS", sep = "/"))
write.csv(climdata_regua, 
          paste(out, "era5_climdata_2024_v2.csv", sep = "/"),
          row.names=FALSE)

# Visualize data:
vegp_reg <- readRDS(paste(out, "vegp_v2.RDS", sep = "/"))
dtm_reg <- rast(paste(out, "dtm_v2.tif", sep = "/"))
soilc_reg <- readRDS(paste(out, "soilc_v2.RDS", sep = "/"))
climdata_reg <- read_csv(paste(out, "era5_climdata_2024_v2.csv", sep = "/"))

# Vegetation
for (name in names(vegp_reg)) {
  vegp_var <- terra::unwrap(vegp_reg[[name]])
  plot(vegp_var, main = name)
}

print(res(vegp_var))

# DTM
plot(dtm_reg, main = "dtm")

# Soil
for (name in names(soilc_reg)) {
  vegp_var <- terra::unwrap(soilc_reg[[name]])
  plot(vegp_var, main = name)
}
print(res(vegp_var))

# Issue with ground radiation
soilc_reg_unpacked <- lapply(soilc_reg, terra::unwrap)
soilc_reg_unpacked$gref <- terra::ifel(is.na(soilc_reg_unpacked$gref),
                                          unique(soilc_reg_unpacked$gref)[[1]],
                                          soilc_reg_unpacked$gref)
soic_reg <- lapply(soilc_reg_unpacked, terra::wrap)
soic_reg <- soic_reg[names(micropoint::groundparams)]
saveRDS(soic_reg, paste(out, "soilc_v2.RDS", sep = "/"))

