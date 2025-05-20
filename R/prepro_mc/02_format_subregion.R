# ------
# Subset the processed data to a 50m x 50m grid and format

library(sf)
library(terra)
library(dplyr)
library(microclimf)


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



vegp_baf <- readRDS("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/vegp.RDS")
dtm_baf <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/soil/dtm.tif")
soilc_baf <- readRDS("/Users/johanna/Uni/masterarbeit/code/data/mc_input/soil/soilc.RDS")
climdata_baf <- readRDS("/Users/johanna/Uni/masterarbeit/code/data/mc_input/climate/era5_climdata_2016.RDS")

# -- REGUA

out <- "/Users/johanna/Uni/masterarbeit/code/data/mc_input/regua"
dir.create(file.path(out))

# Define the center point
lat <- -22.426880
lon <- -42.765096
center_coords <- st_sfc(st_point(c(lon, lat)), crs = 4326)
center_proj <- st_transform(center_coords, crs = 31983)
center_coords_m <- st_coordinates(center_proj)

# Create a 50m x 50m square (i.e., a 50m buffer in all directions)
extent_box <- st_buffer(center_proj, dist = 25, endCapStyle = "SQUARE")
extent_vect <- vect(extent_box)

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
vegp_square <- terra_crop(unwrapped_list, extent_vect, res = 1)
plot(terra::unwrap(vegp_square$hgt))

# - Prep DTM
pai <- deepcopy(terra::unwrap(vegp_square$pai))
dtm_baf_res <- resample(dtm_baf, pai) 
dtm_square <- crop(dtm_baf_res$lyr1, extent_vect)

plot(dtm_square)

# - Prep Soil data
soil_square <- terra_crop(soilc_baf, extent_vect, res = 1)



plot(terra::unwrap(soil_square$soiltype))

# Save data 
saveRDS(vegp_square, paste(out, "vegp.RDS", sep = "/"))
writeRaster(dtm_square, paste(out, "dtm.tif", sep = "/"), 
            filetype = "GTiff", overwrite = TRUE)
saveRDS(soil_square, paste(out, "soilc.RDS", sep = "/"))
write.csv(climdata_regua, 
          paste(out, "era5_climdata_2016.csv", sep = "/"), 
          row.names=FALSE)

# Visualize data:
in_dir <- "/Users/johanna/Uni/masterarbeit/code/data/mc_input/regua"
vegp_reg <- readRDS(paste(in_dir, "vegp.RDS", sep = "/"))
dtm_reg <- rast(paste(in_dir, "dtm.tif", sep = "/"))
soilc_reg <- readRDS(paste(in_dir, "soilc.RDS", sep = "/"))
climdata_reg <- read_csv(paste(in_dir, "era5_climdata_2024.csv", sep = "/"))

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
soilc_reg_unpacked$groundr <- terra::ifel(is.na(soilc_reg_unpacked$groundr), 
                                          0.009858172, soilc_reg_unpacked$groundr)
soic_reg <- lapply(soilc_reg_unpacked, terra::wrap)
saveRDS(soic_reg, paste(out, "soilc.RDS", sep = "/"))

