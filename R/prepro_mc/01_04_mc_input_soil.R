library(terra)
library(httr)

# library(remotes)
# remotes::install_local("/Users/johanna/Uni/masterarbeit/code/microclimdata")

library(microclimdata)
library(micropoint)

# Get soil parameters for the grid from soil type mapping 
get_by_soiltype <- function(soil_raster) {
  # Soil parameters table (already provided)
  soilparams <- micropoint::soilparams
  
  # Exclude specific columns
  exclude_vars <- c("Soil.type")
  target_vars <- setdiff(colnames(soilparams), exclude_vars)
  
  # soil type codes in raster correspond to rows in soilparams:
  soilparams$soil_code <- 1:nrow(soilparams)
  
  # Create a raster for each variable in soilparams
  raster_list <- list()
  for (var in target_vars) {
    # Create a vector to map soil type code to variable value
    map_values <- setNames(soilparams[[var]], soilparams$soil_code)
    
    # Reclassify the raster
    new_raster <- classify(soil_raster, 
                           rcl = cbind(soilparams$soil_code, soilparams[[var]]))
    names(new_raster) <- var # Set layer name
    
    # Add to list
    raster_list[[var]] <- new_raster
  }
  
  return(raster_list)
}


albedo <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/albedo/albedo_modis_2024-01-01.tif")
lai <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/lai_2024_500m.tif")
landcover <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/ESA_WorldCover_2023_10m.tif")

# Create a template raster with 100 m resolution
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
template <- rast(e)
crs(template) <- "EPSG:4326"
template <- project(template, "EPSG:31983", res = 100)

# Resample each raster to 100 m resolution using the template
albedo_100m <- resample(albedo, template, method="bilinear")   
lai_100m <- resample(lai, template, method="bilinear")         
landcover_100m <- resample(landcover, template, method="near")

x <- x_calc(landcover_100m, lctype = "ESA")
x <- mask(x, lai_100m)

# Derive leaf and ground reflectance
refldata <- reflectance_calc(albedo_100m, lai_100m, x, plotprogress=FALSE)

# Plot outputs
par(mfrow = c(1, 2))
plot(refldata$gref, col=gray.colors(255), range = c(0, 1), main = "Ground")
plot(refldata$lref, col=gray.colors(255), range = c(0, 1), main = "Leaf")

writeRaster(refldata$gref, "/Users/johanna/Uni/masterarbeit/data/mc_input/soil/ground_ref_100m.tif", 
            filetype = "GTiff", overwrite = TRUE)
writeRaster(refldata$lref, "/Users/johanna/Uni/masterarbeit/data/mc_input/soil/leaf_ref_100m.tif", 
            filetype = "GTiff", overwrite = TRUE)

r <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/rgb_cir_2024-02-01/cir.tif")[[1]]
soildata <- soildata_download(r, 
                              pathdir = "/Users/johanna/Uni/masterarbeit/data/mc_input/tmp/", 
                              deletefiles = TRUE)


# --- Generate input data for micropoint model

# TODO get the data at this point

# --- Generate data for microclimf (grid model)

soilc <-  create_soilgrid(soildata, refldata, landcover, water = 80)
attributes(soilc)

# Loading data instead of computing the above / comment if necessary
soilc <- readRDS("/Users/johanna/Uni/masterarbeit/data/mc_input/soil/soilc.RDS")

# Adapt name for point model
names(soilc)[names(soilc) == "groundr"] <- "gref"

# Extract missing variables
ptm_soilc <- lapply(soilc, terra::unwrap)
soilp <- get_by_soiltype(ptm_soilc$soiltype)

# Add missing variables: Psie, Smax, Smin etc.
for (var_name in names(soilp)) {
  if (var_name == "psi_e") {
    ptm_soilc[[var_name]] <- -soilp[[var_name]]
  } else {
    ptm_soilc[[var_name]] <- soilp[[var_name]]
  }
}
# Rename for correct structure for pointmodel
names(ptm_soilc)[names(ptm_soilc) == "psi_e"] <- "Psie"

# Get additional data for point model
dtm <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/soil/dtm.tif")
ptm_soilc$slope <- terra::terrain(dtm$lyr1,'slope')
ptm_soilc$aspect <- terra::terrain(dtm$lyr1,'aspect')

# Add emissivity (0.97 is taken from micropoint src code)
em <- 0.97
ptm_soilc$em <- deepcopy(ptm_soilc$slope)
values(ptm_soilc$em) <- em

# Use correct order of variables for the model
ptm_soilc <- ptm_soilc[names(micropoint::groundparams)]

# Wrap all rasters
ptm_soilc_wrp <- lapply(ptm_soilc, terra::wrap)

saveRDS(ptm_soilc_wrp, 
        "/Users/johanna/Uni/masterarbeit/data/mc_input/soil/soilc_ptm.RDS")

# SAve plots
out_dir <- "/Users/johanna/Uni/masterarbeit/figs/mc_input"

png(file.path(out_dir, "soiltype.png"))
plot(terra::unwrap(soilc$soiltype), main = "Soil Types")
dev.off()

png(file.path(out_dir, "groundr.png"))
plot(terra::unwrap(soilc$groundr), main = "Ground Reflectance")
dev.off()

png(file.path(out_dir, "rho.png"))
plot(terra::unwrap(soilc$rho), main = expression("Bulk Density (Mg/m"^3*")"))
dev.off()

png(file.path(out_dir, "Vm.png"))
plot(terra::unwrap(soilc$Vm), main = expression("Volumetric Mineral Fraction (m"^3*"/m"^3*")"))
dev.off()

png(file.path(out_dir, "Vq.png"))
plot(terra::unwrap(soilc$Vq), main = expression("Volumetric Quartz Fraction (m"^3*"/m"^3*")"))
dev.off()

png(file.path(out_dir, "Mc.png"))
plot(terra::unwrap(soilc$Mc), main = expression("Mass Fraction of Clay (kg/kg)"))
dev.off()
