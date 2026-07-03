# Generate soil parameters and ground reflectance

library(terra)
library(httr)
library(microclimdata)
library(micropoint)

# Get soil parameters for the grid from soil type mapping 
get_by_soiltype <- function(soil_raster) {
  # Soil parameters table (already provided)
  soilparams <- micropoint::soilparams
  
  # Exclude specific columns
  exclude_vars <- "Soil.type"
  target_vars <- setdiff(colnames(soilparams), exclude_vars)
  
  # Soil type codes in raster correspond to rows in soilparams:
  soilparams$soil_code <- seq_len(nrow(soilparams))
  
  # Create a raster for each variable in soilparams
  raster_list <- list()
  for (var in target_vars) {

    # Reclassify the raster
    new_raster <- classify(soil_raster, rcl = cbind(soilparams$soil_code, soilparams[[var]]))
    names(new_raster) <- var # Set layer name
    
    # Add to list
    raster_list[[var]] <- new_raster
  }
  
  return(raster_list)
}

# Directories
figs_dir <- file.path("..", "modve_figs", "mc_input")
data_dir <- file.path("..", "modve_data")
soil_dir <- file.path(data_dir, "mc_input", "soil")
veg_dir <- file.path(data_dir, "mc_input", "vegetation")
albedo_dir <- file.path(data_dir, "mc_input", "albedo")

# Check and create directories
if(!dir.exists(soil_dir)) {
  dir.create(soil_dir)
}
if(!dir.exists(figs_dir)) {
  dir.create(figs_dir)
}
if(!dir.exists(veg_dir)) {
  stop("Run 01_vegetation.R to obtain required vegetation data.")
}
if(!dir.exists(albedo_dir)) {
  stop("Run 02_albedo.R to obtain required albedo data.")
}

# Load data
albedo <- rast(file.path(albedo_dir, "albedo_modis_2024-01-01.tif"))
lai <- rast(file.path(veg_dir, "lai_2024_500m.tif"))
landcover <- rast(file.path(veg_dir, "ESA_WorldCover_2023_10m.tif"))

# Create a template raster with 100 m resolution
e <- ext(-42.9, -42.2, -22.8, -22.2)  # xmin, xmax, ymin, ymax
template <- rast(e)
crs(template) <- "EPSG:4326"
template <- project(template, "EPSG:31983", res = 100)

# Resample each raster to 100 m resolution using the template
albedo_100m <- resample(albedo, template, method="bilinear")   
lai_100m <- resample(lai, template, method="bilinear")         
landcover_100m <- resample(landcover, template, method="near")

# Derive leaf and ground reflectance
x <- x_calc(landcover_100m, lctype = "ESA")
x <- mask(x, lai_100m)
refldata <- reflectance_calc(albedo_100m, lai_100m, x, plotprogress=FALSE)

# Save reflectance data
writeRaster(refldata$gref, file.path(soil_dir, "ground_ref_100m.tif"), filetype = "GTiff", overwrite = TRUE)
writeRaster(refldata$lref, file.path(soil_dir, "leaf_ref_100m.tif"), filetype = "GTiff", overwrite = TRUE)

# Download soil data using microclimdata package
r <- rast(file.path(veg_dir, "rgb_cir_2024-02-01", "cir.tif"))
soildata <- soildata_download(r[[1]], pathdir = file.path(data_dir, "tmp"), deletefiles = TRUE)

# --- Generate input data for micropoint model

# Generate data for microclimf (grid model)
soilc <-  create_soilgrid(soildata, refldata, landcover, water = 80)
attributes(soilc)

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

# Get DTM
dtm <- dem_download(r) # from CIR

# Get additional data for point model
ptm_soilc$slope <- terra::terrain(dtm$lyr1,'slope')
ptm_soilc$aspect <- terra::terrain(dtm$lyr1,'aspect')

# Add emissivity (0.97 is taken from micropoint src code)
em <- 0.97
ptm_soilc$em <- deepcopy(ptm_soilc$slope)
values(ptm_soilc$em) <- em

# Use correct order of variables for the model
var_order <- append(names(micropoint::groundparams), "soiltype")
ptm_soilc <- ptm_soilc[var_order]

# Wrap all rasters
ptm_soilc_wrp <- lapply(ptm_soilc, terra::wrap)

# Save soil data
saveRDS(ptm_soilc_wrp, file.path(soil_dir, "soil.RDS"))

# --- Visualize data

png(file.path(figs_dir, "soiltype.png"))
plot(terra::unwrap(ptm_soilc$soiltype), main = "Soil Types")
dev.off()

png(file.path(figs_dir, "dtm.png"))
plot(dtm[[1]])
dev.off()

png(file.path(figs_dir, "groundr.png"))
plot(terra::unwrap(ptm_soilc$gref), col=gray.colors(255), range = c(0, 1), main = "Ground Reflectance")
dev.off()

png(file.path(figs_dir, "rho.png"))
plot(terra::unwrap(ptm_soilc$rho), main = expression("Bulk Density (Mg/m"^3*")"))
dev.off()

png(file.path(figs_dir, "Vm.png"))
plot(terra::unwrap(ptm_soilc$Vm), main = expression("Volumetric Mineral Fraction (m"^3*"/m"^3*")"))
dev.off()

png(file.path(figs_dir, "Vq.png"))
plot(terra::unwrap(ptm_soilc$Vq), main = expression("Volumetric Quartz Fraction (m"^3*"/m"^3*")"))
dev.off()

png(file.path(figs_dir, "Mc.png"))
plot(terra::unwrap(ptm_soilc$Mc), main = expression("Mass Fraction of Clay (kg/kg)"))
dev.off()