library(terra)
library(microclimdata)
library(httr)

albedo <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/albedo/albedo_modis_2024-01-01.tif")
lai <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/lai_2024_500m.tif")
landcover <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/ESA_WorldCover_2023_10m.tif")

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

writeRaster(refldata$gref, "/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/ground_ref_100m.tif", 
            filetype = "GTiff", overwrite = TRUE)
writeRaster(refldata$lref, "/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/leaf_ref_100m.tif", 
            filetype = "GTiff", overwrite = TRUE)

r <- rast("/Users/johanna/Uni/masterarbeit/code/data/mc_input/vegetation/rgb_cir_2024-02-01/cir.tif")[[1]]
soildata <- soildata_download(r, 
                              pathdir = "/Users/johanna/Uni/masterarbeit/code/data/mc_input/tmp", 
                              deletefiles = TRUE)
soilc <-  create_soilgrid(soildata, refldata, landcover, water = 80)
attributes(soilc)

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


