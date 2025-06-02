library(microclimdata)

r <- rast("/Users/johanna/Uni/masterarbeit/data/mc_input/vegetation/rgb_cir_2024-02-01/cir.tif")
dtm <- dem_download(r)

plot(dtm[[1]])

writeRaster(dtm, "/Users/johanna/Uni/masterarbeit/data/mc_input/soil/dtm.tif", 
            filetype = "GTiff", overwrite = TRUE)
