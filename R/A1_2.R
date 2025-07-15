# Add Microclimate to Microhabitat Matrix
source("utils.R")

library("fs")

Mc2MhMatrix <- function(MicrohabDims, MicroclimDims, Microhabitat, Microclimate, indices) {
  newD4 <- MicrohabDims[4] + MicroclimDims[4]  # new size of 4th dimension

  # Initialize the new array
  NewMicrohabitat <- array(NA, dim = c(MicrohabDims[1:3], newD4))
  # Copy existing microhabitat data
  McIdxNames <- c("HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
  for (idx in indices[!(names(indices) %in% McIdxNames)]) {
      NewMicrohabitat[,,,idx] <- Microhabitat[,,,idx]
  }
  # Insert microclimate data into the new array
  McIndices <- indices[names(indices) %in% McIdxNames]
  for (i in seq_along(McIndices)) {  # We use seq_along as Hum/Temp/Wind order won't change
    idx <- McIndices[[i]]
    NewMicrohabitat[,,,idx] <- Microclimate[,,,i]
  }

    return(NewMicrohabitat)
}

config <- parse_config()

# start and end timestep
timeStepStart <- config$timeStepStart
timeStepEnd <- config$timeStepEnd

# Directory paths
DirectoryMicrohabitat <- config$DirectoryMicrohabitat
DirectoryMicroclimate <- config$DirectoryMicroclimate
DirectoryNewMicrohabitat <- config$DirectoryNewMicrohabitat

# Creat directory for new microhabitat if it doesn't exist
if (!dir.exists(DirectoryNewMicrohabitat)) {
  dir.create(DirectoryNewMicrohabitat, recursive=TRUE)
}

# Define option flags as a named list
options_list <- list(
  TotalSurfaceAreaOpt = config$TotalSurfaceAreaOpt,
  SurfaceAreaLossOpt  = config$SurfaceAreaLossOpt,
  LightNicheOpt  = config$LightNicheOpt,
  AverageWeightedAngles = config$AverageWeightedAngles,
  HumNicheOpt    = config$HumNicheOpt,
  TempNicheOpt   = config$TempNicheOpt,
  WindNicheOpt   = config$WindNicheOpt
)

# Only keep active options
active_options <- names(options_list[options_list == 1])

# Assign indices dynamically
indices <- setNames(seq_along(active_options), active_options)

# Extract indieces
HumIndex  <- indices["HumNicheOpt"]
TempIndex <- indices["TempNicheOpt"]
WindIndex <- indices["WindNicheOpt"]

for (TimeStep in int_seq(from=timeStepStart, to=timeStepEnd, by=1)) {

  print(paste("Time step", TimeStep))

  # Load  microhabitat matrix
  MhFile <- paste("MicrohabitatMatrix", TimeStep, ".rds", sep="")
  FileMatrix <- file.path(DirectoryMicrohabitat, MhFile)
  Microhabitat <- readRDS(FileMatrix)

  # Load microclimate matrix
  McFile <- paste("MicroclimateMatrix", "1", ".rds", sep="")
  FileMcMatrix <- file.path(DirectoryMicroclimate, McFile)
  Microclimate <- readRDS(FileMcMatrix)

  # 1. Select relevant microclimate variables
  Microclimate <- Microclimate[,,,c(7, 1, 11)] # mean annual relhum, tair, windspeed

  # 2. Fill up above canopy microclimate given microhabitat height
  MicrohabDims <- dim(Microhabitat)
  MicroclimDims <- dim(Microclimate)
  HeightDiff <- MicrohabDims[3] - MicroclimDims[3]

  if (HeightDiff > 0) { # if microhabitat is taller than microclimate vertical range
    # Fill missing above canopy climate
    Microclimate <- abind::abind( # Combine Mc with repeated last layer
      Microclimate,
      array(  # Repeat last vertical layer as many times as required for Microhabitat height
        rep(Microclimate[,,MicroclimDims[3],], each = HeightDiff),
        dim = c(MicrohabDims[1:2], HeightDiff, MicroclimDims[4])
      ),
      along = 3
    )

  } else if (HeightDiff < 0) {  # Typically this should not happen, but if it does:
    # Trim microclimate to match microhabitat height
    Microclimate <- Microclimate[,,1:MicrohabDims[3],]
  }

  # Ensure dimensions match
  if (!all(dim(Microhabitat)[1:3] == dim(Microclimate)[1:3])) {
    stop("Microhabitat and Microclimate dimensions do not match after height adjustment.")
  }

  # 3. Insert microclimate into microhabitat matrix
  NewMhMatrix <- Mc2MhMatrix(
    MicrohabDims = MicrohabDims,
    MicroclimDims = MicroclimDims,
    Microhabitat = Microhabitat,
    Microclimate = Microclimate,
    indices = indices
  )

  # 4. Update microhabitat matrix
  NewMhFileMatrix <- file.path(DirectoryNewMicrohabitat, MhFile)
  saveRDS(NewMhMatrix, file = NewMhFileMatrix)
}

# Copy necessary forest files
file_copy(
  file.path(DirectoryMicrohabitat, "Forest_param_global.txt"),
  file.path(DirectoryNewMicrohabitat, "Forest_param_global.txt")
)
file_copy(
  file.path(DirectoryMicrohabitat, "Forest_param_pass0.txt"),
  file.path(DirectoryNewMicrohabitat, "Forest_param_pass0.txt")
)
file_copy(
  file.path(DirectoryMicrohabitat, "dimPlot.rds"),
  file.path(DirectoryNewMicrohabitat, "dimPlot.rds")
)