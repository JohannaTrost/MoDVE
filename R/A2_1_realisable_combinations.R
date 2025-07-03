options(warn=-1)  # Suppress warnings
options(digits.secs=3)  # 3 decimal digits for seconds

# Epiphte IBM - Model
# This model simulates the development of the entire epiphyte community
source("utils.R")

library("readr")

main <- function() {
  # Load the configuration file
  config <- parse_config()

  DirectoryOutput <- config$DirectoryOutput
  DirectoryMicrohabitat <- config$DirectoryMicrohabitat
  initialTimeStep <- config$initialTimeStep
  timeSteps <- config$timeSteps
  chunkSize <- config$chunkSize

  # Dynamic handling of microhabitat dimensions
  MicrohabitatVariableFlags <- config$MicrohabitatVariableFlags
  MhVarNames <- c("TotalSurfaceAreaOpt", "SurfaceAreaLossOpt", "LightNicheOpt", "AverageWeightedAngles",
                      "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
  # Only keep active options
  ActiveOpts <- MhVarNames[as.logical(MicrohabitatVariableFlags)]
  Inds <- setNames(seq_along(ActiveOpts), ActiveOpts)
  EnvVarNames <- c("LightNicheOpt", "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
  EnvInds <- Inds[EnvVarNames]

  # --- Load all Microhabitat matrices for the initial time step

  # Load plot dimensions if an artifical theoretical forest is used
  dimPlot <- readRDS(file.path(DirectoryMicrohabitat, "dimPlot.rds"))

  # Initialize empty dataframe with Light Hum Temp and Wind columns
  UniqueEnvVarsComb <- data.frame(
    Light = numeric(0),
    Hum = numeric(0),
    Temp = numeric(0),
    Wind = numeric(0)
  )

  # Looping over chunks of time steps
  firstTimeStep <- initialTimeStep
  lastTimeStep <- initialTimeStep + timeSteps
  starts <- seq(firstTimeStep, lastTimeStep - chunkSize, by = chunkSize)
  ends <- c((starts + chunkSize - 1)[seq_along(starts) -1], lastTimeStep)

  for (i in seq_along(starts)) {

    # Init microhabitat matrix for the current chunk
    McMatrix <- array(NA, dim = c(ends[i] - starts[i] + 1, dimPlot, length(EnvInds)))
    TimeSteps <- seq(starts[i], ends[i])

    # Load microhabitat matrices for the current chunk
    for (j in seq_along(TimeSteps)) {
      timeStep <- TimeSteps[j]

      print(paste0("Time step: ", timeStep))

      FileMhMatrix <- file.path(DirectoryMicrohabitat, paste0("MicrohabitatMatrix", timeStep, ".rds"))
      if (!file.exists(FileMhMatrix)) {
          stop(paste("Microhabitat matrix file does not exist:", FileMhMatrix))
      } else {
          McMatrix[j,,,,] <- readRDS(FileMhMatrix)[,,, EnvInds]
      }
    }

    # -- Extract unique combinations of environmental niches from the microclimate matrix

    # 1. Round variables, humidity/temp/wind to 2 decimals, light to 5 decimals
    McMatrix[,,,,2:4] <- round(McMatrix[,,,,2:4])  # Humidity, Temp, Wind
    McMatrix[,,,,1] <- round(McMatrix[,,,,1], 3)  # Light

    # 2. Extract unique combinations of environmental niches
    UniqueEnvVarsCombChunk <- McMatrix
    dim(UniqueEnvVarsCombChunk) <- c(prod(dim(McMatrix)[1:4]), 4)
    # Convert to dataframe and assign column names
    UniqueEnvVarsCombChunkDf <- as.data.frame(UniqueEnvVarsCombChunk)
    colnames(UniqueEnvVarsCombChunkDf) <- c("Light", "Hum", "Temp", "Wind")

    # Remove rows with any NAs
    UniqueEnvVarsCombChunkDf <- na.omit(UniqueEnvVarsCombChunkDf)

    # Remove duplicate rows
    UniqueEnvVarsCombChunkDf <- unique(UniqueEnvVarsCombChunkDf)

    # 5. Add to the main dataframe
    UniqueEnvVarsComb <- rbind(UniqueEnvVarsComb, UniqueEnvVarsCombChunkDf)
  }

  # Again remove rublicates from all chunks
  UniqueEnvVarsComb <- unique(UniqueEnvVarsComb)

  write_csv(
    UniqueEnvVarsComb,
    file.path(DirectoryOutput, paste0("UniqueEnvVarsComb_t", firstTimeStep, "-t", lastTimeStep, ".csv"))
  )
}

main()




