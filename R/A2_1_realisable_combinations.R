options(warn = -1)          # Suppress warnings
options(digits.secs = 3)    # 3 decimal digits for seconds
setwd("/home/jtrost_ext/MoDVE/R/")

# Epiphyte IBM - Model
source("utils.R")

library(readr)
library(data.table)

main <- function() {
  # Load the configuration file
  config <- parse_config()

  DirectoryOutput <- config$DirectoryOutput
  DirectoryMicrohabitat <- config$DirectoryMicrohabitat
  initialTimeStep <- config$initialTimeStep
  timeSteps <- config$timeSteps
  chunkSize <- config$chunkSize

  if (!dir.exists(DirectoryOutput)) {
    dir.create(DirectoryOutput, recursive = TRUE)
  }

  # Dynamic handling of microhabitat dimensions
  MicrohabitatVariableFlags <- config$MicrohabitatVariableFlags
  MhVarNames <- c(
    "TotalSurfaceAreaOpt", "SurfaceAreaLossOpt", "LightNicheOpt",
    "AverageWeightedAngles", "HumNicheOpt", "TempNicheOpt", "WindNicheOpt"
  )

  ActiveOpts <- MhVarNames[as.logical(MicrohabitatVariableFlags)]
  Inds <- setNames(seq_along(ActiveOpts), ActiveOpts)

  default_env_order <- c("LightNicheOpt", "HumNicheOpt", "TempNicheOpt", "WindNicheOpt")
  env_vars_present <- default_env_order[default_env_order %in% names(Inds)]

  if (length(env_vars_present) == 0) {
    stop("No environmental niche axes are active. At least one of Light/Hum/Temp/Wind must be enabled.")
  }

  EnvInds <- Inds[env_vars_present]

  name_map <- c(
    LightNicheOpt = "Light",
    HumNicheOpt   = "Hum",
    TempNicheOpt  = "Temp",
    WindNicheOpt  = "Wind"
  )
  env_cols <- unname(name_map[env_vars_present])

  digits_map <- c(
    LightNicheOpt = 3L,
    HumNicheOpt   = 0L,
    TempNicheOpt  = 0L,
    WindNicheOpt  = 0L
  )

  # Load plot dimensions
  dimPlot <- readRDS(file.path(DirectoryMicrohabitat, "dimPlot.rds"))

  UniqueEnvVarsCombDT <- as.data.table(
    setNames(
      replicate(length(env_cols), numeric(0), simplify = FALSE),
      env_cols
    )
  )
  UniqueEnvVarsCombDT[, Count := integer()]

  # --- Track raw limits before rounding
  limits <- data.table(
    Variable = env_cols,
    Min = rep(Inf, length(env_cols)),
    Max = rep(-Inf, length(env_cols))
  )

  # Looping over chunks
  firstTimeStep <- initialTimeStep
  lastTimeStep  <- initialTimeStep + timeSteps

  starts <- seq(firstTimeStep, lastTimeStep, by = chunkSize)
  ends   <- pmin(starts + chunkSize - 1L, lastTimeStep)

  for (i in seq_along(starts)) {
    McMatrix <- array(NA_real_, dim = c(ends[i] - starts[i] + 1L, dimPlot, length(EnvInds)))
    TimeSteps <- seq(starts[i], ends[i])

    for (j in seq_along(TimeSteps)) {
      timeStep <- TimeSteps[j]
      message(paste0("Time step: ", timeStep))

      FileMhMatrix <- file.path(DirectoryMicrohabitat, paste0("MicrohabitatMatrix", timeStep, ".rds"))
      if (!file.exists(FileMhMatrix)) {
        stop(paste("Microhabitat matrix file does not exist:", FileMhMatrix))
      } else {
        McMatrix[j, , , , ] <- readRDS(FileMhMatrix)[, , , EnvInds, drop = FALSE]
      }
    }

    # ---- Update raw limits BEFORE rounding
    for (k in seq_along(env_vars_present)) {
      raw_vals <- McMatrix[, , , , k]
      limits[k, Min := min(Min, min(raw_vals, na.rm = TRUE))]
      limits[k, Max := max(Max, max(raw_vals, na.rm = TRUE))]
    }

    # ---- Apply rounding rules
    for (k in seq_along(env_vars_present)) {
      var_name <- env_vars_present[k]
      digits   <- digits_map[[var_name]]
      McMatrix[, , , , k] <- round(McMatrix[, , , , k], digits = digits)
    }

    # Reshape and convert to data.table
    dim(McMatrix) <- c(prod(dim(McMatrix)[1:4]), length(env_vars_present))
    ChunkDT <- as.data.table(McMatrix)
    setnames(ChunkDT, env_cols)
    ChunkDT <- na.omit(ChunkDT)

    # Count unique combos
    ChunkCounts <- ChunkDT[, .N, by = env_cols]
    setnames(ChunkCounts, "N", "Count")

    if (nrow(UniqueEnvVarsCombDT) == 0L) {
      UniqueEnvVarsCombDT <- copy(ChunkCounts)
    } else {
      UniqueEnvVarsCombDT <- merge(
        UniqueEnvVarsCombDT, ChunkCounts,
        by = env_cols,
        all = TRUE,
        suffixes = c(".old", ".new")
      )
      UniqueEnvVarsCombDT[, Count :=
        fifelse(is.na(Count.old), 0L, Count.old) +
        fifelse(is.na(Count.new), 0L, Count.new)
      ]
      UniqueEnvVarsCombDT[, c("Count.old", "Count.new") := NULL]
    }
  }

  setorder(UniqueEnvVarsCombDT, -Count)
  UniqueEnvVarsCombDT[, frequency := Count / sum(Count)]

  # --- Print raw limits summary
  message("\nOverall limits across time and space (before rounding):")
  print(limits)

  # Save results
  out_path <- file.path(
    DirectoryOutput,
    paste0("UniqueEnvVarsComb_t", firstTimeStep, "-t", lastTimeStep, ".csv")
  )
  write_csv(as.data.frame(UniqueEnvVarsCombDT), out_path)
}

main()
