#' Add microclimate data to microhabitat matrices
#'
#' This script integrates selected microclimate variables (humidity, temperature, wind)
#' into time-dependent microhabitat matrices. It is intended to be executed as a standalone script.
#'
#' @details
#' Usage: Rscript model_pipeline/03_add_microclimate_dimensions.R --config config.toml
#'
#' Example config.toml
#'
#' timeStepStart = 1                               # first simulation timestep
#' timeStepEnd = 3                                 # last simulation timestep
#' yearStart = 1980                                # starting calendar year
#' region = "DE"                                   # region identifier for microclimate files
#' Directorymicrohabitat = "/path/to/microhabitat/rep<n>"    # input microhabitat directory
#' Directorymicroclimate = "/path/to/microhabitat_mc/rep<n>" # input microclimate directory
#' Directorynew_microhabitat = "/path/to/output/rep<n>"  # output directory
#' TotalSurfaceAreaOpt = 1                         # {0, 1} include total surface area
#' SurfaceAreaLossOpt = 1                          # {0, 1} include surface area loss
#' LightNicheOpt = 1                               # {0, 1} include light niche variables
#' AverageWeightedAngles = 0                       # {0, 1} include weighted angle averages
#' HumNicheOpt = 1                                 # {0, 1} include humidity niche
#' TempNicheOpt = 1                                # {0, 1} include temperature niche
#' WindNicheOpt = 1                                # {0, 1} include wind niche
#'
#' Output:
#' [MicrohabitatMatrix<timeStepStart>.rds, ..., MicrohabitatMatrix<timeStepEnd>.rds]
#' With matrices xDim X yDim X zDim X nVariables (including microclimate variables)
#' Forest parameter files (Forest_param_global.txt, Forest_param_pass0.txt, dimPlot.rds)
#' are copied to the output directory.
#'
NULL

source("model_pipeline/utils.R")

library(fs)

#' Add microclimate variables to a microhabitat matrix
#'
#' This function merges selected microclimate variables (e.g. humidity,
#' temperature, wind) into an existing microhabitat 4D array by expanding
#' its 4th (variable) dimension. Variables are inserted according to the
#' provided index list, while non-microclimate variables are copied unchanged.
#'
#' @param microhab_dims integer vector giving the dimensions of the original
#' microhabitat array \code{[x, y, height, variables]}
#' @param microclim_dims integer vector giving the dimensions of the
#' microclimate array \code{[x, y, height, variables]}
#' @param microhabitat 4D numeric array containing microhabitat variables
#' @param microclimate 4D numeric array containing microclimate variables
#' @param microhabitat_index_list named integer vector mapping variable names
#' to their positions in the combined microhabitat array
#'
#' @return A 4D numeric array with microclimate variables inserted into the
#' microhabitat matrix
#'
#' @export
add_microclimate_dimensions <- function(microhab_dims,
                                        microclim_dims,
                                        microhabitat,
                                        microclimate,
                                        microhabitat_index_list) {

  newD4 <- microhab_dims[4] + microclim_dims[4]

  # Assert that newD4 matches length of microhabitat_index_list otherwise throw an error
  if (newD4 != length(microhabitat_index_list) + 1) { # PLus 1 because of PAI
      stop("Mismatch between expected new variable dimension and provided index list.")
  }

  new_microhabitat <- array(
    NA_real_,
    dim = c(microhab_dims[1:3], newD4)
  )

  microclimate_var_options <- c( # Possible microclimate index names
    "HumNicheOpt",
    "TempNicheOpt",
    "WindNicheOpt"
  )

  microclimate_index <- 1 # Start index for microclimate variables in microclimate matrix

  for (var_name in names(microhabitat_index_list)) {

    index <- microhabitat_index_list[[var_name]]

    if (!var_name %in% microclimate_var_options) { # If not a microclimate variable, copy as is
      new_microhabitat[,,, index] <- microhabitat[,,, index]
    } else { # If microclimate variable, copy from microclimate matrix
      new_microhabitat[,,, index] <- microclimate[,,, microclimate_index]
      microclimate_index <- microclimate_index + 1
    }
  }

  new_microhabitat
}


#' Adjust microclimate array height to match microhabitat
#'
#' This function fills or trims the vertical (height) dimension of a
#' microclimate array so that it matches the height of a microhabitat array.
#' If the microhabitat is taller, the top microclimate layer is repeated
#' upward. If it is shorter, the microclimate is trimmed.
#'
#' @param microclimate 4D numeric array
#'   \code{[x, y, height, variables]}
#' @param microhab_dims integer vector giving the dimensions of the
#' microhabitat array
#' @param microclim_dims integer vector giving the dimensions of the
#' microclimate array
#'
#' @return A microclimate array with height dimension matching the microhabitat
#'
#' @importFrom abind abind
#' @export
adapt_microclimate_height <- function(microclimate,
                                      microhab_dims,
                                      microclim_dims) {

  height_diff <- microhab_dims[3] - microclim_dims[3] # Difference in vertical dimension

  if (height_diff > 0) { # Fill missing above-canopy climate by repeating top layer

    top_layer <- microclimate[,, microclim_dims[3], , drop = FALSE]

    fill_array <- array(
      rep(top_layer, times = height_diff),
      dim = c(
        microhab_dims[1:2],
        height_diff,
        microclim_dims[4]
      )
    )

    microclimate <- abind::abind(
      microclimate,
      fill_array,
      along = 3
    )

  } else if (height_diff < 0) { # Trim microclimate to match microhabitat height

    microclimate <- microclimate[,, seq_len(microhab_dims[3]), , drop = FALSE]
  }

  microclimate
}

#' Add microclimate data to microhabitat matrices
#'
#' This function loads microhabitat matrices for each timestep, integrates
#' selected microclimate variables (humidity, temperature, wind),
#' adjusts vertical dimensions if needed, and saves updated matrices
#' to a new output directory.
#'
#' @param config A list containing configuration parameters (see script details)
#'
#' @details
#' Microclimate variables are extracted from fixed layers of the
#' microclimate matrix:
#' \itemize{
#' \item Temperature = layer 1
#' \item Humidity = layer 7
#' \item Wind speed = layer 11
#' }
#'
#' The function creates the output directory if necessary and copies
#' required forest parameter files.
#'
#' @return Invisibly returns \code{TRUE} upon successful completion.
#'
#' @import fs
#' @export
add_microclimate_to_microhabitat <- function(config) {

  # Start and end timestep
  timeStepStart <- config$timeStepStart
  timeStepEnd <- config$timeStepEnd
  yearStart <- config$yearStart

  region <- config$region # Region from where we use climate data

  # Directory paths
  Directorymicrohabitat <- config$Directorymicrohabitat
  Directorymicroclimate <- config$Directorymicroclimate
  Directorynew_microhabitat <- config$Directorynew_microhabitat

  # Creat directory for new microhabitat if it doesn't exist
  if (!dir.exists(Directorynew_microhabitat)) {
    dir.create(Directorynew_microhabitat, recursive=TRUE)
  }

  # Define option flags as a named list
  options_list <- list(
    TotalSurfaceAreaOpt   = config$TotalSurfaceAreaOpt,
    SurfaceAreaLossOpt    = config$SurfaceAreaLossOpt,
    LightNicheOpt         = config$LightNicheOpt,
    AverageWeightedAngles = config$AverageWeightedAngles,
    HumNicheOpt           = config$HumNicheOpt,
    TempNicheOpt          = config$TempNicheOpt,
    WindNicheOpt          = config$WindNicheOpt
  )

  # Only keep active options
  active_options <- names(options_list[options_list == 1])

  # Assign microhabitat_index_list
  microhabitat_index_list <- setNames(seq_along(active_options), active_options)

  for (TimeStep in int_seq(from=timeStepStart, to=timeStepEnd, by=1)) {

    # Load  microhabitat matrix
    microhabitat_fname <- paste("MicrohabitatMatrix", TimeStep, ".rds", sep="")
    FileMatrix <- file.path(Directorymicrohabitat, microhabitat_fname)
    if (!file.exists(FileMatrix)) {
      stop("Microhabitat file not found: ", FileMatrix)
    }
    microhabitat <- readRDS(FileMatrix)

    # Load microclimate matrix
    year <- yearStart + TimeStep - timeStepStart
    microclimate_fname <- paste0(year, "_", region, "_mc_matrix.rds")
    FileMcMatrix <- file.path(Directorymicroclimate, microclimate_fname)
    if (!file.exists(FileMcMatrix)) {
      stop("Microclimate file not found: ", FileMcMatrix)
    }
    microclimate <- readRDS(FileMcMatrix)

    # Check if first 2 dimensions of microclimate and microhabitat match
    if (!all(dim(microclimate)[1:2] == dim(microhabitat)[1:2])) {
      stop("Microclimate and microhabitat dimensions do not match in x and y.")
    }

    # 1. Select relevant microclimate variables (humidity=7, temperature=1, windspeed=11)
    microclimate_flags <- c(
      config$HumNicheOpt,
      config$TempNicheOpt,
      config$WindNicheOpt
    ) == 1

    # mean annual humidity (7), temperature (1), wind (11)
    microclimate <- microclimate[,,, c(7, 1, 11)[microclimate_flags], drop = FALSE]

    # 2. Fill up above canopy microclimate given microhabitat matrix height
    microhab_dims <- dim(microhabitat)
    microclim_dims <- dim(microclimate)

    if (any(microclimate_flags)) {
      microclimate <- adapt_microclimate_height(
        microclimate,
        microhab_dims,
        microclim_dims
      )
    }

    # 3. Insert microclimate into microhabitat matrix

    microhabitat_with_climate <- add_microclimate_dimensions(
      microhab_dims,
      microclim_dims,
      microhabitat,
      microclimate,
      microhabitat_index_list
    )

    # 4. Update microhabitat matrix
    Newmicrohabitat_fnameMatrix <- file.path(Directorynew_microhabitat, microhabitat_fname)
    saveRDS(microhabitat_with_climate, file = Newmicrohabitat_fnameMatrix)

    print(paste("Year:", year, "Time step:", TimeStep))
  }

  # Copy necessary forest files
  fs::file_copy(
    file.path(config$Directorymicrohabitat, c(
      "Forest_param_global.txt",
      "Forest_param_pass0.txt",
      "dimPlot.rds"
    )),
    config$Directorynew_microhabitat,
    overwrite = TRUE
  )

  invisible(TRUE)
}

main <- function () {
  config <- parse_config()
  add_microclimate_to_microhabitat(config)
}

if (sys.nframe() == 0) { # If script is run directly (e.g. using source()), execute
  main()
}