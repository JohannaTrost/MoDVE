#' Add microclimate data to microhabitat matrices (script entry point)
#'
#' This script integrates selected microclimate variables (humidity,
#' temperature, wind) into time-dependent microhabitat matrices.
#' It is intended to be executed as a standalone script.
#'
#' @details
#' The script expects one command-line argument:
#' \itemize{
#' \item \code{configFile}: path to a configuration file readable by
#'   \code{read.config()}
#' }
#'
#' @keywords internal
NULL

source("utils.R")

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
      rep(top_layer, each = height_diff),
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

main <- function() {

  config <- parse_config()

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
    microhabitat_fname <- paste("microhabitatMatrix", TimeStep, ".rds", sep="")
    FileMatrix <- file.path(Directorymicrohabitat, microhabitat_fname)
    microhabitat <- readRDS(FileMatrix)

    # Load microclimate matrix
    year <- yearStart + TimeStep - timeStepStart
    microclimate_fname <- paste0(year, "_", region, "_mc_matrix.rds")
    FileMcMatrix <- file.path(Directorymicroclimate, microclimate_fname)
    microclimate <- readRDS(FileMcMatrix)

    # 1. Select relevant microclimate variables (humidity=7, temperature=1, windspeed=11)
    microclimate_flags <- c(
      config$HumNicheOpt,
      config$TempNicheOpt,
      config$WindNicheOpt
    ) == 1

    # mean annual humidity (7), temperature (1), wind (11)
    microclimate <- microclimate <- microclimate[,,, c(7, 1, 11)[microclimate_flags]]

    # 2. Fill up above canopy microclimate given microhabitat matrix height
    microhab_dims <- dim(microhabitat)
    microclim_dims <- dim(microclimate)

    microclimate <- adapt_microclimate_height(
      microclimate,
      microhab_dims,
      microclim_dims
    )

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

main()