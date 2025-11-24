
#' Draw the trait value and location of initial individuals
#'
#' Draw the initial mass, age, occupied surface area of initial individuals
#' and distirbute them along eligible voxels in the habitat matrix.
#'
#' @param distr_params a list containing at least the following parameters:
#' \itemize{
#'    \item IndividualsPerSpecies
#'    \item ScalingPerHa
#'    \item PercentageMaturePerSpecies
#'    \item SurfaceBiomassScaling
#' }
#' @param species_df a data frame containing the species traits, as created by
#' `draw_species_traits`
#' @param microhab_mat a matrix containing the surface area and light conditions
#' in each voxel, as created by `create_microhab_mat()`
#' @param path_to_output a string, where (folder and name) should the output
#' be saved? Must be `.csv`.
#' @param largest_inds_first logical, should larger individuals be allocated
#' first (`TRUE`) to represent competitive advantage, or randomly (`FALSE`)
#' @param largest_voxels_first logical, should individuals first be allocated to
#' the largest voxel available (`TRUE`) or a random one (`FALSE`)?
#'
#' @export
#'
draw_initial_individuals <- function(distr_params, species_df,
                                        microhab_mat, path_to_output,
                                        largest_inds_first = TRUE,
                                        largest_voxels_first = TRUE
                                        ) {

  check_distr_params(distr_params)
  check_species_df(species_df)

  dimPlot <- dim(Microhabitat)[1:3]
  if (ScalingPerHa) {
    distr_params$IndividualsPerSpecies <- distr_params$IndividualsPerSpecies *
      ((dimPlot[1] * dimPlot[2]) / 10000)
  }

  NumberSpecies <- length(SpeciesPool$SpeciesID)
  TotalIndividuals <- NumberSpecies * distr_params$IndividualsPerSpecies
  NumberMaturesPerSpecies <- round(
    distr_params$IndividualsPerSpecies * distr_params$PercentageMaturePerSpecies / 100
  )

  # Initialize individual matrix
  ColumnHeaders <- c("X", "Y", "Z", "Mass", "Status", "IndividualID",
                     "SurfaceAreaOccupied", "Age", "SpeciesID")
  init_ind_mat <- array(
    rep(0, TotalIndividuals * length(ColumnHeaders)),
    dim = c(TotalIndividuals, length(ColumnHeaders))
  )
  # Pre-assign matrix columns to variables
  col_x <- match("X", ColumnHeaders)
  col_y <- match("Y", ColumnHeaders)
  col_z <- match("Z", ColumnHeaders)
  col_mass <- match("Mass", ColumnHeaders)
  col_status <- match("Status", ColumnHeaders)
  col_id <- match("IndividualID", ColumnHeaders)
  col_sa <- match("SurfaceAreaOccupied", ColumnHeaders)
  col_age <- match("Age", ColumnHeaders)
  col_sp <- match("SpeciesID", ColumnHeaders)

  # Individual ID
  init_ind_mat[, col_id] <- seq_len(TotalIndividuals)

  # Species depdt qualities
  for (sp in seq_len(NumberSpecies)) {

    first_row <- (sp - 1) * distr_params$IndividualsPerSpecies + 1
    last_row <- distr_params$IndividualsPerSpecies * sp

    # Species ID
    spID <- SpeciesPool$SpeciesID[sp]
    init_ind_mat[first_row:last_row, col_sp] <- spID

    # Mass
    massMaturity <- SpeciesPool$MassAtMaturity[sp]
    maxMass <- SpeciesPool$MaximumMass[sp]
    row_matures <- first_row + NumberMaturesPerSpecies
    init_ind_mat[first_row:row_matures, col_mass] <-runif(
      NumberMaturesPerSpecies, min = massMaturity, max = maxMass
    )
    init_ind_mat[(row_matures+1):last_row, col_mass] <-runif(
      NumberMaturesPerSpecies, min = 0, max = maxMass
    )

    # Age
    growthRate <- SpeciesPool$GrowthRate[sp]
    init_ind_mat[first_row:last_row, col_age] <- round(
      -log(1 - (init_ind_mat[first_row:last_row, col_mass] / maxMass)) / growthRate
    )
  }

  # Surface area occupied
  init_ind_mat[, col_sa] <- (init_ind_mat[, col_mass]^(2 / 3)) /
    distr_params$SurfaceBiomassScaling

  # Allocate voxels (x, y, z, status)

  # Available surface decreases with each allocated individual
  AvailableSurfaceArea <- Microhabitat[, , , 1]

  # Schedule individual allocation priority
  if (largest_inds_first) # largest individuals win competition
    ind_queue <- order(init_ind_matSub[, col_mass], decreasing = TRUE)
  else # random
    ind_queue <- sample(seq_len(TotalIndividuals), TotalIndividuals, replace=FALSE)

  for (ind in ind_queue) {

    sp <- init_ind_mat$SpeciesID[ind]

    # Find all suitable voxels for this individual
    MinLight <- init_ind_mat[ind, ColMinLight]
    MaxLightInd <- init_ind_mat[ind, ColMaxLight]
    AreaNeededInd <- init_ind_mat[NumIndRand, col_sa]

    # 1. Get the positions of all voxels fullfilling the
    # requirements of the individual (light+area)
    hasEnoughSurface <- AvailableSurfaceArea[, , ] > AreaNeededInd
    hasEnoughLight <- Microhabitat[, , , 3] >= SpeciesPool$MinLight[sp]
    hasEnoughShade <- Microhabitat[, , , 3] <= SpeciesPool$MaxLight[sp]
    SuitableVoxels <- which(hasEnoughLight & hasEnoughShade & hasEnoughSurface)

    if (length(SuitableVoxels) > 0) {

      # Select one voxel
      if (largest_voxels_first) { # voxel with the most available surface area
        whichVoxel <- which(AvailableSurfaceArea[SuitableVoxels] ==
                              max(AvailableSurfaceArea[SuitableVoxels]))[1]
      } else { # random voxel
        whichVoxel <- sample.int(length(SuitableVoxels), size=1)
      }
      ids <- arrayInd(SuitableVoxels[whichVoxel], dim(AvailableSurfaceArea))
      x <- ids[, 1]
      y <- ids[, 2]
      z <- ids[, 3]

      # Update available Surface Area
      AvailableSurfaceArea[x, y, z] <- AvailableSurfaceArea[x, y, z] - AreaNeededInd

      # Update Initial Epiphyte Matrix
      init_ind_mat[NumIndRand, col_x] <- x
      init_ind_mat[NumIndRand, col_y] <- y
      init_ind_mat[NumIndRand, col_z] <- z

      # Set status of individual: status=1 => alive
      init_ind_mat[NumIndRand, col_status] <- 1

    } else { # no suitable voxels

      # Set status of individual: status=2 =>> dead
      init_ind_mat[NumIndRand, col_status] <- 2

      # Set coordinates to 1 (might cause problems in later model if not)
      init_ind_mat[NumIndRand, col_x] <- 1
      init_ind_mat[NumIndRand, col_y] <- 1
      init_ind_mat[NumIndRand, col_z] <- 1
    }
  } # loop individuals

  # Convert to df and save
  init_ind_mat_df <- as.data.frame(init_ind_mat)
  names(init_ind_mat_df) <- ColumnHeaders

  # Save Initial Epiphyte Matrix
  write.csv(init_ind_mat_df, path_to_output, row.names = FALSE)
}
