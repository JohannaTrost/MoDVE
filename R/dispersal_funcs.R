
#' Extract individual coordinates in vector form
#'
#' Reads a count matrix and extracts the x, y and z coordinates of all
#' individuals
#'
#' @param count_matrix a 3d integer matrix containing tallies of individuals
#' present in each cell
#'
#' @returns a list of 3 integer vectors `x`, `y` and `z` containing the
#' corresponding coordinates of all individuals
#' @export
#'
extract_ind_coords <- function(count_matrix) {

  nb_inds <- sum(count_matrix)

  ids <- arrayInd(
    which(count_matrix > 0),
    dim(count_matrix)
  )
  x_inds <- ids[, 1]
  y_inds <- ids[, 2]
  z_inds <- ids[, 3]

  while (nb_inds > length(x_inds)) {

    # recursively distribute coordinates until counts of remaining
    # unprocessed individuals reaches zero
    tmp_ids <- arrayInd(
      which(count_matrix > 0),
      dim(count_matrix)
    )

    # decrement count
    count_matrix[tmp_ids] = count_matrix[tmp_ids] - 1

    tmp_ids <- arrayInd(
      which(count_matrix > 0),
      dim(count_matrix)
    )

    x_inds <- append(x_inds, tmp_ids[, 1])
    y_inds <- append(y_inds, tmp_ids[, 2])
    z_inds <- append(z_inds, tmp_ids[, 3])
  }

  return(list(x  = x_inds, y = y_inds, z = z_inds))
}

#' Resolve the dispersal step of the simulation
#'
#' @param E
#' @param Microhabitat
#' @param SurfaceBiomassScaling
#' @param centralPoint
#' @param InterceptRecruitment
#' @param SlopeRecruitment
#' @param prob_disp_matrix
#' @param SpeciesPool
#' @param max_id
#'
#' @returns a list
#' @export
#'
resolveReproDispersal <- function(E,
                             Microhabitat,
                             SurfaceBiomassScaling,
                             centralPoint,
                             InterceptRecruitment,
                             SlopeRecruitment,
                             prob_disp_matrix,
                             SpeciesPool,
                             max_id) {

  dimPlot <- dim(Microhabitat)[1:3]
  empty_3d_array <- array(
    rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3]),
    dim = c(dimPlot[1], dimPlot[2], dimPlot[3])
  )

  # Store number of individuals at beginning of time step
  NumberOfSpecies <- nrow(SpeciesPool)
  nbIndsBeforeDisp <- array(rep(0, NumberOfSpecies))
  for (sp in seq_len(NumberOfSpecies)) {
    nbIndsBeforeDisp[sp] <- length(which(E$SpeciesID == sp & E$Status == 1))
  }

  # Calculate available surface area per voxel
  avail_sa_matrix <- Microhabitat[, , , 1]
  for (i in seq_len(nrow(E))) {
    sa_needed <- E$Mass[i]^(2/3) / SurfaceBiomassScaling
    avail_sa_matrix[E$X[i], E$Y[i], E$Z[i]] <- max(
      0, avail_sa_matrix[E$X[i], E$Y[i], E$Z[i]] - sa_needed
    )
  }

  # Initialize potential recruitment dataframe
  unique_species <- unique(E$SpeciesID)
  recruitment_df <- data.frame(
    species_index = seq_len(NumberOfSpecies),
    exptd_nb_recruits = numeric(NumberOfSpecies),
    nb_recruits = numeric(NumberOfSpecies)
  )

  # Loop over all species
  for (i in seq_len(NumberOfSpecies)) {

    sp <- unique_species[i]

    # Generate initially empty matrix to store the probabilities for recruitment
    exptd_nb_recruits_matrix <- empty_3d_array

    # Matrix containing all mature individuals of one species
    mature_inds <- E[E$SpeciesID == sp & E$Mass >= E$MassAtMaturity, ]

    minLight <- SpeciesPool$MinLight[i]
    maxLight <- SpeciesPool$MaxLight[i]

    # Probability matrix for each species:
    # Depending on the position of each mature individual,
    # the total probability for the species is calculated.
    #
    # The second part of the equation accounts for the actual size of the individual
    # in relation to the maximum size for which the recruitment per individual is defined
    for (j in seq_len(nrow(mature_inds))) {

      dist_to_center <- c(
        centralPoint[1] - mature_inds$X[j],
        centralPoint[2] - mature_inds$Y[j],
        centralPoint[3] - mature_inds$Z[j]
      )
      x_coords <- seq(dist_to_center[1] + 1, dist_to_center[1] + dimPlot[1])
      y_coords <- seq(dist_to_center[2] + 1, dist_to_center[2] + dimPlot[2])
      z_coords <- seq(dist_to_center[3] + 1, dist_to_center[3] + dimPlot[3])

      # TODO: none of the following terms are documented in the paper?
      # Mass-dependent fecundity?
      factor1 <- (InterceptRecruitment + SlopeRecruitment * mature_inds$Mass[j]) *
        mature_inds$RecruitmentInvestmentRel[j]
      # ??
      factor2 <- (mature_inds$Mass[j] - mature_inds$MassAtMaturity[j]) /
        (mature_inds$MaximumMass[j] - mature_inds$MassAtMaturity[j])
      # ??
      factor3 <- 1 + (mature_inds$RecruitmentInc[j] * factor2)

      # Dispersal probability * fecundity = expected nb offspring in each xyz
      exptd_nb_recruits_matrix <- exptd_nb_recruits_matrix +
        prob_disp_matrix[x_coords, y_coords, z_coords, sp] * factor1 * factor3
    }

    # Store potential normalized number of recruits
    # We will use this to populate sp_output_mat later
    recruitment_df$exptd_nb_recruits[i] <- sum(exptd_nb_recruits_matrix)

    # Matrix containing all voxel for which the light requirements are fulfilled
    pot_hab_matrix <- ifelse(
      Microhabitat[, , , 3] >= minLight &
        Microhabitat[, , , 3] <= maxLight,
      1, 0
    )

    # Disable unsuitable cells and scale with surface area
    exptd_nb_recruits_matrix <- exptd_nb_recruits_matrix *
      pot_hab_matrix * avail_sa_matrix # TODO: confirm fecundity scales with SA?

    # Calculate number of recruits based on final probability matrix
    nb_recruits_matrix <- array(
      rpois(length(exptd_nb_recruits_matrix), exptd_nb_recruits_matrix),
      dim = dim(exptd_nb_recruits_matrix)
    )

    # Increment recruit counts
    totalNbRecruits <- sum(nb_recruits_matrix)
    recruits_tbl$nb_recruits[sp] <- totalNbRecruits

    if (totalNbRecruits > 0) {

      # Add new recruits to epiphyte matrix
      recruit_coords <- extract_ind_coords(nb_recruits_matrix)
      idx_recruits <- seq(nrow(E) + 1, nrow(E) + totalNbRecruits)
      E[idx_recruits, names(SpeciesPool)] <- SpeciesPool[sp, ]
      E$X[idx_recruits] <- recruit_coords$x
      E$Y[idx_recruits] <- recruit_coords$y
      E$Z[idx_recruits] <- recruit_coords$z
      E$Mass[idx_recruits] <- 0  # Initial size
      E$Status[idx_recruits] <- 1  # status 1:alive
      recruits_ids <- seq(max_id + 1, max_id + length(totalNbRecruits))
      E$IndividualID[idx_recruits] <- recruits_ids

      E[is.na(E)] <- 0  # convert all NA to 0 so that the R script matches the Matlab
      # TODO: this is likely to cause bugs

      max_id <- max_id + length(x_recruits)
    } # if any recruits

  } # species loop

  disp_items <- list(
    "nbIndsBeforeDisp" = nbIndsBeforeDisp,
    "E" = E,
    "recruitment_df" = recruitment_df,
    "max_id" = max_id
  )
  return(disp_items)
}

