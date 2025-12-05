#' Title
#'
#' @param centralPoint a numeric vector of length 3 containing the X, Y and Z
#' coordinates of the center
#' @param dimX
#' @param dimY
#' @param dimZ
#' @param NumberOfSpecies
#' @param SpeciesPool
#'
#' @returns
#' @export
#'
#' @examples
calc_prob_disp_matrix <- function(centralPoint, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool) {

  # Calculate distance to central point
  DistanceMatrix <- array(
    rep(0, dimX * dimY * dimZ),
    dim = c(dimX, dimY, dimZ)
  )
  for (i in seq_len(dimX)) {
    for (j in seq_len(dimY)) {
      for (k in seq_len(dimZ)) {
        x1 <- c(i, j, k)
        x2 <- c(centralPoint[1], centralPoint[2], centralPoint[3])
        DistanceMatrix[i, j, k] <- sqrt(sum((x1 - x2)^2))  # call to pdist() in the matlab script
      }
    }
  }

  # Get probabilities to disperse in each voxel
  ProbabilityMatrix <- prob_disp_matrix <- array(
    rep(0, dimX * dimY * dimZ * NumberOfSpecies),
    dim = c(dimX, dimY, dimZ, NumberOfSpecies)
  )

  for (i in seq_len(NumberOfSpecies)) {

    # Negative exponential
    exponentE <- SpeciesPool$DispersalKernel[i]
    ProbabilityMatrix[, , , i] <- exp(-DistanceMatrix * exponentE)

    # Dispersal asymmetry (probability to disperse downwards > upwards)
    dispersalAsymmetry <- SpeciesPool$DispersalKernelAsymmetry[i]
    z_seq_up <- int_seq(from = centralPoint[3], to = dimZ, by = 1)
    z_seq_down <- int_seq(from = 1, to = centralPoint[3] - 1, by = 1)
    ProbabilityMatrix[, , z_seq_up, i] <- ProbabilityMatrix[, , z_seq_up, i] *
      2 * (1 - dispersalAsymmetry)
    ProbabilityMatrix[, , z_seq_down, i] <- ProbabilityMatrix[, , z_seq_down, i] *
      2 * dispersalAsymmetry

    # Normalize
    prob_disp_matrix[, , , i] <- ProbabilityMatrix[, , , i] /
      sum(ProbabilityMatrix[, , , i])
  }
  return(prob_disp_matrix)
}

#' Find a the central point in a three dimension coordinate system
#'
#' Given dimensions X, Y and Z, returns the 3d coordinates of the central point
#'
#' @param dims a length-3 numeric vector containing dimension sizes X, Y and Z
#' of the coordinate system
#'
#' @returns a length-3 numeric vector containing the X, Y and Z coordinates of
#' the central point
#' @export
#'
find_central_point <- function(dims) {
  dimX <- dims[1] * 2 + 1
  dimY <- dims[2] * 2 + 1
  dimZ <- dims[3] * 2 + 1
  return(c(
    floor(dimX/2) + 1,
    floor(dimY/2) + 1,
    floor(dimZ/2) + 1
  ))
}
