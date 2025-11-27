#' Title
#'
#' @param E
#' @param Microhabitat
#' @param CompetitionMethod
#'
#' @returns
#' @export
#'
#' @examples
resolve_competition <- function(E, Microhabitat, CompetitionMethod) {

  # Calculate total surface area occupied by epiphytes per voxel
  dims <- dim(Microhabitat)
  occupiedSA <- array(
    rep(0, dims[1] * dims[2] * dims[3]),
    dim = c(dims[1], dims[2], dims[3])
  )
  for (w in seq_len(nrow(E))) {
    if (E$Status[w] == 1) {
      occupiedSA[E$X[w], E$Y[w], E$Z[w]] <-
        occupiedSA[E$X[w], E$Y[w], E$Z[w]] +
        E$SurfaceAreaOccupied[w]
      # for each cell, sum sa over all inds in this cell
    }
  }

  # Oversaturated voxels: occupied S.A. exceeds available S.A.
  oversat_voxels <- arrayInd(
    which(occupiedSA > Microhabitat[, , , 1]),
    dim(occupiedSA)
  )
  vox_xs <- oversat_voxels[, 1]
  vox_ys <- oversat_voxels[, 2]
  vox_zs <- oversat_voxels[, 3]

  for (i in seq_len(length(vox_xs))) {

    # Get all individuals in thisvoxel
    isInVoxel <- E$X == vox_xs[i] & E$Y == vox_ys[i] & E$Z == vox_zs[i]
    indsInVoxel <- E[isInVoxel & E$Status == 1, ]

    # Sort individuals
    if (CompetitionMethod == 1) { # priority to larger individuals
      ind_seq <- order(indsInVoxel$SurfaceAreaOccupied, decreasing = TRUE)
    } else { # random
      ind_seq <- sample(seq_len(nrow(indsInVoxel)))
    }
    indsInVoxel <- indsInVoxel[ind_seq, ]

    totalSurfAreaOcc <- cumsum(indsInVoxel$SurfaceAreaOccupied)
    availSurfArea <- Microhabitat[vox_xs[i], vox_ys[i], vox_zs[i], 1]

    # The voxel can support this many individuals
    capacity <- length(which(totalSurfAreaOcc <= availSurfArea))

    if (capacity < nrow(indsInVoxel)) {
      # largest n individuals live, the rest die
      seq_beyond_capacity <- int_seq(capacity + 1, nrow(indsInVoxel))
      dead_ids <- indsInVoxel[seq_beyond_capacity, "IndividualID"]
      these_die <- is.element(E$IndividualID, dead_ids)
      E[these_die, "Status"] <- 2
    }
  }
  return(E)
}
