#' Title
#'
#' @param E
#' @param Microhabitat
#' @param SurfaceBiomassScaling
#'
#' @returns
#' @export
#'
#' @examples
resolve_growth <- function(E, Microhabitat, SurfaceBiomassScaling) {

  for (i in seq_len(nrow(E))) {

    vox <- Microhabitat[E$X[i], E$Y[i], E$Z[i], ]

    # maybe it is faster if I do not use the if statement => speed testing
    if (E$Status[i] == 1) {

      # Von Bertalanffy growth function
      growth_term <- E$GrowthRate[i] * (E$MaximumMass[i] - E$Mass[i])

      # Parabolic light response
      light_vox <- vox[3]
      light_term <- E$LightResponseA[i] * light_vox^2 +
        E$LightResponseB[i] * light_vox +
        E$LightResponseC[i]

      E$Mass[i] <- E$Mass[i] + max(0, growth_term * light_term)
    }

    # Add info about the voxel to the epiphyte matrix
    E$SurfaceAreaOccupied[i] <- (E$Mass[i]^(2/3)) / SurfaceBiomassScaling
    # TODO: do we really need individual-level copies of these habitat values?
    E$TotalSurfaceInVoxel[i] <- vox[1]  # Total surface in voxel
    E$SurfaceLossInVoxel[i] <- vox[2]  # Percentage surface loss in this year
    E$LightInVoxel[i] <- vox[3]  # Light conditions in voxel
  }

  return(E)
}
