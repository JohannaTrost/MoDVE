
#' Title
#'
#' @param config
#' \itemize{
#' \item MicrohabitatType  # 1: real GroIMP forest with dynamics
# 2: static GroIMP forest (only forest at timeStepStart is used)
#' \item kL light extinction coefficient
#' \item DistVoxToConsider How many ring around focal voxel to consider in light model (5 voxels in x and y direction)
#' \item TotalSurfaceAreaOpt which of these parameters should be used and added to microhab variables
#' \item SurfaceAreaLossOpt
#' \item LightConditionsOpt
#' \item AverageWeightedAngles
#' }
#' @param shoot_dt
#' @param trunk_dt
#' @param vox_dt only required if LightConditionsOpt is TRUE
#' @param path_to_output
#' @param dead_branches_id
#' @param dead_trees_id
#'
#' @export
#'
create_microhabitat_mat <- function(config, shoot_dt, trunk_dt, vox_dt = NULL,
                                    path_to_output, dead_branches_id,
                                    dead_trees_id) {
  # Inputs are correct
  check_config(config)

  if (is.character(shoot_dt))
    read.table(shoot_dt, sep = "\t",  header = TRUE, skip = 1)
  check_shoot_dt(shoot_dt)

  if (is.character(trunk_dt))
    read.table(trunk_dt, sep = "\t",  header = TRUE, skip = 8)
  check_trunk_dt(trunk_dt)

  if (config$LightConditionsOpt) {
    if (is.character(vox_dt))
      read.table(vox_dt, sep = "\t",  header = TRUE, skip = 1)
    check_vox_dt(vox_dt)
  }

  # Set dimensions
  MaxX <- config$MaxX
  MaxY <- config$MaxY
  MaxZ <- config$MaxZ
  corridor <- config$corridor
  dimPlot <- c(MaxX, MaxY, MaxZ)
  dimX <- MaxX + 2 * corridor
  dimY <- MaxY + 2 * corridor
  dimZ <- MaxZ
  C <- c(0, 0, 1)  # Vector orthogonal to plane of X and Y

  microhab_mat <- array(
    rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3] * 4),
    dim = c(dimPlot[1], dimPlot[2], dimPlot[3], 4)
  )

  # Element indices of the matrix
  sa_elt <- 1
  sa_loss_elt <- 2
  light_elt <- 3
  angle_elt <- 4

  for (s in 1:nrow(shoot_dt)) {

    seg_len <- shoot_dt$length[s]
    seg_diam <- shoot_dt$length[s]
    seg_start <- c(shoot_dt$xbegin[s], shoot_dt$ybegin[s], shoot_dt$zbegin[s])
    seg_end <- c(shoot_dt$xend[s], shoot_dt$yend[s], shoot_dt$zend[s])
    intersectd_voxels <- find_intersecting_voxels(seg_start, seg_end)

    # Calculate total surface area
    seg_surface_area <- seg_len / length(intersectd_voxels) * seg_diam * pi / 2

    if (config$AverageWeightedAngles) {
      V <- seg_end - seg_start
      alpha <- sum(C * V) /
        (sqrt(V[1]^2 + V[2]^2 + V[3]^2) * sqrt(C[1]^2 + C[2]^2 + C[3]^2))
      # ^ sqrt(C) is always 1, remove?
      shoot_angle <- abs(90 - (acos(alpha) / pi * 180))
    }

    for (v in exptd_voxels) {
      x <- v[1]
      y <- v[2]
      z <- v[3]

      if (config$TotalSurfaceAreaOpt)
        microhab_mat[x, y, z, sa_elt] <- microhab_mat[x, y, z, sa_elt] +
        seg_surface_area

      if (shoot_dt$shoot_ID %in% dead_branches_id)
        microhab_mat[x, y, z, sa_loss_elt] <- microhab_mat[x, y, z, sa_loss_elt] +
        seg_surface_area

      if (config$AverageWeightedAngles) {
        V <- seg_end - seg_start
        alpha <- sum(C * V) /
          (sqrt(V[1]^2 + V[2]^2 + V[3]^2) * sqrt(C[1]^2 + C[2]^2 + C[3]^2))
        shoot_angle <- abs(90 - (acos(alpha) / pi * 180))

        # Calculate weighted angle for the voxel
        # ??? source for this?
        tmp1 <- (microhab_mat[x, y, z, sa_elt] - seg_surface_area) /
          microhab_mat[x, y, z, sa_elt] *
          microhab_mat[x, y, z, angle_elt]

        tmp2 <- seg_surface_area / microhab_mat[x, y, z, sa_elt] * ShootAngle

        microhab_mat[x, y, z, angle_elt] <- tmp1 + tmp2
      }

    }
  }

  for (t in 1:nrow(trunk_dt)) {

    x <- ceiling(trunk_dt$x[j])  # X voxel of tree
    y <- ceiling(trunk_dt$y[j])  # Y voxel of tree
    trunk_height <- trunk_dt$height[j]  # trunk_height of tree
    trunk_diameter <- trunk_dt$diameter[j]  # trunk_diameter of tree

    SurfaceAreaTotal <- 0

    z_seq <- rev(seq_len(ceiling(trunk_height))) # top to bottom
    for (z in z_seq) {

      cone_height <- trunk_height - z + 1  # height of cylinder from top to bottom of voxel
      cone_radius <- trunk_diameter / 2  # radius of cylinder at bottom of voxel

      # Calculate total surface area in voxel for trunks
      SurfaceAreaInVoxel <- pi * cone_radius *
        sqrt(cone_radius^2 + cone_height^2) - SurfaceAreaTotal

      # Update total surface area of cylinder so far (to use in next step)
      SurfaceAreaTotal <- SurfaceAreaTotal + SurfaceAreaInVoxel

      if (config$TotalSurfaceAreaOpt) {
        microhab_mat[x, y, z, sa_elt] <- microhab_mat[x, y, z, sa_elt] +
          SurfaceAreaInVoxel
      }

      # If trunk is lost during this time step, add it to lost surface
      if (trunk_dt$tree_id %in% dead_trees_id) {
        microhab_mat[x, y, z, sa_loss_elt] <- microhab_mat[x, y, z, sa_loss_elt] +
          SurfaceAreaInVoxel
      }

      # Update weighted angle for the voxel
      if (config$AverageWeightedAngles) {

        tmp1 <- (microhab_mat[x, y, z, sa_elt] - SurfaceAreaInVoxel) /
          microhab_mat[x, y, z, sa_elt] * microhab_mat[x, y, z, angle_elt]

        tmp2 <- SurfaceAreaInVoxel / microhab_mat[x, y, z, sa_elt] * 90
        # upright 90 degrees angle assumed

        microhab_mat[x, y, z, angle_elt] <- tmp1 + tmp2
      }

    } # z in z seq
  } # t in trunk set

  # Calculate light conditions in voxels (relative light conditions)
  if (config$LightConditionsOpt ) {

    # Total leaf area in each column
    leaf_area_mat <- array(
      rep(0, dimX * dimY * dimZ),
      dim = c(dimX, dimY, dimZ)
      )

    # Store information on leaf area in matrix
    for (vx in 1:nrow(vox_dt)) {
      x <- vox_dt$x[vx]
      y <- vox_dt$y[vx]
      z <- vox_dt$z[vx]
      leaf_area_mat[x, y, z] <- vox_dt$leafarea[vx]
    }

    # Calculate single column light conditions based on leaf area distribution
    for (x in seq_len(dimX)) {
      for (y in seq_len(dimY)) {
        for (z in seq_len(dimZ)) {
          total_leaf_area <- sum(leaf_area_mat[x, y, z:dimZ])
          microhab_mat[x, y, z, light_elt] <- exp(-kL * total_leaf_area / 10000)
        }
      }
    }

    light_matrix_copy <- microhab_mat[, , , light_elt]

    # Calculate final light conditions by accounting for the light
    # conditions in adjacent voxels
    for (x in int_seq(from = corridor, to = dimX - corridor, by = 1)) {
      for (y in int_seq(from = corridor, to = dimY - corridor, by = 1)) {
        for (z in seq_len(dimZ)) {
          TotalContribution <- 0

          # loop over ring surrounding the focal voxel
          for (xx in int_seq(from = x - DistVoxToConsider, to = x + DistVoxToConsider, by = 1)) {
            for (yy in int_seq(from = y - DistVoxToConsider, to = y + DistVoxToConsider, by = 1)) {

              Ring <- max(abs(xx - x), abs(yy - y))
              Contribution <- 1 / (DistVoxToConsider + 1) *
                1 / max(1, (Ring * 8)) * light_matrix_copy[xx, yy, z]

              TotalContribution <- TotalContribution + Contribution
            }
          }

          microhab_mat[x, y, z, light_elt] <- TotalContribution
        } # z
      } # y
    } # x

  }

}
