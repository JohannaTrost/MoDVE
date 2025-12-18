#' Title
#'
#' @param config a list with \itemize{
#' \item MicrohabitatType  # 1: real GroIMP forest with dynamics
# 2: static GroIMP forest (only forest at timeStepStart is used)
#' \item kL light extinction coefficient
#' \item DistVoxToConsider How many ring around focal voxel to consider in light model (5 voxels in x and y direction)
#' \item TotalSurfaceAreaOpt which of these parameters should be used and added to microhab variables
#' \item SurfaceAreaLossOpt
#' \item LightConditionsOpt
#' \item AverageWeightedAngles
#' }
#' @param shoot_dt matrix with branch information
#' @param trunk_dt matrix with trunk information
#' @param vox_dt only required if LightConditionsOpt is TRUE
#' @param path_to_output string, where to save output? Must be an rds file or `NULL`,
#' in which case the result matrix is returned.
#' @param dead_branches_id integer vector containing the IDs of all branches dying this time step
#' @param dead_trees_id integer vector containing the IDs of all trees dying this time step
#'
#' @export
#'
create_microhabitat_mat <- function(config, shoot_dt, trunk_dt, vox_dt = NULL,
                                    path_to_output = NULL, dead_branches_id = NULL,
                                    dead_trees_id = NULL) {
  # Inputs are correct
  #check_config(config)

  if (is.character(shoot_dt))
    read.table(shoot_dt, sep = "\t",  header = TRUE, skip = 1)
  #check_shoot_dt(shoot_dt)

  if (is.character(trunk_dt))
    read.table(trunk_dt, sep = "\t",  header = TRUE, skip = 8)
  #check_trunk_dt(trunk_dt)

  if (config$LightConditionsOpt) {
    if (is.character(vox_dt))
      read.table(vox_dt, sep = "\t",  header = TRUE, skip = 1)
   # check_vox_dt(vox_dt)
  }

  if (!is.null(path_to_output)) {
    dir_output <- dirname(path_to_output)
    if (!dir.exists(dirname(dir_output))) {
      stop(paste0("Output directory ", dir_output, " does not exist."))
    }
    if (!grepl("*.rds$", path_to_output)) {
      stop("path_to_output must be a rds file")
    }
  }

  # Set dimensions
  MaxX <- config$MaxX
  MaxY <- config$MaxY
  MaxZ <- config$MaxZ
  corridor <- config$corridor
  dimPlot <- c(MaxX, MaxY, MaxZ)
  forest_max_x <- MaxX + 2 * corridor
  forest_max_y <- MaxY + 2 * corridor

  microhab_mat <- array(
    rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3] * 4),
    dim = c(dimPlot[1], dimPlot[2], dimPlot[3], 4)
  )

  # Element indices of the matrix
  sa_elt <- 1
  sa_loss_elt <- 2
  light_elt <- 3
  angle_elt <- 4

  for (s in seq_len(nrow(shoot_dt))) {

    seg_len <- shoot_dt$length[s]
    seg_diam <- shoot_dt$diameter[s]
    seg_start <- c(shoot_dt$xbegin[s] - corridor,
                   shoot_dt$ybegin[s] - corridor,
                   shoot_dt$zbegin[s])
    seg_end <- c(shoot_dt$xend[s] - corridor,
                 shoot_dt$yend[s] - corridor,
                 shoot_dt$zend[s])
    intersectd_voxels <- find_intersecting_voxels(seg_start, seg_end)

    # Calculate total surface area and split it evenly across intersected voxels
    seg_surface_area <- seg_len * seg_diam * pi / 2 / length(intersectd_voxels)

    # Drop voxels in the corridor
    is_in_microhab_area <- sapply(
      intersectd_voxels,
      function(v) v[1] <= MaxX && v[2] <= MaxY &&
        v[1] > 0 && v[2] > 0
      )
    intersectd_voxels <- intersectd_voxels[is_in_microhab_area]

    for (v in intersectd_voxels) {
      x <- v[1]
      y <- v[2]
      z <- v[3]
      voxel <- microhab_mat[x, y, z, ]

      if (config$TotalSurfaceAreaOpt) {
        microhab_mat[x, y, z, sa_elt] <- voxel[sa_elt] + seg_surface_area
      }

      if (config$SurfaceAreaLossOpt && shoot_dt$shootID[s] %in% dead_branches_id) {
        microhab_mat[x, y, z, sa_loss_elt] <- voxel[sa_loss_elt] +
        seg_surface_area
      }

      if (config$AverageWeightedAngles) {
        V <- seg_end - seg_start
        alpha <- V[1] / sqrt(V[1]^2 + V[2]^2 + V[3]^2)
        shoot_angle <- abs(90 - (acos(alpha) / pi * 180))

        # Calculate weighted angle for the voxel
        # ??? source for this?
        curr_surface_area <- microhab_mat[x, y, z, sa_elt]
        tmp1 <- (curr_surface_area - seg_surface_area) /
          curr_surface_area * voxel[angle_elt]

        tmp2 <- seg_surface_area / curr_surface_area * shoot_angle

        microhab_mat[x, y, z, angle_elt] <- tmp1 + tmp2
      }

    }
  }

  for (t in seq_len(nrow(trunk_dt))) {

    x <- ceiling(trunk_dt$x[t]) - corridor
    y <- ceiling(trunk_dt$y[t]) - corridor
    trunk_height <- trunk_dt$height[t]
    trunk_diameter <- trunk_dt$diameter[t]

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
      if (trunk_dt$treeID[s] %in% dead_trees_id) {
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
  if (config$LightConditionsOpt) {

    light_range <- config$DistVoxToConsider

    # Total leaf area in each column
    # Must process voxels in the corridor too as they affect neighbouring voxels
    leaf_area_mat <- light_mat <- array(
      rep(0, forest_max_x * forest_max_y * MaxZ),
      dim = c(forest_max_x, forest_max_y, MaxZ)
      )

    # Store information on leaf area in matrix
    for (vx in seq_len(nrow(vox_dt))) {
      x <- vox_dt$x[vx]
      y <- vox_dt$y[vx]
      z <- vox_dt$z[vx]
      leaf_area_mat[x, y, z] <- vox_dt$leafarea[vx]
    }

    # Calculate single column light conditions based on leaf area distribution
    for (x in seq_len(forest_max_x)) {
      for (y in seq_len(forest_max_y)) {
        for (z in seq_len(MaxZ)) {
          total_leaf_area <- sum(leaf_area_mat[x, y, z:MaxZ])
          light_mat[x, y, z] <- exp(-config$kL * total_leaf_area / 10000)
        }
      }
    }

    # Calculate final light conditions by accounting for the light
    # conditions in adjacent voxels
    # x and y are indices in the full matrix including corridors
    for (x in seq(from = corridor + 1, to = forest_max_x - corridor)) {
      for (y in seq(from = corridor + 1, to = forest_max_y - corridor)) {
        for (z in seq_len(MaxZ)) {
          TotalContribution <- 0

          # loop over ring surrounding the focal voxel
          xx_seq <- seq(from = x - light_range, to = x + light_range)
          yy_seq <- seq(from = y - light_range, to = y + light_range)
          for (xx in xx_seq) {
            for (yy in yy_seq) {
              Ring <- max(abs(xx - x), abs(yy - y))
              Contribution <- 1 / (light_range + 1) / max(1, (Ring * 8)) *
                light_mat[xx, yy, z]

              TotalContribution <- TotalContribution + Contribution
            }
          }
          microhab_mat[x - corridor, y - corridor, z, light_elt] <- TotalContribution
        } # z
      } # y
    } # x

  } # lightConditionsOpt

  if (!is.null(path_to_output)) {
    saveRDS(microhab_mat, path_to_output)
  } else {
    return(microhab_mat)
  }
}
