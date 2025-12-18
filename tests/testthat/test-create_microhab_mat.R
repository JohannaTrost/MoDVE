test_that("abuse cases", {
  # here try to enter incorrect input esp. for shoot and trunk tables
})

# Diameter is constant through tests, doesn't influence
# which voxels the branch goes through etc.
shoot_diameter <- 0.01

# Shorthand functions
create_empty_shoot_tbl <- function() {
  return(tibble::tibble(
    "xbegin" = numeric(),
    "ybegin" = numeric(),
    "zbegin" = numeric(),
    "xend" = numeric(),
    "yend" = numeric(),
    "zend" = numeric(),
    "length" = numeric(),
    "diameter" = numeric(),
    "shootID" = numeric()
  ))
}
add_shoot_row <- function(shoot_tbl, begin_coords, end_coords) {
  return(shoot_tbl |>
    tibble::add_row(
      "xbegin" = begin_coords[1], "xend" = end_coords[1],
      "ybegin" = begin_coords[2], "yend" = end_coords[2],
      "zbegin" = begin_coords[3], "zend" = end_coords[3],
      "length" = sqrt((xend - xbegin)^2 + (yend - ybegin)^2 + (zend - zbegin)^2),
      "diameter" = shoot_diameter,
      "shootID" = nrow(shoot_tbl) + 1
    ))
}

create_empty_trunk_tbl <- function() {
  return(tibble::tibble(
    "x" = numeric(),
    "y" = numeric(),
    "height" = numeric(),
    "diameter" = numeric(),
    "treeID" = numeric()
  ))
}

test_that("assign surface area", {

  # A 5*5 landscape with a 1-cell corridor around it
  corridor <- 1
  dim <- 5
  grid_dim <- dim + 2 * corridor
  microhab_extent <- c(corridor, dim + corridor)

  config <- list(
    TotalSurfaceAreaOpt = 1,
    SurfaceAreaLossOpt = 0,
    LightConditionsOpt = 0,
    AverageWeightedAngles = 0,
    MaxX = dim,
    MaxY = dim,
    MaxZ = 1, # 2D
    corridor = corridor
  )

  # Table template
  {
    shoots_dt <- create_empty_shoot_tbl()
    z <- 0.5
    intersctd_voxels <- list()

    # Along x
    # Case 1: along x: segment starts and ends in corridor, beyond microhabitat area
    begin_coords <- c("x" = 0.5, "y" = microhab_extent[1] + 0.5, "z" = z)
    end_coords <- c("x" = grid_dim, "y" = microhab_extent[1] + 0.5, "z" = z)
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[1]] <- find_intersecting_voxels(begin_coords, end_coords)

    # Case 2: along y: segment goes in descending direction
    # intersects the first segment in [1, 1]
    begin_coords <- c("x" = microhab_extent[1] + 0.5, "y" = microhab_extent[2] - 0.5, "z" = z)
    end_coords <- c("x" = microhab_extent[1] + 0.5, "y" = microhab_extent[1] + 0.5, "z" = z)
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[2]] <- find_intersecting_voxels(begin_coords, end_coords)

    # Case 3: segment is smaller than a cell
    # intersects the first segment
    begin_coords <- c("x" = microhab_extent[2] - 0.8, "y" = microhab_extent[2] - 0.8, "z" = z)
    end_coords <- c("x" = microhab_extent[2] - 0.2, "y" = microhab_extent[2] - 0.2, "z" = z)
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[3]] <- find_intersecting_voxels(begin_coords, end_coords)

    # Case 4 and 5: segments on cell limits belong to the cell below
    begin_coords <- c("x" = 1.5, "y" = 1, "z" = z)
    end_coords <- c("x" = 2.5, "y" = 1, "z" = z)
    # this one doesn't contribute to microhabitat
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[4]] <- find_intersecting_voxels(begin_coords, end_coords)
    begin_coords <- c("x" = 2.5, "y" = 6, "z" = z)
    end_coords <- c("x" = 3.5, "y" = 6, "z" = z)
    # but this one does
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[5]] <- find_intersecting_voxels(begin_coords, end_coords)
  }

  # Landscape viz
  shoots_dt |>
    ggplot2::ggplot() +
    ggplot2::geom_rect(
      xmin = microhab_extent[1], xmax = microhab_extent[2],
      ymin = microhab_extent[1], ymax = microhab_extent[2],
      alpha = 0.1
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = xbegin, y = ybegin, xend = xend, yend = yend, colour = as.factor(shootID))
    ) +
    ggplot2::coord_cartesian(xlim = c(0, grid_dim), ylim = c(0, grid_dim)) +
    ggplot2::theme_linedraw()

  # Set expectations
  nb_voxels <- sapply(intersctd_voxels, length)
  surf_area_exptd <- shoots_dt$length * shoots_dt$diameter * pi / 2 / nb_voxels
  expected_mat <- matrix(data = 0, nrow = dim, ncol = dim)
  expected_mat[1:5, 1] <- surf_area_exptd[1]
  expected_mat[1, 1:5] <- expected_mat[1, 1:5] + surf_area_exptd[2]
  expected_mat[5, 5] <- surf_area_exptd[3]
  expected_mat[2:3, 5] <- surf_area_exptd[5]

  # Carry out test
  microhab_mat <- create_microhabitat_mat(
    config = config,
    shoot_dt = shoots_dt,
    trunk_dt = create_empty_trunk_tbl()
  )
  expect_equal(microhab_mat[,,1,1], expected_mat)

})
