
test_that("identify correct voxels", {

  start_point <- c(-1, 1, 1) # x, y, z
  end_point <- c(5, 3, -1)
  exptd_voxels <- list(
    start_point,
    c(0, 1, 1),
    c(1, 2, 0),
    c(2, 2, 0),
    c(3, 2, 0),
    c(4, 3, -1),
    end_point
    )
  voxels <- find_intersecting_voxels(start_point, end_point)
  testthat::expect_setequal(voxels, exptd_voxels)

  start_point <- c(-7, 0, -3) # x, y, z
  end_point <- c(2, -5, -1)
  exptd_voxels <- list(
    start_point,
    c(-6, -1, -3),
    c(-5, -1, -3),
    c(-4, -2, -2),
    c(-3, -2, -2),
    c(-2, -3, -2),
    c(-1, -3, -2),
    c(0, -4, -1),
    c(1, -4, -1),
    end_point
  )
  voxels <- find_intersecting_voxels(start_point, end_point)
  testthat::expect_setequal(voxels, exptd_voxels)

})
