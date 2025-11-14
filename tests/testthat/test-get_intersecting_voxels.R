
test_that("identify correct voxels", {
  # TODO: z is constant here but a robust use case should include change along z
  start_point <- c(1, 1, 1) # x, y, z
  end_point <- c(2, 3, 1)
  exptd_voxels <- list(c(1, 1, 1), c(1, 2, 1), c(2, 2, 1), c(2, 3, 1))
  voxels <- get_intersecting_voxels(start_point, end_point)
  testthat::expect_setequal(voxels, exptd_voxels)
})
