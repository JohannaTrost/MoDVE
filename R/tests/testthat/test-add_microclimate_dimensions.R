library(testthat)
source("02_add_microclimate_dimensions.R", local = TRUE)

test_that("add_microclimate_dimensions expands dimensions correctly", {
  microhab_dims <- c(2, 2, 2, 2)
  microclim_dims <- c(2, 2, 2, 1)

  microhabitat <- array(1, dim = microhab_dims)
  microclimate <- array(2, dim = microclim_dims)

  index_list <- c(
    TotalSurfaceAreaOpt = 1,
    LightNicheOpt = 2,
    HumNicheOpt = 3
  )

  result <- add_microclimate_dimensions(
    microhab_dims,
    microclim_dims,
    microhabitat,
    microclimate,
    index_list
  )

  expect_equal(
    dim(result),
    c(2, 2, 2, 3)
  )
})

test_that("non-microclimate variables are copied unchanged", {
  microhabitat <- array(runif(16), dim = c(2, 2, 2, 2))
  microclimate <- array(99, dim = c(2, 2, 2, 1))

  index_list <- c(
    TotalSurfaceAreaOpt = 1,
    LightNicheOpt = 2,
    HumNicheOpt = 3
  )

  result <- add_microclimate_dimensions(
    dim(microhabitat),
    dim(microclimate),
    microhabitat,
    microclimate,
    index_list
  )

  expect_equal(
    result[,,,c(1, 2)],
    microhabitat[,,,c(1, 2)]
  )
})

test_that("microclimate variable is inserted correctly", {
  microhabitat <- array(0, dim = c(2, 2, 2, 1))
  microclimate <- array(5, dim = c(2, 2, 2, 1))

  index_list <- c(
    TotalSurfaceAreaOpt = 1,
    HumNicheOpt = 2
  )

  result <- add_microclimate_dimensions(
    dim(microhabitat),
    dim(microclimate),
    microhabitat,
    microclimate,
    index_list
  )

  expect_equal(
    result[,,,2],
    microclimate[,,,1]
  )
})

test_that("multiple microclimate variables are inserted in correct order", {
  microhabitat <- array(0, dim = c(2, 2, 1, 3))

  microclimate <- array(
    c(10, 20),
    dim = c(2, 2, 1, 2)
  )

  index_list <- c(
    TotalSurfaceAreaOpt = 1,
    SurfaceAreaLossOpt = 2,
    LightNicheOpt = 3,
    HumNicheOpt = 4,
    TempNicheOpt = 5
  )

  result <- add_microclimate_dimensions(
    dim(microhabitat),
    dim(microclimate),
    microhabitat,
    microclimate,
    index_list
  )

  expect_equal(result[,,,4], microclimate[,,,1])
  expect_equal(result[,,,5], microclimate[,,,2])
})

test_that("no microclimate variables copies all data from microhabitat", {
  microhabitat <- array(runif(8), dim = c(2, 2, 2, 1))
  microclimate <- array(runif(20), dim = c(2, 2, 2, 0))

  index_list <- c(TotalSurfaceAreaOpt = 1)

  result <- add_microclimate_dimensions(
    dim(microhabitat),
    dim(microclimate),
    microhabitat,
    microclimate,
    index_list
  )

  expect_equal(
    result[,,,1],
    microhabitat[,,,1]
  )
})

test_that("NA values are preserved", {
  microhabitat <- array(NA_real_, dim = c(2, 2, 1, 1))
  microclimate <- array(NA_real_, dim = c(2, 2, 1, 1))

  index_list <- c(TotalSurfaceAreaOpt = 1, HumNicheOpt = 2)

  result <- add_microclimate_dimensions(
    dim(microhabitat),
    dim(microclimate),
    microhabitat,
    microclimate,
    index_list
  )

  expect_true(all(is.na(result)))
})
