library(testthat)
source("02_add_microclimate_dimensions.R", local = TRUE)

test_that("no height difference returns identical array", {
  microclimate <- array(runif(8), dim = c(2, 2, 2, 1))

  result <- adapt_microclimate_height(
    microclimate,
    microhab_dims = c(2, 2, 2, 1),
    microclim_dims = dim(microclimate)
  )

  expect_equal(result, microclimate)
})

test_that("microclimate height is expanded by repeating top layer", {
  microclimate <- array(1:8, dim = c(2, 2, 2, 1))

  result <- adapt_microclimate_height(
    microclimate,
    microhab_dims = c(2, 2, 4, 1),
    microclim_dims = dim(microclimate)
  )

  expect_equal(dim(result)[3], 4)
  expect_equal(
    result[,,3,],
    microclimate[,,2,]
  )
  expect_equal(
    result[,,4,],
    microclimate[,,2,]
  )
})

test_that("microclimate height is trimmed when too tall", {
  microclimate <- array(runif(16), dim = c(2, 2, 4, 1))

  result <- adapt_microclimate_height(
    microclimate,
    microhab_dims = c(2, 2, 2, 1),
    microclim_dims = dim(microclimate)
  )

  expect_equal(dim(result)[3], 2)
  expect_equal(
    result,
    microclimate[,,1:2,, drop = FALSE]
  )
})

test_that("single-layer microclimate repeats correctly", {
  microclimate <- array(42, dim = c(2, 2, 1, 1))

  result <- adapt_microclimate_height(
    microclimate,
    microhab_dims = c(2, 2, 3, 1),
    microclim_dims = dim(microclimate)
  )

  expect_equal(dim(result)[3], 3)
  expect_true(all(result == 42))
})

test_that("multiple variables are preserved during height adaptation", {
  microclimate <- array(
    runif(16),
    dim = c(2, 2, 1, 2)
  )

  result <- adapt_microclimate_height(
    microclimate,
    microhab_dims = c(2, 2, 2, 2),
    microclim_dims = dim(microclimate)
  )

  expect_equal(dim(result), c(2, 2, 2, 2))
  expect_equal(result[,,2,], microclimate[,,1,])
})
