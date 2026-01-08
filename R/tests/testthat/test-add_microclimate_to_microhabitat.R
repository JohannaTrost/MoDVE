library(testthat)
source("02_add_microclimate_dimensions.R", local = TRUE)

test_that("add_microclimate_to_microhabitat runs for one timestep with provided config", {

  # ---- temporary directories ----
  tmp_microhabitat <- file.path(tempdir(), "microhab")
  tmp_microclimate <- file.path(tempdir(), "microclim")
  tmp_output       <- file.path(tempdir(), "output")

  dir.create(tmp_microhabitat, showWarnings = FALSE)
  dir.create(tmp_microclimate, showWarnings = FALSE)
  dir.create(tmp_output, showWarnings = FALSE)

  # ---- synthetic microhabitat ----
  microhabitat <- array(
    1,
    dim = c(2, 2, 2, 1)  # x, y, height, variables
  )

  saveRDS(
    microhabitat,
    file.path(tmp_microhabitat, "microhabitatMatrix1.rds")
  )

  # ---- synthetic microclimate ----
  # Needs at least 11 variables because indices 1, 7, 11 are accessed
  microclimate <- array(
    2,
    dim = c(2, 2, 1, 11)
  )

  saveRDS(
    microclimate,
    file.path(tmp_microclimate, "2000_test_mc_matrix.rds")
  )

  # ---- required forest files ----
  file.create(file.path(tmp_microhabitat, "Forest_param_global.txt"))
  file.create(file.path(tmp_microhabitat, "Forest_param_pass0.txt"))
  saveRDS(1, file.path(tmp_microhabitat, "dimPlot.rds"))

  # ---- config passed directly to function ----
  config <- list(
    timeStepStart = 1,
    timeStepEnd   = 1,
    yearStart     = 2000,
    region        = "test",

    Directorymicrohabitat     = tmp_microhabitat,
    Directorymicroclimate     = tmp_microclimate,
    Directorynew_microhabitat = tmp_output,

    TotalSurfaceAreaOpt   = 1,
    SurfaceAreaLossOpt    = 0,
    LightNicheOpt         = 0,
    AverageWeightedAngles = 0,
    HumNicheOpt           = 1,
    TempNicheOpt          = 0,
    WindNicheOpt          = 0
  )

  # ---- run function ----
  result <- add_microclimate_to_microhabitat(config)

  # ---- expectations ----
  expect_true(isTRUE(result))

  # output file exists
  out_file <- file.path(tmp_output, "microhabitatMatrix1.rds")
  expect_true(file.exists(out_file))

  # output content sanity checks
  out <- readRDS(out_file)

  # height adapted from 1 -> 2
  expect_equal(dim(out)[3], 2)

  # original microhabitat vars (1) + microclimate vars (1)
  expect_equal(dim(out)[4], 2)
})

test_that("function fails if microhabitat file is missing", {
  tmp_microhabitat <- file.path(tempdir(), "missing_microhab")
  tmp_microclimate <- file.path(tempdir(), "microclim")
  tmp_output       <- file.path(tempdir(), "out_missing")

  dir.create(tmp_microclimate, showWarnings = FALSE)
  dir.create(tmp_output, showWarnings = FALSE)

  config <- list(
    timeStepStart = 1,
    timeStepEnd   = 1,
    yearStart     = 2000,
    region        = "test",

    Directorymicrohabitat     = tmp_microhabitat, # directory does NOT exist
    Directorymicroclimate     = tmp_microclimate,
    Directorynew_microhabitat = tmp_output,

    TotalSurfaceAreaOpt   = 1,
    SurfaceAreaLossOpt    = 0,
    LightNicheOpt         = 0,
    AverageWeightedAngles = 0,
    HumNicheOpt           = 1,
    TempNicheOpt          = 0,
    WindNicheOpt          = 0
  )

  expect_error(
    add_microclimate_to_microhabitat(config),
    regexp = "Microhabitat file not found"
  )
})


test_that("function fails if microclimate file is missing", {
  tmp_microhabitat <- file.path(tempdir(), "microhab_ok")
  tmp_microclimate <- file.path(tempdir(), "missing_microclim")
  tmp_output       <- file.path(tempdir(), "out_missing2")

  dir.create(tmp_microhabitat, showWarnings = FALSE)
  dir.create(tmp_output, showWarnings = FALSE)

  # create a minimal microhabitat
  microhab <- array(1, dim = c(2,2,2,1))
  saveRDS(microhab, file.path(tmp_microhabitat, "microhabitatMatrix1.rds"))

  config <- list(
    timeStepStart = 1,
    timeStepEnd   = 1,
    yearStart     = 2000,
    region        = "test",

    Directorymicrohabitat     = tmp_microhabitat,
    Directorymicroclimate     = tmp_microclimate, # missing
    Directorynew_microhabitat = tmp_output,

    TotalSurfaceAreaOpt   = 1,
    SurfaceAreaLossOpt    = 0,
    LightNicheOpt         = 0,
    AverageWeightedAngles = 0,
    HumNicheOpt           = 1,
    TempNicheOpt          = 0,
    WindNicheOpt          = 0
  )

  expect_error(
    add_microclimate_to_microhabitat(config),
    regexp = "Microclimate file not found"
  )
})


test_that("function fails if microclimate and microhabitat have incompatible dimensions", {
  tmp_microhabitat <- file.path(tempdir(), "microhab_dim")
  tmp_microclimate <- file.path(tempdir(), "microclim_dim")
  tmp_output       <- file.path(tempdir(), "out_dim")

  dir.create(tmp_microhabitat, showWarnings = FALSE)
  dir.create(tmp_microclimate, showWarnings = FALSE)
  dir.create(tmp_output, showWarnings = FALSE)

  # microhabitat is 2x2x2
  microhab <- array(1, dim = c(2,2,2,1))
  saveRDS(microhab, file.path(tmp_microhabitat, "microhabitatMatrix1.rds"))

  # microclimate is incompatible (3x2x1)
  microclim <- array(2, dim = c(3,2,1,1))
  saveRDS(microclim, file.path(tmp_microclimate, "2000_test_mc_matrix.rds"))

  config <- list(
    timeStepStart = 1,
    timeStepEnd   = 1,
    yearStart     = 2000,
    region        = "test",

    Directorymicrohabitat     = tmp_microhabitat,
    Directorymicroclimate     = tmp_microclimate,
    Directorynew_microhabitat = tmp_output,

    TotalSurfaceAreaOpt   = 1,
    HumNicheOpt           = 1,
    SurfaceAreaLossOpt    = 0,
    LightNicheOpt         = 0,
    AverageWeightedAngles = 0,
    TempNicheOpt          = 0,
    WindNicheOpt          = 0
  )

  expect_error(
    add_microclimate_to_microhabitat(config),
    regexp = "Microclimate and microhabitat dimensions do not match in x and y."
  )
})

test_that("function works with zero microclimate variables selected", {
  tmp_microhabitat <- file.path(tempdir(), "microhab_no_clim")
  tmp_microclimate <- file.path(tempdir(), "microclim_no_clim")
  tmp_output       <- file.path(tempdir(), "out_no_clim")

  dir.create(tmp_microhabitat, showWarnings = FALSE)
  dir.create(tmp_microclimate, showWarnings = FALSE)
  dir.create(tmp_output, showWarnings = FALSE)

  microhab <- array(1, dim = c(2,2,2,1))
  saveRDS(microhab, file.path(tmp_microhabitat, "microhabitatMatrix1.rds"))

  microclim <- array(2, dim = c(2,2,3,11))
  saveRDS(microclim, file.path(tmp_microclimate, "2000_test_mc_matrix.rds"))

  # ---- required forest files ----
  file.create(file.path(tmp_microhabitat, "Forest_param_global.txt"))
  file.create(file.path(tmp_microhabitat, "Forest_param_pass0.txt"))
  saveRDS(1, file.path(tmp_microhabitat, "dimPlot.rds"))

  config <- list(
    timeStepStart = 1,
    timeStepEnd   = 1,
    yearStart     = 2000,
    region        = "test",

    Directorymicrohabitat     = tmp_microhabitat,
    Directorymicroclimate     = tmp_microclimate,
    Directorynew_microhabitat = tmp_output,

    TotalSurfaceAreaOpt   = 1,
    SurfaceAreaLossOpt    = 0,
    LightNicheOpt         = 0,
    AverageWeightedAngles = 0,
    HumNicheOpt           = 0,
    TempNicheOpt          = 0,
    WindNicheOpt          = 0
  )

  result <- add_microclimate_to_microhabitat(config)
  expect_true(isTRUE(result))

  out_file <- file.path(tmp_output, "microhabitatMatrix1.rds")
  expect_true(file.exists(out_file))

  out <- readRDS(out_file)
  # should have only original microhabitat variable
  expect_equal(dim(out)[4], 1)
})
