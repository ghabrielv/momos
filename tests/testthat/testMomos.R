context("Momos")

test_that("Get params without params", {
  expect_type(get_params(), "double")
})

test_that("Get params with params", {
  params <- NA
  expect_type(get_params(params), "logical")
})

test_that("Calculate momos without params", {
  params <- NA
  expect_type(calculate_momos(params), "logical")
})

test_that("Calibrate momos without params", {
  params <- NA
  expect_type(calibrate_momos(params), "logical")
})

test_that("Momos within params", {
  params <- list()
  expect_type(momos(params), "list")
})

test_that("Momos with params", {
  params <- list(from = 0, to = 60)
  expect_type(momos(params), "list")
})

test_that("Momos with non-numerical parameters", {
  params <- list(from = "0")
  expect_type(momos(params), "logical")
})

test_that("Momos with non-numerical parameters", {
  params <- c(0, 30)
  expect_type(momos(params), "logical")
})

test_that("Momos with calibrate", {
  params <- list()
  expect_type(momos(params), "list")
})
