context("Momos")

test_that("Momos within params", {
  expect_type(momos(), "list")
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
