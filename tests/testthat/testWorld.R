context("World")

test_that("World", {
  expect_equal(world("hello", "world"), "hello world")
})
