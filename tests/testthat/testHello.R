context("Hello")

test_that("Hello", {
  expect_equal(hello("hello", "world"), "hello world")
})
