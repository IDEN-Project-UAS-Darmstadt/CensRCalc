test_that("integral_with_discr integrates a simple continuous function", {
  res <- integral_with_discr(
    fun = function(x) x,
    bounds = list(x = c(0, 1))
  )

  expect_type(res, "list")
  expect_named(res, c("value", "error", "int_objects", "errors"))

  expect_equal(res$value, 0.5, tolerance = 1e-6)
  expect_true(res$error >= 0)
})


test_that("integral_with_discr handles discrete variables", {
  res <- integral_with_discr(
    fun = function(x) x,
    bounds = list(x = list(1, 2, 3))
  )
  expect_equal(res$value, 6)
})


test_that("integral_with_discr handles mixed continuous, discrete variables", {
  res <- integral_with_discr(
    fun = function(x, y) x * y,
    bounds = list(
      x = c(0, 1),
      y = list(1, 2, 3)
    )
  )

  # sum_y integral_0^1 x*y dx
  # = (1+2+3) * 1/2
  expect_equal(res$value, 3, tolerance = 1e-6)
})


test_that("integral_with_discr supports multiple continuous dimensions", {
  res <- integral_with_discr(
    fun = function(x, y) x + y,
    bounds = list(
      x = c(0, 1),
      y = c(0, 1)
    )
  )

  # Integral over unit square:
  # int int (x+y) dx dy = 1
  expect_equal(res$value, 1, tolerance = 1e-6)
})


test_that("integral_with_discr works with vectorized functions", {
  res <- integral_with_discr(
    fun = function(x) x^2,
    bounds = list(x = c(0, 1)),
    vectorize = TRUE
  )

  expect_equal(res$value, 1 / 3, tolerance = 1e-6)
})


test_that("integral_with_discr respects numerical tolerance", {
  res_low_tol <- integral_with_discr(
    fun = function(x) exp(x),
    bounds = list(x = c(0, 1)),
    abs_tol = 1e-4
  )

  res_high_tol <- integral_with_discr(
    fun = function(x) exp(x),
    bounds = list(x = c(0, 1)),
    abs_tol = 1e-8
  )

  expect_equal(
    res_low_tol$value,
    exp(1) - 1,
    tolerance = 1e-4
  )

  expect_equal(
    res_high_tol$value,
    exp(1) - 1,
    tolerance = 1e-8
  )
})


test_that("integral_with_discr rejects invalid bounds", {
  expect_error(
    integral_with_discr(
      fun = function(x) x,
      bounds = list(x = c(0, 1, 2))
    )
  )

  expect_error(
    integral_with_discr(
      fun = function(x) x,
      bounds = list(x = "invalid")
    )
  )
})


test_that("integral_with_discr handles constant functions", {
  res <- integral_with_discr(
    fun = function(x) 2,
    bounds = list(x = c(0, 3))
  )

  expect_equal(res$value, 6, tolerance = 1e-6)
})


test_that("integral_with_discr returns non-negative errors", {
  res <- integral_with_discr(
    fun = function(x) sin(x),
    bounds = list(x = c(0, pi))
  )

  expect_true(all(res$errors >= 0))
  expect_true(res$error >= 0)
})
