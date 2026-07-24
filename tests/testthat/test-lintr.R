test_that("Code is lint free", {
  skip_on_cran()
  if (requireNamespace("lintr", quietly = TRUE)) {
    lintr::expect_lint_free()
  }
})
