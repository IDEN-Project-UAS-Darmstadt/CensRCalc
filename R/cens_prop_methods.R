#' Print and plot methods for censoring proportion functions
#'
#' @param x A `"cens_prop_fun"` object.
#' @param ... Additional arguments.
#' @return The input object, invisibly.
#' @rdname cens_prop_methods
#' @aliases cens_prop_methods
#' @name S3 methods for cens_prop_fun
NULL

#' @rdname cens_prop_methods
#' @export
print.cens_prop_fun <- function(x, ...) {
  spec <- environment(x)$spec

  cat("<Expected censoring proportion function>\n")
  cat("Target parameter :", spec$target, "\n")
  cat("Target range     :", paste(spec$target_bounds, collapse = " - "), "\n")
  cat("Administrative censoring:", spec$time_admin_cens, "\n")
  cat("Accrual period   :", spec$time_accrual, "\n")

  invisible(x)
}

#' @rdname cens_prop_methods
#' @param n Number of evaluation points.
#' @param tau Optional time point for evaluating the expected censoring
#'  proportion. If `NULL`, the censoring proportion is evaluated overall,
#'  otherwise it is evaluated at the specified time point.
#' @export
plot.cens_prop_fun <- function(x, n = 100, tau = NULL, ...) {
  spec <- environment(x)$spec

  grid <- seq(
    spec$target_bounds[1],
    spec$target_bounds[2],
    length.out = n
  )

  values <- vapply(grid, x, numeric(1), tau = tau)

  ylab <- "Expected censoring proportion"
  if (!is.null(tau)) {
    ylab <- paste0(ylab, " at tau = ", tau)
  }

  plot(
    grid,
    values,
    type = "l",
    xlab = spec$target,
    ylab = ylab,
    ...
  )

  invisible(x)
}
