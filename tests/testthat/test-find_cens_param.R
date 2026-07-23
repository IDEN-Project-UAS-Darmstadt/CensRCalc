test_that("The analytical example in the vignette returns the right answer", {
  tau_acc <- 1
  tau_adm <- 3
  p_c <- 0.3

  event_surv <- function(t, x) 1 - punif(t, min = 0, max = x)
  cens_surv <- function(t, lambda_c) 1 - punif(t, min = 0, max = lambda_c)
  cens_dens <- function(t, lambda_c) dunif(t, min = 0, max = lambda_c)
  cov_dens <- function(x) dunif(x, min = 1, max = 2)
  cov_bounds <- list(x = c(1, 2))
  bounds <- c(0.1, 10)

  result <- find_cens_param(p_c, event_surv, cov_dens, cov_bounds,
    time_accrual = tau_acc,
    time_admin_cens = tau_adm,
    cens_surv = cens_surv,
    cens_dens = cens_dens,
    target_bounds = bounds
  )
  expect_equal(result$parameter, 2.5, tolerance = 1e-5)
})
