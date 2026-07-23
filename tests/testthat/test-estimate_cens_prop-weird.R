event_surv_valid <- function(t) exp(-0.2 * t)

cens_surv_valid <- function(t, lambda_c) {
  exp(-lambda_c * t)
}

cens_dens_valid <- function(t, lambda_c) {
  lambda_c * exp(-lambda_c * t)
}

test_that("estimate_cens_prop rejects missing required arguments", {
  expect_error(
    estimate_cens_prop(
      event_survival = NULL
    )
  )
})


test_that("estimate_cens_prop checks event survival function signature", {
  # Missing t argument
  expect_error(
    estimate_cens_prop(
      event_survival = function(x) exp(-x)
    )
  )

  # t not first argument
  expect_error(
    estimate_cens_prop(
      event_survival = function(lambda_c, t) exp(-t)
    )
  )
})


test_that("estimate_cens_prop checks censoring survival signature", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      cens_survival = function(x, lambda_c) exp(-lambda_c * x),
      cens_density = function(t, lambda_c) lambda_c * exp(-lambda_c * t)
    )
  )

  # Missing target parameter
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      cens_survival = function(t) exp(-t),
      cens_density = function(t, lambda_c) {
        lambda_c * exp(-lambda_c * t)
      }
    )
  )
})


test_that("estimate_cens_prop checks censoring density signature", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      cens_survival = function(t, lambda_c) exp(-lambda_c * t),
      cens_density = function(t) exp(-t)
    )
  )
})


test_that("estimate_cens_prop validates target parameter", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      cens_survival = function(t, lambda_c) exp(-lambda_c * t),
      cens_density = function(t, lambda_c) {
        lambda_c * exp(-lambda_c * t)
      },
      target = "wrong_name"
    )
  )
})


test_that("estimate_cens_prop validates target argument consistency", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      cens_survival = function(t, lambda) exp(-lambda * t),
      cens_density = function(t, lambda) lambda * exp(-lambda * t),
      target = "lambda_c"
    )
  )
})


test_that("estimate_cens_prop validates target bounds", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      cens_survival = function(t, lambda_c) exp(-lambda_c * t),
      cens_density = function(t, lambda_c) {
        lambda_c * exp(-lambda_c * t)
      },
      target_bounds = 1
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      cens_survival = function(t, lambda_c) exp(-lambda_c * t),
      cens_density = function(t, lambda_c) {
        lambda_c * exp(-lambda_c * t)
      },
      target_bounds = c(2, 1)
    )
  )
})


test_that("estimate_cens_prop validates time limits", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      t_min = -1
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      t_max = -1
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      t_min = 10,
      t_max = 1
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      time_admin_cens = -1
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      time_accrual = -1
    )
  )
})


test_that("estimate_cens_prop validates tolerance", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      abs_tol = 0
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      abs_tol = 1
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      abs_tol = -1e-5
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t) exp(-t),
      abs_tol = 0.5
    )
  )
})


test_that("estimate_cens_prop validates covariate arguments", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t, x) exp(-t),
      covariate_density = function(x) 1,
      covariate_bounds = NULL
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = function(t, x) exp(-t),
      covariate_density = NULL,
      covariate_bounds = list(x = c(0, 1))
    )
  )
})


test_that("estimate_cens_prop rejects invalid covariate density functions", {
  bad_functions <- list(
    non_function = 1,
    null_return = function(x) NULL,
    character_return = function(x) "invalid",
    negative = function(x) -1
  )

  for (f in bad_functions) {
    expect_error(
      estimate_cens_prop(
        event_survival = function(t, x) exp(-t),
        covariate_density = f,
        covariate_bounds = list(x = c(0, 1))
      )
    )
  }
})


test_that("estimate_cens_prop returns a callable object", {
  f <- estimate_cens_prop(
    event_survival = event_surv_valid
  )

  expect_true(is.function(f))
  expect_s3_class(f, "cens_prop_fun")
  expect_true(is.function(print.cens_prop_fun))

  expect_length(formals(f), 2)
  expect_named(formals(f), c("lambda_c", "tau"))
})


test_that("estimate_cens_prop rejects invalid event survival functions", {
  bad_functions <- list(
    non_function = 1,
    null_return = function(t) NULL,
    character_return = function(t) "x",
    list_return = function(t) list(t),
    wrong_length = function(t) 1,
    missing_values = function(t) rep(NA_real_, length(t)),
    nan_values = function(t) rep(NaN, length(t)),
    infinite_values = function(t) rep(Inf, length(t)),
    below_zero = function(t) rep(-0.1, length(t)),
    above_one = function(t) rep(1.1, length(t)),
    non_decreasing = function(t) 1 - t,
    wrong_start = function(t) 0.5 + 0.5 * exp(-t)
  )

  for (f in bad_functions) {
    expect_error(
      estimate_cens_prop(
        event_survival = f,
        cens_survival = cens_surv_valid,
        cens_density = cens_dens_valid
      )
    )
  }
})


test_that("estimate_cens_prop checks event survival vectorization", {
  expect_error(
    estimate_cens_prop(
      event_survival = function(t) {
        if (length(t) > 1) {
          stop("not vectorized")
        }
        exp(-t)
      }
    )
  )
})


test_that("estimate_cens_prop rejects invalid censor survival functions", {
  bad_functions <- list(
    non_function = 1,
    null_return = function(t, lambda_c) NULL,
    wrong_type = function(t, lambda_c) rep("x", length(t)),
    wrong_length = function(t, lambda_c) 1,
    missing_values = function(t, lambda_c) rep(NA_real_, length(t)),
    negative = function(t, lambda_c) rep(-1, length(t)),
    above_one = function(t, lambda_c) rep(2, length(t))
  )

  for (f in bad_functions) {
    expect_error(
      estimate_cens_prop(
        event_survival = event_surv_valid,
        cens_survival = f,
        cens_density = cens_dens_valid
      )
    )
  }
})


test_that("estimate_cens_prop rejects invalid censor density functions", {
  bad_functions <- list(
    non_function = 1,
    null_return = function(t, lambda_c) NULL,
    wrong_type = function(t, lambda_c) rep("x", length(t)),
    wrong_length = function(t, lambda_c) 1,
    missing_values = function(t, lambda_c) rep(NA_real_, length(t)),
    negative = function(t, lambda_c) rep(-1, length(t))
  )

  for (f in bad_functions) {
    expect_error(
      estimate_cens_prop(
        event_survival = event_surv_valid,
        cens_survival = cens_surv_valid,
        cens_density = f
      )
    )
  }
})


test_that("estimate_cens_prop checks censor functions vectorization", {
  expect_error(
    estimate_cens_prop(
      event_survival = event_surv_valid,
      cens_survival = function(t, lambda_c) {
        if (length(t) > 1) stop("not vectorized")
        exp(-lambda_c * t)
      }
    )
  )

  expect_error(
    estimate_cens_prop(
      event_survival = event_surv_valid,
      cens_density = function(t, lambda_c) {
        if (length(t) > 1) stop("not vectorized")
        lambda_c * exp(-lambda_c * t)
      }
    )
  )
})


test_that("returned function evaluates with target and tau", {
  f <- estimate_cens_prop(event_survival = event_surv_valid)

  expect_no_error(
    f(lambda_c = 0.2)
  )

  expect_no_error(
    f(lambda_c = 0.2, tau = 5)
  )
})
