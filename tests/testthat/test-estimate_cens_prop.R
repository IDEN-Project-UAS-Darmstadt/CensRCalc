test_that("Symmetric U distributions yield exactly 0.5 censoring proportion", {
  event_surv <- \(t) pmax(0, 1 - t / 10)
  # pmax ensures it stays >= 0 for vectorized t

  cens_surv <- \(t, lambda_c) pmax(0, 1 - t / lambda_c)

  cens_dens <- \(t, lambda_c) ifelse(t >= 0 & t <= lambda_c, 1 / lambda_c, 0)


  f_obj <- estimate_cens_prop(
    event_survival = event_surv,
    cens_survival = cens_surv,
    cens_density = cens_dens,
    time_admin_cens = Inf,
    time_accrual = 0,
    target = "lambda_c",
    target_bounds = c(1, 20)
  )

  expect_equal(f_obj(10), 0.5, tolerance = 1e-5)
})

test_that("Admin censoring combined with random censoring (hard cutoff)", {
  event_surv <- \(t) pmax(0, 1 - t / 10)
  cens_surv <- \(t, lambda_c) pmax(0, 1 - t / lambda_c)
  cens_dens <- \(t, lambda_c) ifelse(t >= 0 & t <= lambda_c, 1 / lambda_c, 0)

  f_obj <- estimate_cens_prop(
    event_survival = event_surv,
    cens_survival = cens_surv,
    cens_density = cens_dens,
    time_admin_cens = 5,
    time_accrual = 0,
    t_max = Inf,
    target = "lambda_c",
    target_bounds = c(5, 15)
  )
  # Math:
  # Random censored: integral_0^5 (1/10) * (1 - t/10) dt = 0.375
  # Admin censored: S_T(5) * S_C(5) = 0.5 * 0.5 = 0.25
  # Total is = 0.625
  expect_equal(1 - f_obj(10), 0.625, tolerance = 1e-4)
})


test_that("Admin censoring with uniform accrual (3-way symmetry)", {
  event_surv <- \(t) pmax(0, 1 - t / 10)
  cens_surv <- \(t, lambda_c) pmax(0, 1 - t / lambda_c)
  cens_dens <- \(t, lambda_c) ifelse(t >= 0 & t <= lambda_c, 1 / lambda_c, 0)

  f_obj <- estimate_cens_prop(
    event_survival = event_surv,
    cens_survival = cens_surv,
    cens_density = cens_dens,
    time_admin_cens = 10,
    time_accrual = 10,
    t_max = Inf,
    target = "lambda_c",
    target_bounds = c(5, 15)
  )

  # T, C_rnd, and C_adm are all U(0, 10).
  # Probability that T is smallest is exactly 1/3.
  # Censoring proportion = 2/3.
  expect_equal(f_obj(10), 2 / 3, tolerance = 1e-4)
})

test_that("Print function works (default params)", {
  event_surv <- \(t) 1 - pexp(t, rate = 0.1)
  estimator <- estimate_cens_prop(event_surv)
  expect_output(print(estimator), "lambda_c")
  expect_output(print(estimator), "Administrative censoring")
})

test_that("Plot function creates a plot (default params)", {
  event_surv <- \(t) 1 - pexp(t, rate = 0.1)
  estimator <- estimate_cens_prop(event_surv)
  tmp_file <- tempfile(fileext = ".png")
  png(tmp_file)
  expect_silent(plot(estimator))
  dev.off()
  # is the file created with some content?
  expect_true(file.exists(tmp_file))
  expect_true(file.info(tmp_file)$size > 100)
})

test_that("Plot function creates a plot (tau)", {
  event_surv <- \(t) 1 - pexp(t, rate = 0.1)
  estimator <- estimate_cens_prop(event_surv)
  tmp_file <- tempfile(fileext = ".png")
  png(tmp_file)
  expect_silent(plot(estimator, tau = 5))
  dev.off()
  # is the file created with some content?
  expect_true(file.exists(tmp_file))
  expect_true(file.info(tmp_file)$size > 100)
})
