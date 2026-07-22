marginal_function <- function(fun, cov_bounds, cov_dens, abs_tol = 1e-6) {
  single_t <- function(t) {
    res <- integral_with_discr(
      \(...) fun(t = t, ...) * cov_dens(...),
      bounds = cov_bounds,
      abs_tol = abs_tol,
      vectorize = TRUE
    )
    logger::log_trace(
      "Marginal function evaluated at t = ", t, ": ",
      res$value, " (integral error estimate: ", res$error, ")"
    )
    res$value
  }
  return(Vectorize(single_t, vectorize.args = "t"))
}

check_surv_fun <- function(
  surv_fun, t_min = 0, t_max = Inf,
  n_points = 100L, abs_tol = 1e-6,
  abs_tol_buffer = 1e-2,
  cov_bounds = NULL,
  cov_dens = NULL
) {
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_function(surv_fun, add = coll)
  checkmate::assert_numeric(t_min,
    len = 1, any.missing = FALSE, finite = TRUE, lower = 0, add = coll
  )
  checkmate::assert_numeric(t_max,
    len = 1, any.missing = FALSE, finite = FALSE,
    lower = t_min, add = coll
  )
  checkmate::assert_int(n_points,
    na.ok = FALSE, lower = 2, add = coll
  )
  checkmate::assert_numeric(abs_tol,
    len = 1, any.missing = FALSE, finite = TRUE,
    lower = 0, add = coll
  )
  checkmate::assert_numeric(abs_tol_buffer,
    len = 1, any.missing = FALSE, finite = TRUE,
    lower = 0, add = coll
  )
  digits_to_show <- max(1, ceiling(-log10(abs_tol_buffer)))
  if (!coll$isEmpty()) {
    return(paste(coll$getMessages(), collapse = "; "))
  }

  m_surv_fun <- surv_fun
  # if we have covariates, we look at the marginal survival function instead
  if (!is.null(cov_bounds)) {
    m_surv_fun <- marginal_function(surv_fun, cov_bounds, cov_dens,
      abs_tol = abs_tol
    )
  }

  if (!is.finite(t_max)) {
    logger::log_trace(
      "t_max is infinite, trying to find a suitable t_max"
    )


    num_t_max <- locate_tail_bound(
      m_surv_fun,
      0,
      t_min = t_min,
      abs_tol = abs_tol
    )
    logger::log_trace("Setting t_max to ", num_t_max)
  } else {
    num_t_max <- t_max
  }
  t_values <- seq(t_min, num_t_max, length.out = n_points)
  surv_values <- m_surv_fun(t_values)
  if (length(surv_values) != n_points) {
    return(paste0(
      "Survival function failed vectorization check: got a length-",
      length(surv_values),
      " output, expected length ", n_points, "."
    ))
  }

  surv_values <- surv_values[!is.na(surv_values)]
  if (any(surv_values < -abs_tol_buffer) ||
    any(surv_values > 1 + abs_tol_buffer)) {
    return(paste0(
      "Survival function is out of bounds: values must fall within [0, 1] ",
      "(got range [", round(min(surv_values), digits_to_show), ", ",
      round(max(surv_values), digits_to_show), "])."
    ))
  }

  diffs <- diff(surv_values)
  if (any(diffs > abs_tol_buffer)) {
    return(paste0(
      "Survival function violates monotonicity: probabilities cannot increase",
      "over time (got a maximum increase of ", round(
        max(diffs),
        digits_to_show
      ),
      ")."
    ))
  }

  val_upper <- m_surv_fun(num_t_max)
  if (!is.finite(val_upper) || abs(val_upper) > abs_tol_buffer) {
    return(paste0(
      "Survival function tail does not converge to 0 near upper bound (got ",
      round(val_upper, digits_to_show), ", expected 0)."
    ))
  }

  val_t_min <- m_surv_fun(t_min)
  if (abs(val_t_min - 1) > abs_tol_buffer) {
    return(paste0(
      "Survival function baseline condition failed: S(t_min) must equal 1 ",
      "(got ", round(val_t_min, digits_to_show), ", expected 1)."
    ))
  }
  return(TRUE)
}
assert_surv_fun <- checkmate::makeAssertionFunction(check_surv_fun)


check_event_density <- function(
  dens_fun,
  t_min = 0,
  t_max = Inf,
  n_points = 100L,
  abs_tol = 1e-6,
  num_t_max = NA,
  abs_tol_buffer = 1e-2,
  cov_bounds = NULL,
  cov_dens = NULL
) {
  coll <- checkmate::makeAssertCollection()

  checkmate::assert_function(dens_fun, add = coll)
  checkmate::assert_numeric(t_min,
    len = 1, any.missing = FALSE, finite = TRUE,
    lower = 0, add = coll
  )
  checkmate::assert_numeric(t_max,
    len = 1, any.missing = FALSE, finite = FALSE,
    lower = t_min, add = coll
  )
  checkmate::assert_int(n_points, lower = 2, add = coll)
  checkmate::assert_numeric(abs_tol,
    len = 1, any.missing = FALSE,
    finite = TRUE, lower = 0, add = coll
  )
  checkmate::assert_numeric(num_t_max,
    len = 1, any.missing = TRUE,
    finite = TRUE, lower = t_min, add = coll
  )
  checkmate::assert_numeric(abs_tol_buffer,
    len = 1, any.missing = FALSE,
    finite = TRUE, lower = 0, add = coll
  )
  digits_to_show <- max(1, ceiling(-log10(abs_tol_buffer)))

  if (!coll$isEmpty()) {
    return(paste(coll$getMessages(), collapse = "; "))
  }

  # if we have covariates, we look at the marginal density function instead
  m_dens_fun <- dens_fun
  if (!is.null(cov_bounds)) {
    m_dens_fun <- marginal_function(dens_fun, cov_bounds, cov_dens,
      abs_tol = abs_tol
    )
  }


  # Resolve infinite time horizons safely
  if (!is.finite(t_max) && is.na(num_t_max)) {
    logger::log_trace(paste0(
      "t_max is infinite, searching for a suitable ",
      "truncation bound start position"
    ))

    # we start to go right to skip the first low part of the density
    t_start_search <- t_min + 1e-3
    found_support <- FALSE
    for (k in seq_len(20L)) {
      val <- m_dens_fun(t_start_search)
      checkmate::assert_numeric(
        val,
        len = 1, any.missing = FALSE, finite = TRUE,
        .var.name = "dens_fun(t) output during pre-search"
      )
      if (val > abs_tol) {
        found_support <- TRUE
        break
      }
      t_start_search <- t_start_search * 3
    }

    if (!found_support) {
      logger::log_warn(
        "Pre-search did not find dens_fun(t) > abs_tol within 20 iterations ",
        "(t reached ", t_start_search, "). Falling back to locate_tail_bound ",
        "with this horizon; result may be unreliable."
      )
    }

    num_t_max <- locate_tail_bound(
      m_dens_fun,
      sup_lim = 0,
      t_min = t_min,
      t_start = t_start_search * 2,
      abs_tol = abs_tol
    )
    logger::log_trace(
      "Setting numerical validation t_max horizon to ",
      num_t_max
    )
  } else if (is.finite(t_max)) {
    num_t_max <- t_max
  }

  # Generate evaluation grid
  t_values <- seq(t_min, num_t_max, length.out = n_points)
  dens_values <- m_dens_fun(t_values)

  if (length(dens_values) != n_points) {
    return(paste0(
      "Density function failed vectorization check: got a length-",
      length(dens_values),
      " output, expected length ", n_points, "."
    ))
  }
  if (any(dens_values < -abs_tol_buffer)) {
    return(paste0(
      "Density function violates probability constraints: negative density ",
      "values detected (got a minimum value of ", round(
        min(dens_values),
        digits_to_show
      ),
      ")."
    ))
  }
  if (any(!is.finite(dens_values))) {
    return(paste0(
      "Density function returns non-finite (NA/Inf) values within its active ",
      "support (", sum(!is.finite(dens_values)), " of ", n_points,
      " evaluated points)."
    ))
  }

  dens_integral <- integral_with_discr(
    m_dens_fun,
    bounds = list(t = c(t_min, t_max)),
    abs_tol = abs_tol,
    vectorize = TRUE
  )

  logger::log_debug(
    "Density total integral calculated as: ",
    dens_integral$value
  )

  if (abs(dens_integral$value - 1) > abs_tol_buffer) {
    return(paste0(
      "Density function structural failure: total area under curve does not ",
      "equal 1 (got ", round(dens_integral$value, digits_to_show),
      ", expected 1)."
    ))
  }

  TRUE
}

assert_event_density <- checkmate::makeAssertionFunction(check_event_density)

check_surv_and_density_funs <- function(
  surv_fun, dens_fun, t_min = 0,
  t_max = Inf, n_points = 100L, abs_tol = 1e-6, num_t_max = NA,
  abs_tol_buffer = 1e-2,
  cov_bounds = NULL, cov_dens = NULL
) {
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_function(surv_fun, add = coll)
  checkmate::assert_function(dens_fun, add = coll)
  checkmate::assert_numeric(t_min,
    len = 1, any.missing = FALSE, finite = TRUE,
    lower = 0, add = coll
  )
  checkmate::assert_numeric(t_max,
    len = 1, any.missing = FALSE, finite = FALSE,
    lower = t_min, add = coll
  )
  checkmate::assert_int(
    n_points,
    lower = 2, add = coll
  )
  checkmate::assert_numeric(abs_tol,
    len = 1, any.missing = FALSE, finite = TRUE,
    lower = 0, add = coll
  )
  checkmate::assert_numeric(num_t_max,
    len = 1, any.missing = TRUE, finite = TRUE,
    lower = t_min, add = coll
  )
  checkmate::assert_numeric(abs_tol_buffer,
    len = 1, any.missing = FALSE, finite = TRUE,
    lower = 0, add = coll
  )
  digits_to_show <- max(1, ceiling(-log10(abs_tol_buffer)))
  if (!coll$isEmpty()) {
    return(paste(coll$getMessages(), collapse = "; "))
  }

  m_dens_fun <- dens_fun
  m_surv_fun <- surv_fun
  if (!is.null(cov_bounds)) {
    m_dens_fun <- marginal_function(dens_fun, cov_bounds, cov_dens,
      abs_tol = abs_tol
    )
    m_surv_fun <- marginal_function(surv_fun, cov_bounds, cov_dens,
      abs_tol = abs_tol
    )
  }

  if (!is.finite(t_max) && is.na(num_t_max)) {
    logger::log_trace("t_max is infinite, trying to find a suitable t_max")
    num_t_max <- locate_tail_bound(
      m_surv_fun,
      t_min = t_min, sup_lim = 0,
      abs_tol = abs_tol
    )
    logger::log_trace("Setting t_max to ", num_t_max)
  } else if (is.finite(t_max)) {
    num_t_max <- t_max
  }
  res_surv_fun <- check_surv_fun(m_surv_fun,
    t_min = t_min,
    t_max = num_t_max, n_points = n_points, abs_tol = abs_tol,
    abs_tol_buffer = abs_tol_buffer
  )
  if (res_surv_fun != TRUE) {
    return(res_surv_fun)
  }

  res_dens_fun <- check_event_density(m_dens_fun,
    t_min = t_min,
    t_max = t_max, n_points = n_points, abs_tol = abs_tol,
    num_t_max = num_t_max, abs_tol_buffer = abs_tol_buffer,
    cov_bounds = NULL, cov_dens = NULL
  )
  if (res_dens_fun != TRUE) {
    return(res_dens_fun)
  }


  t_values <- seq(t_min + abs_tol, num_t_max - abs_tol, length.out = n_points)
  surv_values <- m_surv_fun(t_values)
  for (i in seq_along(t_values)) {
    t <- t_values[i]
    surv_val <- surv_values[i]
    dens_integral <- integral_with_discr(m_dens_fun,
      bounds = list(t = c(t_min, t)), abs_tol = abs_tol, vectorize = TRUE
    )
    if (abs(surv_val - (1 - dens_integral$value)) > abs_tol_buffer) {
      return(paste0(
        "Survival and density functions do not match at t = ",
        round(t, digits_to_show),
        " (got survival value ", round(surv_val, digits_to_show), ", expected ",
        round(1 - dens_integral$value, digits_to_show),
        " from integrated density)."
      ))
    }
  }
  TRUE
}
assert_surv_and_density_funs <-
  checkmate::makeAssertionFunction(check_surv_and_density_funs)

check_covariate_funs <- function(
  event_survival, covariate_density,
  covariate_bounds
) {
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_function(event_survival, add = coll)
  checkmate::assert_function(covariate_density, add = coll)
  checkmate::assert_list(covariate_bounds,
    types = c("numeric", "list"),
    any.missing = FALSE, min.len = 1, names = "unique", add = coll
  )
  if (!coll$isEmpty()) {
    return(paste(coll$getMessages(), collapse = "; "))
  }
  # Check that the functions have the same arguments
  # (except t in event_survival)
  ev_surv_args <- partial_formalargs(event_survival)
  if (!"t" %in% ev_surv_args) {
    return("event_survival must have an argument `t`.")
  }
  ev_surv_args <- setdiff(ev_surv_args, "t")
  cov_dens_args <- partial_formalargs(covariate_density)
  if (!checkmate::test_set_equal(ev_surv_args, cov_dens_args)) {
    return(paste0(
      "event_survival and covariate_density must have the same arguments",
      "(except `t`)."
    ))
  }
  if (!checkmate::test_set_equal(ev_surv_args, names(covariate_bounds))) {
    return(paste0(
      "covariate_bounds must have the same names as the arguments of ",
      "event_survival (except `t`)."
    ))
  }
  TRUE
}
assert_covariate_funs <- checkmate::makeAssertionFunction(check_covariate_funs)

check_covariate_density_fun <- function(
  covariate_density, covariate_bounds,
  abs_tol = 1e-6,
  abs_tol_buffer = 1e-2
) {
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(abs_tol,
    len = 1, any.missing = FALSE, finite = TRUE,
    lower = 0, add = coll
  )
  checkmate::assert_numeric(abs_tol_buffer,
    len = 1, any.missing = FALSE, finite = TRUE,
    lower = 0, add = coll
  )
  digits_to_show <- max(1, ceiling(-log10(abs_tol_buffer)))
  if (!coll$isEmpty()) {
    return(paste(coll$getMessages(), collapse = "; "))
  }
  int_res <- integral_with_discr(covariate_density, covariate_bounds,
    abs_tol = abs_tol,
    vectorize = TRUE
  )
  if (abs(int_res$value - 1) > (abs_tol_buffer)) {
    return(paste0(
      "Covariate density function does not integrate to 1 (got ",
      round(int_res$value, digits_to_show), ", expected 1)."
    ))
  }
  return(TRUE)
}
assert_covariate_density_fun <-
  checkmate::makeAssertionFunction(check_covariate_density_fun)


locate_tail_bound <- function(
  fun,
  sup_lim,
  t_min = 0,
  t_start = t_min + abs_tol,
  factor = 2,
  max_iter = 100L,
  abs_tol = 1e-6
) {
  coll <- checkmate::makeAssertCollection()

  checkmate::assert_function(fun, add = coll)
  checkmate::assert_numeric(sup_lim,
    len = 1, any.missing = FALSE,
    finite = TRUE, add = coll
  )
  checkmate::assert_numeric(t_min,
    len = 1, any.missing = FALSE,
    finite = TRUE, lower = 0, add = coll
  )
  checkmate::assert_numeric(t_start,
    len = 1, any.missing = FALSE,
    finite = TRUE, lower = t_min, add = coll
  )
  checkmate::assert_numeric(factor,
    len = 1, any.missing = FALSE,
    finite = TRUE, lower = 1.001, add = coll
  )
  checkmate::assert_int(max_iter, lower = 1, add = coll)
  checkmate::assert_numeric(abs_tol,
    len = 1, any.missing = FALSE,
    finite = TRUE, lower = 0, add = coll
  )
  # this is not a check function, so we throw an error if the assertions fail
  checkmate::reportAssertions(coll)

  t <- t_start
  for (i in seq_len(max_iter)) {
    val <- fun(t)

    checkmate::assert_numeric(
      val,
      len = 1,
      any.missing = FALSE,
      finite = TRUE,
      .var.name = "fun(t) output"
    )

    if (abs(val - sup_lim) < abs_tol) {
      logger::log_trace("Found t_max = ", t, " after ", i, " iterations.")
      return(t)
    }

    t <- t * factor
  }

  logger::log_warn(
    "Tail search exceeded max_iter (", max_iter, "). ",
    "Function value at t = ", t, " is ", val, " (target: ", sup_lim, "). ",
    "Returning unconverged horizon."
  )

  return(t)
}

survival_fun_safety_wrap <- function(surv_fun, abs_tol_buffer = 1e-2) {
  function(...) {
    res <- surv_fun(...)
    if (!is.numeric(res)) {
      stop("Survival function returned invalid value.")
    }
    not_missing_res <- res[!is.na(res)]
    if (any(not_missing_res < -abs_tol_buffer) ||
      any(not_missing_res > 1 + abs_tol_buffer)) {
      stop("Survival function returned out-of-bounds value.")
    }
    res
  }
}

density_fun_safety_wrap <- function(dens_fun, abs_tol_buffer = 1e-2) {
  function(...) {
    res <- dens_fun(...)
    if (!is.numeric(res)) {
      stop("Density function returned invalid value.")
    }
    not_missing_res <- res[!is.na(res)]
    if (any(not_missing_res < -abs_tol_buffer)) {
      stop("Density function returned negative value.")
    }
    res
  }
}
