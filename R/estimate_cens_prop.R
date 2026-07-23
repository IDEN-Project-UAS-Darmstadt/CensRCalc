#' Calculate expected censoring proportion
#'
#' @description
#' \loadmathjax
#' Evaluates the expected censoring proportion under
#' independent censoring by integrating survival and censoring functions
#' over time and covariates. Administrative censoring and accrual periods
#' are supported.
#'
#' @details
#' This function returns the expected censoring proportion over the
#' covariates
#' \mjeqn{P(\delta=0\mid\lambda_{C_{\text{rnd}}})=E_X\Bigg[P(\delta=0\mid
#' \mathbf{X},\lambda_{C_{\text{rnd}}})\Bigg]}{P(delta=0|lambda_C_rnd)=
#' E_X[P(delta=0|X,lambda_C_rnd)}
#' obtained by numerical integration of survival and censoring functions
#' over time and finally integrating over the covariate distributions. The
#' mathematical background and the full expression for
#' \mjeqn{P(\delta=0\mid\lambda_{C_{\text{rnd}}})}{P(delta=0|lambda_C_rnd)} are
#' provided in the package description, see [CensRCalc].
#'
#' Functions supplied for time must be vectorized in their first argument
#' (time). The survival functions and densities of the event time and
#' censoring time variables, respectively, must be consistent, i.e., they
#' correspond to the same distribution.
#'
#' @param event_survival Function. Event-time survival function
#'   \mjeqn{S_{T^*\mid\mathbf{X}}(t)}{S_T*|X(t)}. The first argument must
#'   be `t` and further arguments are considered covariates. Must be vectorized
#'   in `t`.
#' @param covariate_density Optional function. Joint density of the
#'   covariate vector over `covariate_bounds`. If provided, it must
#'   integrate to one over the specified domain. If `NULL`, no covariates
#'   are integrated.
#' @param covariate_bounds Named list or `NULL`. Bounds for covariates
#'   used both to integrate `covariate_density`. Use numeric length-2
#'   vectors for continuous bounds and lists of numeric values for
#'   discrete variables. Names must match the non-`t` arguments of
#'   `event_survival`. Must be `NULL` when `covariate_density` is `NULL`.
#' @param cens_survival Function. Censoring-time survival function
#'   \mjeqn{S_{C_{\text{rnd}}\mid\lambda_{C_{\text{rnd}}}}(t)}{
#'   S_C_rnd|lambda_C_rnd(t)}
#'   with first argument `t` and the second argument named as in `target`.
#'   Must be vectorized in `t`.
#' @param cens_density Function. Censoring-time density
#'   \mjeqn{f_{C_{\text{rnd}}\mid\lambda_{C_{\text{rnd}}}}(t)}{
#'   f_C_rnd|lambda_C_rnd(t)}
#'   with first argument `t` and the second argument named as in `target`.
#'   Must be vectorized in `t` and be consistent with `cens_survival`.
#' @param time_admin_cens Numeric scalar. Administrative censoring time
#'   \mjeqn{\tau_{\text{adm}}}{tau_adm}. Use `Inf` for no administrative
#'   censoring.
#' @param time_accrual Numeric scalar. Accrual time span
#'   \mjeqn{\tau_{\text{acc}}}{tau_acc}. Use `0` for no accrual.
#' @param target Character scalar. Name of the censoring parameter inside
#'   `cens_survival` and `cens_density` (e.g., `"lambda_c"`).
#' @param target_bounds Numeric length-2 vector. Lower and upper bounds of
#'   the validation interval for the parameter named in `target`.
#' @param t_min Numeric scalar. Left boundary of the time support used for
#'   validation and numerical checks. Must be non-negative.
#' @param t_max Numeric scalar. Right boundary of the time support used
#'   for validation and numerical checks. Use `Inf` for an unbounded right
#'   tail.
#' @param abs_tol Numeric scalar. Absolute tolerance used for integration
#'   (on the probability scale). Must be positive and not too large
#'   (e.g., < 0.1).
#'
#' @return
#' A callable function of two numeric scalars, one named by `target` and
#' the other `tau`, which returns the expected censoring proportion under the
#' specified model components. When `tau` is specified, the expected
#' censoring proportion is evaluated at the time point `tau` instead of
#' overall time.
#'
#' The returned function has class `"cens_prop_fun"` and supports
#' S3 methods such as [print()] and [plot()] for summarizing and
#' visualizing the censoring model.
#'
#' @import mathjaxr
#' @importFrom Rdpack reprompt
#'
#' @examples
#' # Exponential event time without covariates
#' event_surv <- function(t) exp(-0.2 * t)
#'
#' cens_surv <- function(t, lambda_c) exp(-lambda_c * t)
#' cens_dens <- function(t, lambda_c) lambda_c * exp(-lambda_c * t)
#'
#' f_obj <- estimate_cens_prop(
#'   event_survival = event_surv,
#'   covariate_density = NULL,
#'   covariate_bounds = NULL,
#'   cens_survival = cens_surv,
#'   cens_density = cens_dens,
#'   time_admin_cens = 10,
#'   time_accrual = 0,
#' )
#'
#' print(f_obj)
#' f_obj(0.3)
#' f_obj(0.3, tau = 5)
#' @export
estimate_cens_prop <- function(
  event_survival,
  covariate_density = NULL,
  covariate_bounds = NULL,
  cens_survival = function(t, lambda_c) {
    1 - stats::pexp(t, rate = lambda_c)
  },
  cens_density = function(t, lambda_c) {
    stats::dexp(t, rate = lambda_c)
  },
  time_admin_cens = t_max,
  time_accrual = t_min,
  target = "lambda_c",
  target_bounds = c(1e-5, 5),
  t_min = 0,
  t_max = Inf,
  abs_tol = 1e-6
) {
  # The functions need to be vectorized over t
  coll <- checkmate::makeAssertCollection()
  logger::log_trace("Validating arguments...")

  checkmate::assert_function(event_survival, add = coll)
  checkmate::assert_function(covariate_density, null.ok = TRUE, add = coll)
  checkmate::assert_list(
    covariate_bounds,
    types = c("numeric", "list"),
    any.missing = FALSE,
    null.ok = is.null(covariate_density),
    add = coll
  )
  checkmate::assert_character(
    target,
    len = 1,
    any.missing = FALSE,
    add = coll
  )
  checkmate::assert_function(
    cens_survival,
    args = c("t", target),
    nargs = 2,
    add = coll
  )
  checkmate::assert_function(
    cens_density,
    args = c("t", target),
    nargs = 2,
    add = coll
  )
  checkmate::assert_numeric(
    time_accrual,
    lower = t_min, upper = time_admin_cens,
    finite = TRUE,
    len = 1,
    any.missing = FALSE,
    add = coll
  )
  checkmate::assert_numeric(
    time_admin_cens,
    lower = time_accrual, upper = t_max,
    finite = FALSE,
    len = 1,
    any.missing = FALSE,
    add = coll
  )
  checkmate::assert_numeric(
    target_bounds,
    finite = TRUE,
    len = 2,
    any.missing = FALSE,
    sorted = TRUE,
    add = coll
  )
  checkmate::assert_numeric(
    t_min,
    len = 1,
    any.missing = FALSE,
    finite = TRUE,
    lower = 0,
    add = coll
  )
  checkmate::assert_numeric(
    t_max,
    len = 1,
    any.missing = FALSE,
    finite = FALSE,
    lower = t_min,
    add = coll
  )
  checkmate::assert_numeric(
    abs_tol,
    len = 1,
    any.missing = FALSE,
    finite = TRUE,
    lower = 0,
    upper = 0.1,
    add = coll
  )
  checkmate::reportAssertions(coll)
  coll <- checkmate::makeAssertCollection()
  logger::log_trace("Validating functions...")
  assert_surv_fun(
    event_survival,
    t_min = t_min,
    t_max = t_max,
    add = coll,
    abs_tol = abs_tol,
    cov_bounds = covariate_bounds,
    cov_dens = covariate_density
  )
  if (!is.null(covariate_density)) {
    assert_covariate_density_fun(
      covariate_density,
      covariate_bounds,
      abs_tol = abs_tol,
      add = coll
    )
    assert_covariate_funs(
      event_survival,
      covariate_density,
      covariate_bounds,
      add = coll
    )
  }
  checkmate::assert_function(event_survival,
    args = c("t", names(covariate_bounds)),
    nargs = 1 + length(covariate_bounds),
    add = coll
  )
  param_span <- target_bounds[2] - target_bounds[1]
  param_margin <- max(param_span * 0.01, 1e-4)

  check_target_vals <- seq(
    target_bounds[1] + param_margin,
    target_bounds[2] - param_margin,
    length.out = 5
  )

  for (check_target_val in check_target_vals) {
    args <- stats::setNames(list(check_target_val), target)
    partial_cens_survival <- partial(cens_survival, args)
    partial_cens_density <- partial(cens_density, args)
    add_name <- paste0("(", target, "=", check_target_val, ")")

    assert_surv_fun(
      partial_cens_survival,
      t_min = t_min,
      t_max = t_max,
      add = coll,
      .var.name = paste0("cens_survival", add_name),
      abs_tol = abs_tol
    )
    assert_event_density(
      partial_cens_density,
      t_min = t_min,
      t_max = t_max,
      add = coll,
      .var.name = paste0("cens_density", add_name),
      abs_tol = abs_tol
    )
    assert_surv_and_density_funs(
      partial_cens_survival,
      dens_fun = partial_cens_density,
      t_min = t_min,
      t_max = t_max,
      add = coll,
      .var.name = paste0(
        "cens_survival and cens_density ",
        add_name
      ),
      abs_tol = abs_tol
    )
  }
  checkmate::reportAssertions(coll)

  t_start_admin <- time_admin_cens - time_accrual

  surv_admin <- \(t) {
    if (is.infinite(time_admin_cens)) {
      rep_len(1, length(t))
    } else {
      1 - stats::punif(
        t,
        min = t_start_admin,
        max = time_admin_cens
      )
    }
  }

  dens_admin <- \(t) {
    stats::dunif(
      t,
      min = t_start_admin,
      max = time_admin_cens
    )
  }

  integrate_random_cens <- function(event_survival,
                                    cens_density_only_t,
                                    t_upper_rand,
                                    cov_dens = \(...) 1,
                                    cov_bounds = list(),
                                    abs_tol = 1e-6) {
    bounds <- list(
      t = c(t_min, t_upper_rand)
    )
    bounds <- c(bounds, cov_bounds)
    int_res <- integral_with_discr(
      \(t, ...) {
        event_survival(t, ...) *
          cens_density_only_t(t) *
          surv_admin(t) *
          cov_dens(...)
      },
      bounds = bounds,
      abs_tol = abs_tol,
      vectorize = TRUE
    )
    logger::log_debug(
      "Integral during random censoring(bounds: ",
      paste(sapply(bounds, function(x) paste(x, collapse = ", ")),
        collapse = "; "
      ),
      "): ",
      int_res$value,
      " with error ",
      int_res$error
    )
    list(value = int_res$value, error = int_res$error)
  }

  integrate_admin_cens <- function(event_survival,
                                   cens_survival_only_t,
                                   t_upper_rand,
                                   cov_dens = \(...) 1,
                                   cov_bounds = list(),
                                   abs_tol = 1e-6) {
    # Admin censoring only has an impact from:
    if (is.infinite(time_admin_cens)) {
      # With no admin censoring, no one will be admin censored
      list(value = 0, error = 0)
    } else if (t_start_admin >= t_upper_rand) {
      # Returns 0 if tau is before the admin censoring window,
      # or if time_accrual == 0 (where t_start_admin == time_admin_cens)
      list(value = 0, error = 0)
    } else {
      bounds <- list(
        t = c(t_start_admin, t_upper_rand)
      )
      bounds <- c(bounds, cov_bounds)

      int_res <- integral_with_discr(
        \(t, ...) {
          event_survival(t, ...) *
            cens_survival_only_t(t) *
            dens_admin(t) *
            cov_dens(...)
        },
        bounds = bounds,
        abs_tol = abs_tol,
        vectorize = TRUE
      )

      logger::log_debug(
        "Integral during administrative censoring(bounds: ",
        paste(sapply(bounds, function(x) paste(x, collapse = ", ")),
          collapse = "; "
        ),
        "): ",
        int_res$value,
        " with error ",
        int_res$error
      )
      list(value = int_res$value, error = int_res$error)
    }
  }


  combine_censoring_mechanisms <- function(event_survival,
                                           cens_survival_only_t,
                                           cens_density_only_t,
                                           t_upper_rand,
                                           cov_dens = \(...) 1,
                                           cov_bounds = list()) {
    event_survival_w <- survival_fun_safety_wrap(event_survival)
    cens_survival_only_t_w <- survival_fun_safety_wrap(cens_survival_only_t)
    cens_density_only_t_w <- density_fun_safety_wrap(cens_density_only_t)
    cov_dens_w <- density_fun_safety_wrap(cov_dens)
    first_tol <- abs_tol / 2
    irc <- integrate_random_cens(
      event_survival_w, cens_density_only_t_w, t_upper_rand, cov_dens_w,
      cov_bounds,
      abs_tol = first_tol
    )
    remaining_tol <- max(.Machine$double.eps, abs_tol - irc$error)
    iac <- integrate_admin_cens(event_survival_w, cens_survival_only_t_w,
      t_upper_rand, cov_dens_w, cov_bounds,
      abs_tol = remaining_tol
    )
    value <- irc$value + iac$value
    error <- irc$error + iac$error
    if (error > abs_tol) {
      logger::log_warn(
        "Estimated error ({error}) exceeds requested tolerance ({abs_tol})."
      )
    }
    logger::log_debug(
      "Total censoring proportion: ",
      irc$value + iac$value,
      " with error ", irc$error + iac$error
    )
    irc$value + iac$value
  }

  evaluate_cens_prop <- function(target_val, tau = NULL) {
    checkmate::assert_numeric(
      target_val,
      len = 1,
      any.missing = FALSE,
      finite = TRUE,
      lower = target_bounds[1],
      upper = target_bounds[2]
    )
    checkmate::assert_numeric(
      tau,
      len = 1,
      any.missing = FALSE,
      finite = TRUE,
      lower = t_min,
      upper = t_max,
      null.ok = TRUE,
      add = coll
    )
    if (is.null(tau)) {
      # Random censoring only happens up to:
      t_upper_rand <- time_admin_cens
    } else {
      t_upper_rand <- min(tau, time_admin_cens)
    }
    args <- stats::setNames(list(target_val), target)
    cens_density_only_t <- partial(cens_density, args)
    cens_survival_only_t <- partial(cens_survival, args)
    if (is.null(covariate_density)) {
      cov_dens <- \(...) 1
      cov_bounds <- list()
    } else {
      cov_dens <- covariate_density
      cov_bounds <- covariate_bounds
    }
    est_cens_prop <- combine_censoring_mechanisms(
      event_survival,
      cens_survival_only_t,
      cens_density_only_t,
      t_upper_rand,
      cov_dens,
      cov_bounds
    )

    logger::log_debug(
      "Estimated censoring proportion at ",
      target,
      " = ",
      target_val,
      ": ",
      est_cens_prop,
      " (tau = ",
      ifelse(is.null(tau), "NULL", tau),
      ")"
    )
    est_cens_prop
  }
  obj <- evaluate_cens_prop

  # rename the first argument of the returned function to the target name
  old_name <- "target_val"
  f_args <- formals(obj)
  names(f_args)[names(f_args) == old_name] <- target
  formals(obj) <- f_args
  body(obj) <- do.call(substitute, list(
    body(obj),
    stats::setNames(list(as.symbol(target)), old_name)
  ))

  environment(obj)$spec <- list(
    event_survival = event_survival,
    covariate_density = covariate_density,
    covariate_bounds = covariate_bounds,
    cens_survival = cens_survival,
    cens_density = cens_density,
    time_admin_cens = time_admin_cens,
    time_accrual = time_accrual,
    target = target,
    target_bounds = target_bounds,
    t_min = t_min,
    t_max = t_max,
    abs_tol = abs_tol
  )

  class(obj) <- c("cens_prop_fun", class(obj))

  obj
}
