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
#' \mjeqn{P(\delta=0\mid\lambda_{C_{rnd}})=E_X\Bigg[P(\delta=0\mid
#' \mathbf{X},\lambda_{C_{rnd}})\Bigg]}{P(delta=0|lambda_C_rnd)=
#' E_X[P(delta=0|X,lambda_C_rnd)}
#' obtained by numerical integration of survival and censoring functions
#' over time and finally integrating over the covariate distributions. The
#' mathematical background and the full expression for
#' \mjeqn{P(\delta=0\mid\lambda_{C_{rnd}})}{P(delta=0|lambda_C_rnd)} are
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
#'   \mjeqn{S_{C_{rnd}\mid\lambda_{C_{rnd}}}(t)}{S_C_rnd|lambda_C_rnd(t)}
#'   with first argument `t` and the second argument named as in `target`.
#'   Must be vectorized in `t`.
#' @param cens_density Function. Censoring-time density
#'   \mjeqn{f_{C_{rnd}\mid\lambda_{C_{rnd}}}(t)}{f_C_rnd|lambda_C_rnd(t)}
#'   with first argument `t` and the second argument named as in `target`.
#'   Must be vectorized in `t` and be consistent with `cens_survival`.
#' @param time_admin_cens Numeric scalar. Administrative censoring time
#'   \mjeqn{\tau_{adm}}{tau_adm}. Use `Inf` for no administrative
#'   censoring.
#' @param time_accrual Numeric scalar. Accrual time span
#'   \mjeqn{\tau_{acc}}{tau_acc}. Use `0` for no accrual.
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
#' @return A function of one numeric scalar (named by `target`) that
#'   returns the expected censoring proportion in \mjeqn{[0,1]}{[0,1]}
#'   under the supplied model components.
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
#' f_obj(0.3)
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

  surv_admin <- \(t) {
    if (is.infinite(time_admin_cens)) {
      rep_len(1, length(t))
    } else {
      1 - stats::punif(
        t,
        min = time_admin_cens - time_accrual,
        max = time_admin_cens
      )
    }
  }

  dens_admin <- \(t) {
    if (is.infinite(time_admin_cens) || time_accrual == 0) {
      rep_len(0, length(t))
    } else {
      stats::dunif(
        t,
        min = time_admin_cens - time_accrual,
        max = time_admin_cens
      )
    }
  }

  integrate_random_cens <- function(event_survival,
                                    cens_density_only_t,
                                    cov_dens = \(...) 1,
                                    cov_bounds = list()) {
    bounds <- list(
      t = c(0, time_admin_cens)
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
      "Integral during random censoring: ",
      int_res$value,
      " with error ",
      int_res$error
    )
    int_res$value
  }

  integrate_admin_cens <- function(event_survival,
                                   cens_survival_only_t,
                                   cov_dens = \(...) 1,
                                   cov_bounds = list()) {
    if (is.infinite(time_admin_cens)) {
      # With no admin censoring, no one will be admin censored
      0
    } else {
      to_integrate <- \(t, ...) {
        event_survival(t, ...) *
          cens_survival_only_t(t) *
          dens_admin(t) *
          cov_dens(...)
      }
      if (time_accrual == 0 && length(cov_bounds) > 0) {
        # With no accrual, everyone is recruited at time 0
        int_res <- integral_with_discr(
          \(...) {
            to_integrate(time_admin_cens, ...)
          },
          bounds = cov_bounds,
          abs_tol = abs_tol,
          vectorize = TRUE
        )
      } else if (time_accrual == 0 && length(cov_bounds) == 0) {
        # With no accrual and no covariates, everyone is recruited at time 0
        int_res <- list(
          value = event_survival(time_admin_cens) *
            cens_survival_only_t(time_admin_cens) *
            dens_admin(time_admin_cens),
          error = 0
        )
      } else {
        bounds <- list(
          t = c(time_admin_cens - time_accrual, time_admin_cens)
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
      }
      logger::log_debug(
        "Integral during administrative censoring: ",
        int_res$value,
        " with error ",
        int_res$error
      )
      int_res$value
    }
  }


  combine_censoring_mechanisms <- function(event_survival,
                                           cens_survival_only_t,
                                           cens_density_only_t,
                                           cov_dens = \(...) 1,
                                           cov_bounds = list()) {
    event_survival_w <- survival_fun_safety_wrap(event_survival)
    cens_survival_only_t_w <- survival_fun_safety_wrap(cens_survival_only_t)
    cens_density_only_t_w <- density_fun_safety_wrap(cens_density_only_t)
    cov_dens_w <- density_fun_safety_wrap(cov_dens)
    integrate_random_cens(
      event_survival_w, cens_density_only_t_w, cov_dens_w, cov_bounds
    ) +
      integrate_admin_cens(
        event_survival_w, cens_survival_only_t_w, cov_dens_w,
        cov_bounds
      )
  }

  evaluate_cens_prop <- function(target_val) {
    checkmate::assert_numeric(
      target_val,
      len = 1,
      any.missing = FALSE,
      finite = TRUE,
      lower = target_bounds[1],
      upper = target_bounds[2]
    )
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
      cov_dens,
      cov_bounds
    )

    logger::log_debug(
      "Estimated censoring proportion at ",
      target,
      " = ",
      target_val,
      ": ",
      est_cens_prop
    )
    est_cens_prop
  }
  evaluate_cens_prop
}
