#' Determine a censoring parameter for a target censoring proportion
#'
#' @description
#' \loadmathjax
#' Solves for the parameter of a random censoring distribution such that the
#' expected censoring proportion equals a predefined target proportion.
#' This is designed for simulation studies where target censoring proportions
#' are required. Administrative censoring and accrual periods are
#' supported.
#'
#' @details
#' The parameter \mjeqn{\lambda_{C_{\text{rnd}}}}{lambda_C_rnd} is determined by
#' solving the root finding problem
#'
#' \mjdeqn{P(\delta=0\mid\lambda_{C_{\text{rnd}}}) - p_C = 0}{
#' P(delta=0|lambda_C_rnd)- p_C = 0}
#'
#' where \mjeqn{P(\delta=0\mid\lambda_{C_{\text{rnd}}})}{
#' P(delta=0|lambda_C_rnd)}
#' is the expected censoring proportion over the covariate distributions and
#' \mjeqn{p_C\in[0,1]}{p_C in [0,1]} is the target censoring proportion.
#' The expected censoring proportion \mjeqn{P(\delta=0\mid
#' \lambda_{C_{\text{rnd}}})=
#' E_X\left[P(\delta=0\mid \mathbf{X},\lambda_{C_{\text{rnd}}})\right]}{
#' P(delta=0|lambda_C_rnd)=E_X[P(delta=0|X,lambda_C_rnd)]} is
#' evaluated by [estimate_cens_prop()] and its mathematical background is
#' detailed in the package description ([CensRCalc]). A
#' bracketing root finder is used with the interval given by `target_bounds`.
#'
#' Convergence depends on the chosen bounds covering a sign change of the
#' objective function. The tolerance `abs_tol` is used for
#' both integration accuracy and the root-finding convergence criterion.
#'
#' Progress reporting is supported via the \pkg{progressr} package. When a
#' progressr handler is enabled (for example,
#' `progressr::with_progress()`), the iterations of the bracketing solver
#' emit progress updates.
#'
#' @param cens_prop Numeric scalar in \mjeqn{[0,1]}{[0,1]}. Expected target
#'   censoring proportion in the simulated data.
#' @inheritParams estimate_cens_prop
#' @param tau Numeric scalar. Optional time point at which the expected
#'  censoring proportion is evaluated. If `NULL`, the expected censoring
#'  proportion is evaluated overall, otherwise it is evaluated at the specified
#'  time point.
#' @param target_bounds Numeric length-2 vector. Lower and upper bounds of
#'   the search interval for the parameter named in `target`.
#'
#' @return A list with components:
#'   - `parameter`: the value of the censoring parameter.
#'   - `cens_prop`: the achieved expected censoring proportion at the
#'     calculated parameter.
#'
#' @seealso [estimate_cens_prop()] for the analytical evaluation of the
#'   expected censoring proportion, and [CensRCalc] for
#'   mathematical details.
#'
#' @import mathjaxr
#' @importFrom Rdpack reprompt
#'
#' @examples
#' # Target 30% censoring under exponential models without covariates
#' event_surv <- function(t) exp(-0.2 * t)
#'
#' sol <- find_cens_param(
#'   cens_prop = 0.3,
#'   event_survival = event_surv
#' )
#'
#' sol$parameter
#' @export
find_cens_param <- function(
  cens_prop,
  event_survival,
  covariate_density = NULL,
  covariate_bounds = NULL,
  cens_survival = \(t, lambda_c) {
    1 - stats::pexp(t, rate = lambda_c)
  },
  cens_density = \(t, lambda_c) {
    stats::dexp(t, rate = lambda_c)
  },
  tau = NULL,
  time_admin_cens = t_max,
  time_accrual = t_min,
  target = "lambda_c",
  target_bounds = c(1e-5, 5),
  t_min = 0,
  t_max = Inf,
  abs_tol = 1e-6
) {
  checkmate::assert_numeric(
    cens_prop,
    lower = 0, upper = 1,
    finite = TRUE,
    len = 1,
    any.missing = FALSE
  )
  cens_prop_fun_no_tau <- estimate_cens_prop(
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
  cens_prop_fun <- partial(cens_prop_fun_no_tau, list(tau = tau))
  result <- uniroot_with_progress(
    fun = cens_prop_fun,
    target = cens_prop,
    bounds = target_bounds,
    abs_tol = abs_tol
  )

  list(
    parameter = result$root,
    # because f.root = P(delta=0|lambda_C_rnd) - cens_prop
    cens_prop = result$f.root + cens_prop
  )
}
