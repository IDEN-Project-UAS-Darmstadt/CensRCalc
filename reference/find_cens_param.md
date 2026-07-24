# Determine a censoring parameter for a target censoring proportion

Solves for the parameter of a random censoring distribution such that
the expected censoring proportion equals a predefined target proportion.
This is designed for simulation studies where target censoring
proportions are required. Administrative censoring and accrual periods
are supported.

## Usage

``` r
find_cens_param(
  cens_prop,
  event_survival,
  covariate_density = NULL,
  covariate_bounds = NULL,
  cens_survival = function(t, lambda_c) {
1 - stats::pexp(t, rate = lambda_c)
 },
  cens_density = function(t, lambda_c) {
     stats::dexp(t, rate = lambda_c)
 },
  tau = NULL,
  time_admin_cens = t_max,
  time_accrual = t_min,
  target = "lambda_c",
  target_bounds = c(1e-05, 5),
  t_min = 0,
  t_max = Inf,
  abs_tol = 1e-06
)
```

## Arguments

- cens_prop:

  Numeric scalar in \\\[0,1\]\\. Expected target censoring proportion in
  the simulated data.

- event_survival:

  Function. Event-time survival function \\S\_{T^\*\mid\mathbf{X}}(t)\\.
  The first argument must be `t` and further arguments are considered
  covariates. Must be vectorized in `t`.

- covariate_density:

  Optional function. Joint density of the covariate vector over
  `covariate_bounds`. If provided, it must integrate to one over the
  specified domain. If `NULL`, no covariates are integrated.

- covariate_bounds:

  Named list or `NULL`. Bounds for covariates used both to integrate
  `covariate_density`. Use numeric length-2 vectors for continuous
  bounds and lists of numeric values for discrete variables. Names must
  match the non-`t` arguments of `event_survival`. Must be `NULL` when
  `covariate_density` is `NULL`.

- cens_survival:

  Function. Censoring-time survival function
  \\S\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(t)\\ with first
  argument `t` and the second argument named as in `target`. Must be
  vectorized in `t`.

- cens_density:

  Function. Censoring-time density
  \\f\_{C\_{\text{rnd}}\mid\lambda\_{C\_{\text{rnd}}}}(t)\\ with first
  argument `t` and the second argument named as in `target`. Must be
  vectorized in `t` and be consistent with `cens_survival`.

- tau:

  Numeric scalar. Optional time point at which the expected censoring
  proportion is evaluated. If `NULL`, the expected censoring proportion
  is evaluated overall, otherwise it is evaluated at the specified time
  point.

- time_admin_cens:

  Numeric scalar. Administrative censoring time \\\tau\_{\text{adm}}\\.
  Use `Inf` for no administrative censoring.

- time_accrual:

  Numeric scalar. Accrual time span \\\tau\_{\text{acc}}\\. Use `0` for
  no accrual (default).

- target:

  Character scalar. Name of the censoring parameter inside
  `cens_survival` and `cens_density` (e.g., `"lambda_c"`).

- target_bounds:

  Numeric length-2 vector. Lower and upper bounds of the search interval
  for the parameter named in `target`.

- t_min:

  Numeric scalar. Left boundary of the time support used for validation
  and numerical checks. Must be non-negative.

- t_max:

  Numeric scalar. Right boundary of the time support used for validation
  and numerical checks. Use `Inf` for an unbounded right tail.

- abs_tol:

  Numeric scalar. Absolute tolerance used for integration (on the
  probability scale). Must be positive and not too large (e.g., \< 0.1).

## Value

A list with components:

- `parameter`: the value of the censoring parameter.

- `cens_prop`: the achieved expected censoring proportion at the
  calculated parameter.

## Details

The parameter \\\lambda\_{C\_{\text{rnd}}}\\ is determined by solving
the root finding problem

\\P(\delta=0\mid\lambda\_{C\_{\text{rnd}}}) - p_C = 0\\

where \\P(\delta=0\mid\lambda\_{C\_{\text{rnd}}})\\ is the expected
censoring proportion over the covariate distributions and
\\p_C\in\[0,1\]\\ is the target censoring proportion. The expected
censoring proportion \\P(\delta=0\mid \lambda\_{C\_{\text{rnd}}})=
E_X\left\[P(\delta=0\mid
\mathbf{X},\lambda\_{C\_{\text{rnd}}})\right\]\\ is evaluated by
[`estimate_cens_prop()`](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/estimate_cens_prop.md)
and its mathematical background is detailed in the package description
([CensRCalc](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/CensRCalc.md)).
A bracketing root finder is used with the interval given by
`target_bounds`.

Convergence depends on the chosen bounds covering a sign change of the
objective function. The tolerance `abs_tol` is used for both integration
accuracy and the root-finding convergence criterion.

Progress reporting is supported via the progressr package. When a
progressr handler is enabled (for example,
[`progressr::with_progress()`](https://progressr.futureverse.org/reference/with_progress.html)),
the iterations of the bracketing solver emit progress updates.

## See also

[`estimate_cens_prop()`](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/estimate_cens_prop.md)
for the analytical evaluation of the expected censoring proportion, and
[CensRCalc](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/CensRCalc.md)
for mathematical details.

## Examples

``` r
# Target 30% censoring under exponential models without covariates
event_surv <- function(t) exp(-0.2 * t)

sol <- find_cens_param(
  cens_prop = 0.3,
  event_survival = event_surv
)

sol$parameter
#> [1] 0.08571428
```
