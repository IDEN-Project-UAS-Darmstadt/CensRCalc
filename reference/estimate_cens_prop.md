# Calculate expected censoring proportion

Evaluates the expected censoring proportion under independent censoring
by integrating survival and censoring functions over time and
covariates. Administrative censoring and accrual periods are supported.

## Usage

``` r
estimate_cens_prop(
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
  time_accrual = 0,
  target = "lambda_c",
  target_bounds = c(1e-05, 5),
  t_min = 0,
  t_max = Inf,
  abs_tol = 1e-06
)
```

## Arguments

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

  Numeric length-2 vector. Lower and upper bounds of the validation
  interval for the parameter named in `target`.

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

A callable function of two numeric scalars, one named by `target` and
the other `tau`, which returns the expected censoring proportion under
the specified model components. When `tau` is specified, the expected
censoring proportion is evaluated at the time point `tau` instead of
overall time.

The returned function has class `"cens_prop_fun"` and supports S3
methods such as
[`print()`](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/cens_prop_methods.md)
and
[`plot()`](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/cens_prop_methods.md)
for summarizing and visualizing the censoring model.

## Details

This function returns the expected censoring proportion over the
covariates
\\P(\delta=0\mid\lambda\_{C\_{\text{rnd}}})=E_X\left\[P(\delta=0\mid
\mathbf{X},\lambda\_{C\_{\text{rnd}}})\right\]\\ obtained by numerical
integration of survival and censoring functions over time and finally
integrating over the covariate distributions. The mathematical
background and the full expression for
\\P(\delta=0\mid\lambda\_{C\_{\text{rnd}}})\\ are provided in the
package description, see
[CensRCalc](https://iden-project-uas-darmstadt.github.io/CensRCalc/reference/CensRCalc.md).

Functions supplied for time must be vectorized in their first argument
(time). The survival functions and densities of the event time and
censoring time variables, respectively, must be consistent, i.e., they
correspond to the same distribution.

## Examples

``` r
# Exponential event time without covariates
event_surv <- function(t) exp(-0.2 * t)

cens_surv <- function(t, lambda_c) exp(-lambda_c * t)
cens_dens <- function(t, lambda_c) lambda_c * exp(-lambda_c * t)

f_obj <- estimate_cens_prop(
  event_survival = event_surv,
  covariate_density = NULL,
  covariate_bounds = NULL,
  cens_survival = cens_surv,
  cens_density = cens_dens,
  time_admin_cens = 10,
  time_accrual = 0,
)

print(f_obj)
#> <Expected censoring proportion function>
#> Target parameter : lambda_c 
#> Target range     : 1e-05 - 5 
#> Administrative censoring: 10 
#> Accrual period   : 0 
f_obj(0.3)
#> [1] 0.5959572
f_obj(0.3, tau = 5)
#> [1] 0.550749
```
