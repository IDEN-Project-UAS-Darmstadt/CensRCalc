# Print and plot methods for censoring proportion functions

Print and plot methods for censoring proportion functions

## Usage

``` r
# S3 method for class 'cens_prop_fun'
print(x, ...)

# S3 method for class 'cens_prop_fun'
plot(x, n = 100, tau = NULL, ...)
```

## Arguments

- x:

  A `"cens_prop_fun"` object.

- ...:

  Additional arguments.

- n:

  Number of evaluation points.

- tau:

  Optional time point for evaluating the expected censoring proportion.
  If `NULL`, the censoring proportion is evaluated overall, otherwise it
  is evaluated at the specified time point.

## Value

The input object, invisibly.
