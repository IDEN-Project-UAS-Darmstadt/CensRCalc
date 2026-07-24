# CensRCalc (development version)

Release candidate of 0.1.1.900

## Major changes

None

## Minor improvements and bug fixes

* Fixed a bug in the `integral_with_discr()` function, where the error was
  relative instead of absolute.
* Fixed a bug in the `progressr_uniroot()` function, where the last reported
  progress was updated incorrectly.

# CensRCalc 0.1.0

## Major changes

* Added parameter `tau`, which allows evaluating the censoring proportion at a
  specific time point $\tau$.
* Added `plot()` and `print()` methods for the expected censoring proportion
  function.
* Added tests.
* Better handling and guarantee for `abs_tol` for numerical integration.

## Minor improvements and bug fixes

* Added evaluation of the passed event survival function, using the marginal.
* Removed renv usage.
* Development workflow improvements and fixes.
* Documentation improvements and fixes.

# CensRCalc 0.0.1

## Major changes

* Initial package functionality

## Minor improvements and bug fixes

* None
