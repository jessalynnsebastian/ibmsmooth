# Summarize function and derivative estimation performance

Computes mean absolute deviation (MAD) of the posterior median, mean
credible interval width (MCIW), empirical coverage percentage, and the
empirical continuous ranked probability score (CRPS), pointwise and
averaged across evaluation locations.

## Usage

``` r
summarize_performance(
  ibmfit,
  truth,
  new_t = NULL,
  level = 0.95,
  n_samples = 1000,
  seed = NULL,
  pointwise = FALSE
)
```

## Arguments

- ibmfit:

  An `ibmfit` object.

- truth:

  Named list containing `f` and/or `fprime`. Each element can be a
  numeric vector evaluated at `new_t` or a function that accepts
  `new_t`. For non-Gaussian families these are truths on the latent
  predictor scale.

- new_t:

  Optional finite evaluation locations. By default, uses the fitted
  prediction grid.

- level:

  Credible interval confidence level strictly between zero and one.

- n_samples:

  Maximum number of posterior draws.

- seed:

  Optional seed for conditional bridge draws at new locations.

- pointwise:

  If `TRUE`, return both aggregate and location-specific results;
  otherwise return only the aggregate data frame.

## Value

A data frame with one row per supplied truth, or a list containing
`aggregate` and `pointwise` data frames when `pointwise = TRUE`.
