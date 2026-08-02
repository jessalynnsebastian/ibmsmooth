# Summarize posterior curves

Summarize posterior curves

## Usage

``` r
get_curve_summary(
  ibmfit,
  n_samples = 1000,
  probs = c(0.025, 0.5, 0.975),
  param = c("both", "f", "fprime"),
  new_t = NULL,
  seed = NULL
)
```

## Arguments

- ibmfit:

  An `ibmfit` object.

- n_samples:

  Maximum number of draws.

- probs:

  Posterior probabilities.

- param:

  One of `"f"`, `"fprime"`, or `"both"`.

- new_t:

  Optional finite locations at which to summarize either curve.

- seed:

  Optional seed for conditional bridge draws at new locations.

## Value

A data frame with pointwise summaries.
