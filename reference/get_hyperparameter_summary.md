# Summarize posterior hyperparameters

Summarize posterior hyperparameters

## Usage

``` r
get_hyperparameter_summary(
  ibmfit,
  n_samples = 1000,
  probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
  natural = TRUE,
  include_internal = FALSE
)
```

## Arguments

- ibmfit:

  An object of class `ibmfit`.

- n_samples:

  Number of samples to draw for INLA fits, or maximum number of Stan
  draws to keep.

- probs:

  Numeric vector of posterior probabilities.

- natural:

  Logical. If `TRUE`, summarize hyperparameters on the original
  response/time scale when possible.

- include_internal:

  Logical. If `TRUE`, include internal backend hyperparameters in
  addition to natural-scale quantities.

## Value

A data frame with one row per hyperparameter component.
