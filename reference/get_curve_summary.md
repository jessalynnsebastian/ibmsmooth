# Summarize posterior curves

Compute pointwise posterior summaries for the latent function and
derivative.

## Usage

``` r
get_curve_summary(
  ibmfit,
  n_samples = 1000,
  probs = c(0.025, 0.1, 0.5, 0.9, 0.975)
)
```

## Arguments

- ibmfit:

  An object of class `ibmfit`.

- n_samples:

  Number of samples to use for INLA fits, or maximum number of MCMC
  draws to keep for Stan fits.

- probs:

  Numeric vector of posterior probabilities.

## Value

A data frame with one row per parameter and time point.
