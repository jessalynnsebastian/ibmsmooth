# Summarize posterior curves

Summarize posterior curves

## Usage

``` r
get_curve_summary(ibmfit, n_samples = 1000, probs = c(0.025, 0.5, 0.975))
```

## Arguments

- ibmfit:

  An `ibmfit` object.

- n_samples:

  Maximum number of draws.

- probs:

  Posterior probabilities.

## Value

A data frame with pointwise summaries.
