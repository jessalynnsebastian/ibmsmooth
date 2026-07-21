# Summarize an ibmfit object

Compute pointwise posterior summaries for the latent function, its
derivative, and the available model hyperparameters. Unlike a scalar
parameter summary, the curve summaries are returned at every point of
the latent time grid.

## Usage

``` r
# S3 method for class 'ibmfit'
summary(object, ..., n_samples = 1000, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
```

## Arguments

- object:

  An object of class `"ibmfit"`.

- ...:

  Additional arguments. Currently passed to curve-summary methods.

- n_samples:

  Number of posterior samples to draw for INLA summaries. For Stan fits,
  existing MCMC draws are used; if `n_samples` is smaller than the
  number of available draws, the first `n_samples` draws are used.

- probs:

  Numeric vector of posterior probabilities to summarize.

## Value

A list with components `curve` and `hyperparameters`.
