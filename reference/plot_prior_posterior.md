# Plot prior and posterior distributions

Plot prior and posterior distributions

## Usage

``` r
plot_prior_posterior(
  ibmfit,
  prior_fit = NULL,
  parameters = NULL,
  n_samples = 1000,
  ...
)
```

## Arguments

- ibmfit:

  A posterior `ibmfit`.

- prior_fit:

  Optional prior-only fit from
  [`sample_prior()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/sample_prior.md).
  When omitted, a prior fit is sampled automatically.

- parameters:

  Hyperparameters to display. The default shows an observation
  scale/shape parameter when present and the model's global smoothing
  scale (`tau` or `gamma`).

- n_samples:

  Maximum draws from each distribution.

- ...:

  Arguments passed to
  [`sample_prior()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/sample_prior.md)
  when `prior_fit` is omitted.

## Value

A faceted ggplot object.
