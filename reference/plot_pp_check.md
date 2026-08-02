# Visual posterior predictive check

Visual posterior predictive check

## Usage

``` r
plot_pp_check(
  ibmfit,
  type = c("interval", "density"),
  n_samples = 500,
  probs = c(0.05, 0.5, 0.95),
  n_density = 30,
  seed = NULL
)
```

## Arguments

- ibmfit:

  A posterior `ibmfit`.

- type:

  Either `"interval"` for observation-level predictive intervals or
  `"density"` for observed and replicated-data densities.

- n_samples:

  Number of replicated data sets.

- probs:

  Predictive interval probabilities for `type = "interval"`.

- n_density:

  Number of replicated densities to draw.

- seed:

  Optional random seed.

## Value

A ggplot object.
