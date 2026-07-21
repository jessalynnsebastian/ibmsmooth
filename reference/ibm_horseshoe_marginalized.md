# Fit a marginalized Gaussian horseshoe IBM

Fits the same Gaussian-data model as the current noncentered horseshoe
Stan model, but integrates the two-dimensional IBM state out with a
Kalman filter during HMC. Exact conditional state trajectories are
reconstructed afterward with a simulation smoother.

## Usage

``` r
ibm_horseshoe_marginalized(
  t,
  y,
  infer_at = NULL,
  log_sigma = list(mu = -1, sd = 1),
  zeta = 0.1,
  initial_sd = 5,
  n_state_draws = 1000L,
  seed = NULL,
  iter = 2000,
  chains = 4,
  cores = getOption("mc.cores", chains),
  max_treedepth = 12,
  adapt_delta = 0.9,
  ...
)
```

## Arguments

- t:

  Numeric vector of observation times.

- y:

  Numeric vector of observations.

- infer_at:

  Optional additional times at which state draws are required.

- log_sigma:

  Prior specification used by the existing Stan backend.

- zeta:

  Global horseshoe scale.

- initial_sd:

  Scale of the correlated unit-time IBM initial-state prior.

- n_state_draws:

  Number of conditional state trajectories to reconstruct. Set to zero
  to retain only hyperparameter draws.

- seed:

  Optional random seed used by Stan and the simulation smoother.

- iter, chains, cores:

  Stan sampling controls.

- max_treedepth, adapt_delta:

  Stan HMC controls.

- ...:

  Additional arguments passed to
  [`rstan::stan()`](https://mc-stan.org/rstan/reference/stan.html).

## Value

An object of class `ibmfit` containing the marginalized Stan fit and
cached posterior state draws.

## Details

The model is also available through
[`ibm()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm.md)
by selecting the marginalized horseshoe engine.
