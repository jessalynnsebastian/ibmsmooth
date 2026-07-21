# Unified interface for integrated Brownian motion smoothing

Fit an integrated Brownian motion smoother using either Stan or INLA.

## Usage

``` r
ibm(
  t,
  y,
  infer_at = NULL,
  method = c("stan", "inla"),
  adaptive = FALSE,
  stan_adaptive_method = c("horseshoe", "baseline_horseshoe", "rw", "rhs", "bridge"),
  stan_horseshoe_engine = c("joint", "marginalized"),
  ...
)
```

## Arguments

- t:

  Numeric vector of time points.

- y:

  Numeric vector of observations at times `t`.

- infer_at:

  Optional numeric vector of additional time points at which inference
  should be returned.

- method:

  Computational backend, either `"stan"` or `"inla"`.

- adaptive:

  Logical. If `FALSE`, use a single smoothness parameter. If `TRUE`, fit
  a locally adaptive model.

- stan_adaptive_method:

  Stan-only adaptive prior. The default is `"horseshoe"`. Use
  `"baseline_horseshoe"` for independent horseshoe excess roughness
  above a positive global IBM baseline.

- stan_horseshoe_engine:

  Stan engine used for the adaptive horseshoe model. `"joint"` retains
  the existing sampler; `"marginalized"` uses the Gaussian-data Kalman
  marginalization and a simulation smoother.

- ...:

  Additional arguments passed to
  [`ibm_smooth()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm_smooth.md)
  or
  [`ibm_inla_fit()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm_inla_fit.md).

## Value

An object of class `ibmfit`.
