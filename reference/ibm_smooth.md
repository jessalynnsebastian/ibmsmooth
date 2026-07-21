# Fit an integrated Brownian motion (IBM) smoother using Stan

Fit a Bayesian integrated Brownian motion smoothing model to
observations (t, y) using Stan. The function rescales time and response
for more efficient sampling, constructs the data list for Stan, and runs
rstan::stan() on a centered or non-centered Stan implementation.

## Usage

``` r
ibm_smooth(
  t = NULL,
  y = NULL,
  infer_at = NULL,
  get_code = FALSE,
  write_to_file = NULL,
  adaptive = "nonadaptive",
  log_sigma = list(mu = -1, sd = 1),
  log_tau = list(mu = -2, sd = 0.5),
  zeta = 0.1,
  slab_scale = 2,
  alpha = 0.5,
  initial_sd = 5,
  centered = NULL,
  noncentered = NULL,
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

  Numeric vector of time points.

- y:

  Numeric vector of observations at times `t`. Must be same length as
  `t`.

- infer_at:

  Optional numeric vector of time points to do function inference at (in
  addition to just the observed time points). If not `NULL` these will
  be included in the Stan data and the fitted function will be returned
  at these points as well.

- get_code:

  Logical; if `TRUE` return the Stan model code as a single string
  (after potentially writing to `write_to_file`). If `get_code = TRUE`
  and both `t` and `y` are `NULL` the function will return the Stan code
  without attempting to fit a model.

- write_to_file:

  Optional character path. If provided and `get_code = TRUE` the Stan
  code is written to this file.

- adaptive:

  Character; specifies if and how to use adaptive smoothing. Options are
  `"nonadaptive"` (default), `"rw"`, `"horseshoe"`,
  `"baseline_horseshoe"`, `"rhs"`, or `"bridge"`.

- log_sigma:

  List with elements `mu` and `sd` giving the prior mean and standard
  deviation for `log(sigma)`. Defaults to `list(mu = -1, sd = 1)`.

- log_tau:

  List with elements `mu` and `sd` giving the prior mean and standard
  deviation for `log(tau)`. Used for `adaptive = "nonadaptive"`, `"rw"`,
  and `"baseline_horseshoe"`. For the latter this is the prior on the
  positive baseline process scale `tau0`. Defaults to
  `list(mu = -2, sd = 0.5)`.

- zeta:

  Numeric; global shrinkage scale used by horseshoe, baseline_horseshoe,
  rhs and bridge priors. In the baseline-horseshoe model it controls
  only excess roughness. Default `0.1`.

- slab_scale:

  Numeric; slab scale for the regularized horseshoe (rhs). Default `2`.

- alpha:

  Numeric in (0,2\]; exponent for bridge prior (only used when
  `adaptive = "bridge"`). Default `0.5`.

- initial_sd:

  Prior scale for the correlated initial state. The Stan backend now
  uses the same unit-time IBM initial covariance as the INLA backend:
  for `c(f_1, fprime_1)`, the covariance is
  `initial_sd^2 * matrix(c(1/3, 1/2, 1/2, 1), 2, 2)` on the scaled
  response/time scale.

- centered:

  Logical; if `TRUE` use the centered Stan parameterization. If `NULL` a
  sensible default is chosen: centered is `TRUE` when
  `adaptive == "nonadaptive"` (except that `"bridge"` forces centered =
  TRUE). You can alternatively pass `noncentered` (logical) to
  explicitly request the non-centered form.

- noncentered:

  Deprecated alias for specifying non-centered parameterization.

- iter:

  Integer; total number of iterations for each chain (including warmup).
  Default `2000`.

- chains:

  Integer; number of MCMC chains. Default `4`.

- cores:

  Integer; number of parallel processes to use for Stan chains. Defaults
  to `getOption("mc.cores", chains)`.

- max_treedepth:

  Integer; maximum tree depth for Stan's NUTS sampler. Default `12`.

- adapt_delta:

  Numeric in (0,1); target average acceptance probability for Stan's
  NUTS. Default `0.85`.

- ...:

  Additional arguments passed on to
  [`rstan::stan()`](https://mc-stan.org/rstan/reference/stan.html).

## Value

An object of class `ibmfit` (a list) with components:

- stanfit:

  The `stanfit` object returned by
  [`rstan::stan()`](https://mc-stan.org/rstan/reference/stan.html).

- data:

  A list containing the original data and scaling constants used:
  `y_raw`, `y_mean`, `y_sd`, `t_raw`, `t_min`, and `dt_mean`.

## Details

This function supports several adaptive smoothing options:

- "nonadaptive" (default): single global smoothing parameter (tau).

- "rw": random-walk prior on log(tau) across segments.

- "horseshoe": horseshoe prior on local IBM process scales.

- "baseline_horseshoe": horseshoe excess roughness above a positive IBM
  baseline.

- "rhs": regularized horseshoe prior on local IBM process scales.

- "bridge": bridge prior on whitened IBM transition innovations
  (requires centered parameterization).

The function rescales `y` to zero mean and unit variance and rescales
time so that the mean time step is near 1. A flag `regular` is computed
internally to indicate whether time points are equally spaced. The
transition covariance for interval `i` is `tau_i^2 S_i`, where `S_i` is
the IBM covariance determined by `deltat[i]`. INLA uses the equivalent
precision parameterization `lambda_i = 1 / tau_i^2`.

If `get_code = TRUE` the Stan code used for the chosen combination of
`adaptive` and `centered` is returned as a single string; if
`write_to_file` is provided the code will also be written to that path.

## See also

[`stan`](https://mc-stan.org/rstan/reference/stan.html)
