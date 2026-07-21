# Fit an integrated Brownian motion smoother using INLA

Fit a global or locally adaptive integrated Brownian motion smoother
using INLA. The latent field is the state vector `(f, fprime)` on the
latent time grid. The observation model is Gaussian.

## Usage

``` r
ibm_inla_fit(
  t,
  y,
  infer_at = NULL,
  adaptive = FALSE,
  n_knots = 5L,
  bs_degree = 3L,
  adaptive_deviation_sd = 0.5,
  log_process_sd = list(mu = -2, sd = 0.5),
  log_sigma = list(mu = -1, sd = 1),
  theta_mu = NULL,
  theta_sd = NULL,
  sigma_mu = NULL,
  sigma_sd = NULL,
  initial_sd = 5,
  initial_slope_sd = NULL,
  control.compute = list(config = TRUE),
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

- adaptive:

  Logical. If `FALSE`, use one global diffusion-precision parameter. If
  `TRUE`, model the log diffusion precision as a global log precision
  plus centered B-spline deviations shrunk toward zero.

- n_knots:

  Number of internal knots for the B-spline basis when
  `adaptive = TRUE`.

- bs_degree:

  Degree of the B-spline basis.

- adaptive_deviation_sd:

  Prior standard deviation for the centered B-spline deviation
  coefficients in adaptive INLA fits. The adaptive model uses
  `log(lambda_i) = alpha + B_i beta`, where `alpha` is the global log
  diffusion precision and the components of `beta` have independent
  `N(0, adaptive_deviation_sd^2)` priors. Smaller values shrink the
  model more strongly toward the global-precision IBM.

- log_process_sd:

  List with elements `mu` and `sd` giving the prior mean and standard
  deviation for the log IBM process standard deviation on the scaled
  model scale. The default matches the Stan backend:
  `list(mu = -2, sd = 0.5)`.

- log_sigma:

  List with elements `mu` and `sd` giving the prior mean and standard
  deviation for the log Gaussian observation standard deviation on the
  scaled model scale. The default matches the Stan backend:
  `list(mu = -1, sd = 1)`.

- theta_mu:

  Optional advanced override for the mean of the internal INLA global
  log diffusion-precision parameter. If `NULL`, this is set from
  `log_process_sd`.

- theta_sd:

  Optional advanced override for the standard deviation of the internal
  INLA global log diffusion-precision parameter. If `NULL`, this is set
  from `log_process_sd`.

- sigma_mu:

  Optional advanced override for the mean of the Gaussian likelihood log
  precision. If `NULL`, this is set from `log_sigma`.

- sigma_sd:

  Optional advanced override for the standard deviation of the Gaussian
  likelihood log precision. If `NULL`, this is set from `log_sigma`.

- initial_sd:

  Prior scale for the joint initial state `c(f_1, fprime_1)`. The
  initial covariance is
  `initial_sd^2 * matrix(c(1/3, 1/2, 1/2, 1), 2, 2)`, the unit-time
  integrated Brownian motion covariance for the level and derivative.

- initial_slope_sd:

  Deprecated. Retained only for backward compatibility; ignored by the
  INLA backend.

- control.compute:

  List passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) as
  `control.compute`. Defaults to `list(config = TRUE)` so that posterior
  sampling works for plotting, summaries, and
  [`get_samples_inla()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/get_samples_inla.md).

- ...:

  Additional arguments passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

## Value

An object of class `ibmfit`.
