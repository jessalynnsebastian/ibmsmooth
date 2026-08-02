# Fit an IBM smoother using Stan

`ibm_smooth()` is an alias for
[`ibm()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm.md)
retained as the package's lower-level fitting name.

## Usage

``` r
ibm_smooth(
  t = NULL,
  y = NULL,
  infer_at = NULL,
  adaptive = FALSE,
  family = c("gaussian", "student_t", "bernoulli", "binomial", "poisson",
    "negative_binomial", "lognormal", "gamma", "beta", "exponential"),
  trials = NULL,
  exposure = 1,
  student_df = 4,
  parameterization = c("noncentered", "centered"),
  smoothing_prior = "tau ~ lognormal(-2, 0.5);",
  global_prior = c("half_cauchy", "half_normal"),
  global_upper = 1,
  global_alpha = 0.05,
  adaptive_prior = c("horseshoe", "regularized_horseshoe"),
  slab_scale = 1,
  slab_df = 4,
  regularized_retry = c("ask", "never"),
  log_sigma = list(mu = -1, sd = 1),
  log_phi = list(mu = 0, sd = 1),
  initial_sd = 5,
  iter = 2000,
  chains = 4,
  cores = getOption("mc.cores", chains),
  init = "random",
  max_treedepth = 12,
  adapt_delta = 0.9,
  get_code = FALSE,
  print_code = FALSE,
  stan_file = NULL,
  ...
)
```

## Arguments

- t:

  Numeric observation locations.

- y:

  Numeric observations. Replicates at the same location are allowed.

- infer_at:

  Optional additional locations at which to predict the state after
  sampling. These locations do not enlarge the Stan latent state.

- adaptive:

  Logical; use locally adaptive IBM when `TRUE`.

- family:

  Conditionally independent observation family. Available families are
  `"gaussian"`, `"student_t"`, `"bernoulli"`, `"binomial"`, `"poisson"`,
  `"negative_binomial"`, `"lognormal"`, `"gamma"`, `"beta"`, and
  `"exponential"`.

- trials:

  Binomial trial counts. Required for `family = "binomial"`.

- exposure:

  Positive observation exposures for Poisson and negative binomial
  models. A scalar is recycled. The default is one.

- student_df:

  Positive degrees of freedom for the Student-t likelihood.

- parameterization:

  Either `"noncentered"` (the default) or `"centered"`. Both define the
  same IBM model with different HMC coordinates.

- smoothing_prior:

  A character string containing valid Stan code for the prior on the
  positive ordinary-IBM smoothing parameter `tau`. For example,
  `"tau ~ normal(0, 0.5);"` or `"tau ~ student_t(3, 0, 1);"`. This is
  used only when `adaptive = FALSE`.

- global_prior:

  Prior family for the adaptive global scale `gamma`: either
  `"half_cauchy"` or `"half_normal"`.

- global_upper:

  Positive upper bound `U` used to calibrate the global prior through
  `Pr(gamma / sqrt(3) > U) = global_alpha`. It is in
  standardized-response units for Gaussian and Student-t models and
  latent predictor units otherwise.

- global_alpha:

  Tail probability used in the global-prior calibration.

- adaptive_prior:

  Adaptive local-scale prior. `"horseshoe"` is the default;
  `"regularized_horseshoe"` adds a finite Student-t slab.

- slab_scale:

  Positive slab scale for the regularized horseshoe, on the internal
  latent derivative-diffusion scale.

- slab_df:

  Positive slab degrees of freedom. The default is 4.

- regularized_retry:

  Either `"ask"` (the default) or `"never"`. After a problematic
  interactive ordinary-horseshoe fit, `"ask"` offers to rerun the same
  model and HMC settings with a regularized horseshoe. Non-interactive
  sessions never start a refit automatically.

- log_sigma:

  List with `mu` and positive `sd` for the log observation standard
  deviation for Gaussian, Student-t, and lognormal likelihoods.

- log_phi:

  List with `mu` and positive `sd` for log `phi`: the negative binomial
  shape, gamma shape, or beta precision parameter.

- initial_sd:

  Positive prior standard deviation for the initial function value and
  derivative on the internal latent scale.

- iter, chains, cores:

  Stan sampling controls.

- init:

  Initial values passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

- max_treedepth, adapt_delta:

  Stan HMC controls.

- get_code:

  If `TRUE`, return the selected Stan program rather than fit.

- print_code:

  If `TRUE`, print the selected Stan program. This works without
  supplying `t` or `y`.

- stan_file:

  Optional path at which to save the selected Stan program. This also
  works without supplying data.

- ...:

  Additional arguments passed to
  [`rstan::stan()`](https://mc-stan.org/rstan/reference/stan.html).

## Value

An object of class `ibmfit`, or Stan code when `get_code = TRUE`.
