# Fit an integrated Brownian motion smoother

Fits either ordinary integrated Brownian motion (IBM) with one diffusion
scale or the locally adaptive IBM described by
\\df'(t)=\gamma\xi(t)dW(t)\\, \\df(t)=f'(t)dt\\. In the adaptive model,
each interval has a half-Cauchy local scale \\\xi_j\\.

## Usage

``` r
ibm(
  t = NULL,
  y = NULL,
  infer_at = NULL,
  adaptive = FALSE,
  smoothing_prior = "tau ~ lognormal(-2, 0.5);",
  global_scale = 0.1,
  log_sigma = list(mu = -1, sd = 1),
  initial_sd = 5,
  iter = 2000,
  chains = 4,
  cores = getOption("mc.cores", chains),
  max_treedepth = 12,
  adapt_delta = 0.9,
  get_code = FALSE,
  ...
)
```

## Arguments

- t:

  Numeric observation locations.

- y:

  Numeric observations. Replicates at the same location are allowed.

- infer_at:

  Optional additional locations at which to infer the state.

- adaptive:

  Logical; use locally adaptive IBM when `TRUE`.

- smoothing_prior:

  A character string containing valid Stan code for the prior on the
  positive ordinary-IBM smoothing parameter `tau`. For example,
  `"tau ~ normal(0, 0.5);"` or `"tau ~ student_t(3, 0, 1);"`. This is
  used only when `adaptive = FALSE`.

- global_scale:

  Positive half-Cauchy scale for `gamma` in the adaptive global–local
  prior.

- log_sigma:

  List with `mu` and positive `sd` for the log observation standard
  deviation on the internally standardized response scale.

- initial_sd:

  Positive prior standard deviation for the initial function value and
  derivative on the internally standardized scale.

- iter, chains, cores:

  Stan sampling controls.

- max_treedepth, adapt_delta:

  Stan HMC controls.

- get_code:

  If `TRUE`, return the selected Stan program rather than fit.

- ...:

  Additional arguments passed to
  [`rstan::stan()`](https://rdrr.io/pkg/rstan/man/stan.html).

## Value

An object of class `ibmfit`, or Stan code when `get_code = TRUE`.

## Details

Time is shifted and divided by its mean grid spacing and the response is
standardized before fitting. Consequently, `smoothing_prior` is
expressed on that stable internal scale. Natural-scale draws are
returned by the extraction helpers.

The adaptive transition over an interval of length \\\Delta_j\\ is
exactly bivariate normal with mean \\(f'\_j, f_j+\Delta_j f'\_j)\\ and
covariance \\\gamma^2\xi_j^2 Q(\Delta_j)\\, where
\\Q(\Delta)=((\Delta,\Delta^2/2),(\Delta^2/2,\Delta^3/3))\\.
