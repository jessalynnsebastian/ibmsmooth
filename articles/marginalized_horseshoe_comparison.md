# Prototype: marginalized Gaussian horseshoe IBM

This vignette benchmarks a marginalized Stan model that analytically
integrates out the Gaussian IBM state. It is intentionally separate from
the public
[`ibm()`](https://jessalynnsebastian.github.io/ibmsmooth/reference/ibm.md)
interface. The existing model samples the function, derivative, and
horseshoe scales jointly; the prototype uses an exact two-dimensional
Kalman filter and asks HMC to sample only the observation noise and
horseshoe scales.

The two models therefore target the same hyperparameter posterior. The
marginalized fit can now recover exact conditional draws of the function
and derivative afterward with a simulation smoother.

## Check the marginal likelihood

For a small problem, the exported validation helper independently
constructs the dense observation covariance and compares its
multivariate Gaussian log likelihood with the Kalman-filter result. The
difference should be close to machine precision.

``` r

set.seed(4)
t_check <- sort(runif(8))
y_check <- rnorm(8)
tau_check <- exp(rnorm(7, -2, 0.3))

validate_marginalized_horseshoe(
  t_check, y_check, sigma = 0.25, tau = tau_check
)
```

    ##        kalman         dense    difference 
    ## -5.660236e+01 -5.660236e+01  4.646949e-12

## Common initial-state convention

The centered, noncentered, and marginalized horseshoe models now use the
same correlated unit-time IBM prior:

``` math
\begin{pmatrix}f_1\\f'_1\end{pmatrix}
\sim N\left\{0,\ \texttt{initial\_sd}^2
\begin{pmatrix}1/3&1/2\\1/2&1\end{pmatrix}\right\}.
```

This is also the convention used by the INLA implementation.
Standardizing the prior removes an earlier difference between the
centered and noncentered Stan files and makes their posterior comparison
meaningful.

## Simulate data

``` r

set.seed(11)
n <- 80
t <- sort(runif(n, 0, 10))
truth <- sin(t) + 0.7 * plogis(12 * (t - 5))
y <- truth + rnorm(n, 0, 0.2)
```

## Prepare the same scaled data used by `ibm_smooth()`

The following helper is local to this vignette because it prepares data
for the low-level Stan comparison rather than the ordinary package
interface.

``` r

make_horseshoe_stan_data <- function(t, y, infer_at = NULL,
                                     log_sigma = list(mu = -1, sd = 1),
                                     zeta = 0.1, initial_sd = 5) {
  stopifnot(length(t) == length(y), length(unique(t)) >= 2)

  grid_raw <- sort(unique(c(t, infer_at)))
  dt_mean <- mean(diff(grid_raw))
  grid <- (grid_raw - min(grid_raw)) / dt_mean
  t_scaled <- (t - min(grid_raw)) / dt_mean
  ord <- order(match(t_scaled, grid))

  list(
    N_obs = length(y),
    T = length(grid),
    obs_time_idx = match(t_scaled, grid)[ord],
    deltat = diff(grid),
    y_obs = as.numeric(scale(y))[ord],
    regular = as.integer(all(abs(diff(grid) - diff(grid)[1]) <
                               sqrt(.Machine$double.eps))),
    log_sigma_mu = log_sigma$mu,
    log_sigma_sd = log_sigma$sd,
    zeta = zeta,
    initial_sd = initial_sd
  )
}

stan_data <- make_horseshoe_stan_data(t, y)
```

## Fit both models

These chunks are not evaluated while building the package vignette
because compiling and sampling two Stan models is expensive. Run them
interactively from the package source directory.

``` r

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

existing_file <- system.file(
  "stan", "ibm_adaptive_horseshoe.stan", package = "ibmsmooth"
)
marginal_file <- system.file(
  "stan", "ibm_adaptive_horseshoe_gaussian_marginalized.stan",
  package = "ibmsmooth"
)

# When running directly from a source checkout rather than an installed package:
if (!nzchar(existing_file))
  existing_file <- "inst/stan/ibm_adaptive_horseshoe.stan"
if (!nzchar(marginal_file))
  marginal_file <- "inst/stan/ibm_adaptive_horseshoe_gaussian_marginalized.stan"

existing_model <- stan_model(existing_file)
marginal_model <- stan_model(marginal_file)

set.seed(2026)
existing_time <- system.time({
  existing_fit <- sampling(
    existing_model, data = stan_data,
    chains = 4, iter = 2000, seed = 2026,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
})

set.seed(2026)
marginal_time <- system.time({
  marginal_fit <- sampling(
    marginal_model,
    data = stan_data[names(stan_data) != "regular"],
    chains = 4, iter = 2000, seed = 2026,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
})
```

For an ordinary marginalized fit with reconstructed curve draws, use the
package wrapper instead of compiling the model manually:

``` r

fit_marginalized <- ibm_horseshoe_marginalized(
  t, y,
  n_state_draws = 1000,
  chains = 4,
  iter = 2000,
  seed = 2026,
  adapt_delta = 0.95
)

plots <- plot(fit_marginalized)
plots$function_plot
plots$derivative_plot

f_draws <- get_samples(fit_marginalized, "f")
fprime_draws <- get_samples(fit_marginalized, "fprime")
```

It is also available as an explicit opt-in engine through the unified
interface. The existing joint sampler remains the default.

``` r

fit_marginalized <- ibm(
  t, y,
  method = "stan",
  adaptive = TRUE,
  stan_adaptive_method = "horseshoe",
  stan_horseshoe_engine = "marginalized",
  n_state_draws = 1000,
  chains = 4,
  iter = 2000,
  seed = 2026,
  adapt_delta = 0.95
)
```

## Compare posterior agreement and sampling efficiency

The hyperparameter names are shared by the two models. Agreement in
`sigma`, `gamma`, and representative `tau` values checks the target
distribution; effective sample size per elapsed second measures whether
marginalization is actually useful.

``` r

pars <- c("sigma", "gamma", "tau[1]", paste0("tau[", stan_data$T - 1, "]"))

posterior_table <- function(fit, label) {
  s <- summary(fit, pars = pars)$summary
  data.frame(
    model = label,
    parameter = rownames(s),
    mean = s[, "mean"],
    sd = s[, "sd"],
    n_eff = s[, "n_eff"],
    Rhat = s[, "Rhat"],
    row.names = NULL
  )
}

comparison <- rbind(
  posterior_table(existing_fit, "existing noncentered"),
  posterior_table(marginal_fit, "marginalized")
)
comparison

elapsed <- c(
  existing = unname(existing_time["elapsed"]),
  marginalized = unname(marginal_time["elapsed"])
)
elapsed

aggregate(n_eff ~ model, comparison, sum)
```

Also inspect divergences rather than relying on runtime alone:

``` r

count_divergences <- function(fit) {
  sum(vapply(
    get_sampler_params(fit, inc_warmup = FALSE),
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))
}

c(
  existing = count_divergences(existing_fit),
  marginalized = count_divergences(marginal_fit)
)
```

The prototype is promising only if posterior summaries agree,
diagnostics do not deteriorate, and effective sample size per second
improves on datasets of the size the package is meant to handle.
Marginalization reduces the sampled dimension, but differentiating
through the Kalman recursion at every leapfrog step has its own cost, so
elapsed time by itself is not a sufficient comparison.

## Repeat the benchmark across scenarios

The first result can be unusually favorable or unfavorable. The
following scenarios exercise a globally smooth curve, a sharp change,
and locally varying roughness at several grid sizes. Keep the
data-generating seed distinct from the Stan seed so the experiment can
be repeated systematically.

``` r

simulate_scenario <- function(n, scenario, seed) {
  set.seed(seed)
  t <- sort(runif(n, 0, 10))
  truth <- switch(
    scenario,
    smooth = sin(t),
    jump = sin(t) + ifelse(t >= 5, 1, 0),
    varying = sin(t) + ifelse(t > 4 & t < 7, 0.3 * sin(12 * t), 0)
  )
  list(t = t, y = truth + rnorm(n, 0, 0.2), truth = truth)
}

benchmark_grid <- expand.grid(
  n = c(30, 100, 300),
  scenario = c("smooth", "jump", "varying"),
  data_seed = 1:3,
  stringsAsFactors = FALSE
)

# Use the fitting and diagnostic code above inside a loop over benchmark_grid.
# Save one row per fit containing elapsed time, divergences, maximum R-hat, and
# bulk/tail ESS per second. Curve RMSE and interval coverage can be computed
# from get_samples() for the marginalized wrapper and the existing ibmfit.
```

The simulation smoother itself is tested against the corresponding dense
conditional Gaussian posterior, including a case with replicated
observations and grid points introduced only through `infer_at`. This
separates correctness of state reconstruction from HMC performance.
