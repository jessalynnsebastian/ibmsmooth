# ibmsmooth

`ibmsmooth` fits integrated Brownian motion smoothers for conditionally
independent observations using Stan. The package contains two latent models:

- ordinary IBM with one positive smoothing parameter and a user-supplied Stan
  prior for that parameter; and
- locally adaptive IBM with horseshoe global--local shrinkage on the derivative
  diffusion rate.

For locations \(t_j<t_{j+1}\), \(\Delta_j=t_{j+1}-t_j\), and state
\(x_j=(f'(t_j),f(t_j))^\mathsf T\), both models use the exact transition

\[
x_{j+1}\mid x_j,\lambda_j \sim
N_2\left[
\begin{pmatrix}f'(t_j)\\f(t_j)+\Delta_jf'(t_j)\end{pmatrix},
\lambda_j
\begin{pmatrix}
\Delta_j & \Delta_j^2/2\\
\Delta_j^2/2 & \Delta_j^3/3
\end{pmatrix}
\right].
\]

Ordinary IBM sets \(\lambda_j=\tau^2\). The adaptive model sets
\(\lambda_j=\gamma^2\xi_j^2\), with
\(\xi_j\sim\mathrm{C}^+(0,1)\) and
reference-process-calibrated half-Cauchy or half-normal prior on \(\gamma\).

After mapping time to \([0,1]\), `global_upper = U` and
`global_alpha = alpha` set the prior scale through
\(\Pr\{\gamma/\sqrt{3}>U\}=\alpha\). This uses the maximum reference-process
standard deviation and is independent of the observation grid.

## Installation

```r
remotes::install_github("jessalynnsebastian/ibmsmooth")
```

## Examples

```r
library(ibmsmooth)

fit <- ibm(t, y)

# Any valid Stan prior statement involving the positive parameter `tau`.
fit_half_t <- ibm(
  t, y,
  smoothing_prior = "tau ~ student_t(3, 0, 0.5);"
)

fit_adaptive <- ibm(
  t, y,
  adaptive = TRUE,
  global_prior = "half_cauchy",
  global_upper = 1,
  global_alpha = 0.05
)

# Centered implementation of the same model
fit_adaptive_centered <- ibm(
  t, y,
  adaptive = TRUE,
  parameterization = "centered",
  global_prior = "half_normal",
  global_upper = 1,
  global_alpha = 0.05
)

# Optional finite-slab regularized horseshoe
fit_adaptive_regularized <- ibm(
  t, y,
  adaptive = TRUE,
  adaptive_prior = "regularized_horseshoe",
  slab_scale = 1,
  slab_df = 4
)

# The default regularized_retry = "ask" also detects the specific combination
# of poor HMC diagnostics and local scales pinned near the half-Cauchy boundary.
# In an interactive session it offers to run the regularized fit above.

plot(fit_adaptive)
plot_diffusion(fit_adaptive)

# Visual model checks
plot_pp_check(fit_adaptive)
plot_pp_check(fit_adaptive, type = "density")
plot_prior_posterior(fit_adaptive, seed = 123)

# Timing, ESS, R-hat, divergences, treedepth, acceptance, and E-BFMI
diagnostics <- get_diagnostics(fit_adaptive)
diagnostics$overview
```

## Observation families

Use `family` to select a likelihood. Available families are `"gaussian"`,
`"student_t"`, `"bernoulli"`, `"binomial"`, `"poisson"`,
`"negative_binomial"`, `"lognormal"`, `"gamma"`, `"beta"`, and
`"exponential"`.

```r
fit_counts <- ibm(t, counts, adaptive = TRUE, family = "poisson",
                  exposure = person_time)
fit_binary <- ibm(t, events, family = "binomial", trials = sample_sizes)
```

Gaussian and Student-t responses are standardized internally and `f` is
returned on the original response scale. For other families, `f` is the latent
predictor: logit scale for Bernoulli, binomial, and beta means, and log-mean
scale for Poisson, negative-binomial, gamma, and exponential observations.
For lognormal observations, `f` is the mean on the log scale. Consequently,
`global_upper` is interpreted on the latent-predictor scale for these families.

The generated Stan program can be used without supplying data:

```r
code <- ibm(family = "poisson", adaptive = TRUE, get_code = TRUE)
ibm(family = "poisson", adaptive = TRUE, print_code = TRUE)
ibm(family = "poisson", adaptive = TRUE, stan_file = "ibm_poisson.stan")
```

When the truth is known, such as in a simulation study, function and
derivative performance use the same interface:

```r
summarize_performance(
  fit_adaptive,
  truth = list(f = truth_f, fprime = truth_fprime),
  new_t = evaluation_grid,
  level = 0.95
)
```

This reports mean absolute deviation of the posterior median (MAD), mean
credible interval width (MCIW), coverage percentage, and empirical CRPS.
Truth entries may instead be functions of `t`.

Time is standardized internally. Gaussian and Student-t responses are also
standardized, while other families retain their latent link scale. Custom
smoothing priors are specified on that internal model scale. Extraction
functions transform applicable posterior draws back to the original response
and time scales. Locations supplied through
`infer_at` are not added to the Stan state vector; `predict_curve()` samples
them afterward using exact conditional IBM bridges.

The default `parameterization = "noncentered"` is generally preferable when
most interval innovations are weakly informed or strongly shrunk.
`parameterization = "centered"` can be preferable when the states and one or
more large changes are strongly identified by the likelihood. Both use the
same exact IBM transitions and priors; they differ only in their HMC
coordinates.

The adaptive default remains `adaptive_prior = "horseshoe"`. The optional
regularized horseshoe retains shrinkage near zero but caps the effective local
diffusion scale with a Student-t slab. `slab_scale` is specified on the internal
latent derivative-diffusion scale.

The post-fit recommendation can also be inspected directly:

```r
check_regularized_horseshoe(fit_adaptive)
```

Interactive retries return the regularized fit and store the triggering
ordinary-horseshoe diagnostic summary in `fit$retry_info`. In scripts and other
non-interactive sessions, the package prints the recommendation but never
starts another fit automatically. Use `regularized_retry = "never"` to disable
the check and prompt entirely.

## License

This project is released under the MIT license. See [LICENSE](LICENSE).
