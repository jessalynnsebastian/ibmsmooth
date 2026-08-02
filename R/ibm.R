#' Fit an integrated Brownian motion smoother
#'
#' Fits either ordinary integrated Brownian motion (IBM) with one diffusion
#' scale or the locally adaptive IBM described by
#' \eqn{df'(t)=\gamma\xi(t)dW(t)}, \eqn{df(t)=f'(t)dt}.  In the adaptive model,
#' each interval has a half-Cauchy local scale \eqn{\xi_j}.
#'
#' @param t Numeric observation locations.
#' @param y Numeric observations. Replicates at the same location are allowed.
#' @param family Conditionally independent observation family. Available
#'   families are `"gaussian"`, `"student_t"`, `"bernoulli"`, `"binomial"`,
#'   `"poisson"`, `"negative_binomial"`, `"lognormal"`, `"gamma"`, `"beta"`,
#'   and `"exponential"`.
#' @param trials Binomial trial counts. Required for `family = "binomial"`.
#' @param exposure Positive observation exposures for Poisson and negative
#'   binomial models. A scalar is recycled. The default is one.
#' @param student_df Positive degrees of freedom for the Student-t likelihood.
#' @param infer_at Optional additional locations at which to predict the state
#'   after sampling. These locations do not enlarge the Stan latent state.
#' @param adaptive Logical; use locally adaptive IBM when `TRUE`.
#' @param parameterization Either `"noncentered"` (the default) or
#'   `"centered"`. Both define the same IBM model with different HMC
#'   coordinates.
#' @param smoothing_prior A character string containing valid Stan code for the
#'   prior on the positive ordinary-IBM smoothing parameter `tau`. For example,
#'   `"tau ~ normal(0, 0.5);"` or `"tau ~ student_t(3, 0, 1);"`. This is used
#'   only when `adaptive = FALSE`.
#' @param global_prior Prior family for the adaptive global scale `gamma`:
#'   either `"half_cauchy"` or `"half_normal"`.
#' @param global_upper Positive upper bound `U` used to calibrate the global
#'   prior through `Pr(gamma / sqrt(3) > U) = global_alpha`. It is in
#'   standardized-response units for Gaussian and Student-t models and latent
#'   predictor units otherwise.
#' @param global_alpha Tail probability used in the global-prior calibration.
#' @param adaptive_prior Adaptive local-scale prior. `"horseshoe"` is the
#'   default; `"regularized_horseshoe"` adds a finite Student-t slab.
#' @param slab_scale Positive slab scale for the regularized horseshoe, on the
#'   internal latent derivative-diffusion scale.
#' @param slab_df Positive slab degrees of freedom. The default is 4.
#' @param regularized_retry Either `"ask"` (the default) or `"never"`.
#'   After a problematic interactive ordinary-horseshoe fit, `"ask"` offers
#'   to rerun the same model and HMC settings with a regularized horseshoe.
#'   Non-interactive sessions never start a refit automatically.
#' @param log_sigma List with `mu` and positive `sd` for the log observation
#'   standard deviation for Gaussian, Student-t, and lognormal likelihoods.
#' @param log_phi List with `mu` and positive `sd` for log `phi`: the negative
#'   binomial shape, gamma shape, or beta precision parameter.
#' @param initial_sd Positive prior standard deviation for the initial function
#'   value and derivative on the internal latent scale.
#' @param iter,chains,cores Stan sampling controls.
#' @param init Initial values passed to [rstan::sampling()].
#' @param max_treedepth,adapt_delta Stan HMC controls.
#' @param get_code If `TRUE`, return the selected Stan program rather than fit.
#' @param print_code If `TRUE`, print the selected Stan program. This works
#'   without supplying `t` or `y`.
#' @param stan_file Optional path at which to save the selected Stan program.
#'   This also works without supplying data.
#' @param ... Additional arguments passed to [rstan::stan()].
#'
#' @details
#' Time is mapped to the unit interval. Gaussian and Student-t responses are
#' centered and scaled to unit standard deviation; for other families, `f` is
#' represented on its usual link scale (logit for Bernoulli/binomial and beta
#' means; log mean for count, gamma, and exponential observations; and mean-log
#' for lognormal observations).
#' Locations in `infer_at` are evaluated after fitting
#' using exact conditional IBM bridges, avoiding additional Stan parameters.
#'
#' For an adaptive fit, the unit-diffusion IBM departure from its initial
#' linear trajectory has standard deviation `t^(3/2) / sqrt(3)` on the
#' standardized domain. Its maximum over `[0, 1]` is therefore the
#' grid-independent reference scale `sigma_ref = 1 / sqrt(3)`. The selected
#' global prior is assigned a scale `s_gamma` such that
#' `Pr(gamma / sqrt(3) > global_upper) = global_alpha`. The fitted object
#' retains `sigma_ref` and `s_gamma` as `data$reference_sd` and
#' `fit_spec$global_scale`.
#'
#' The non-centered parameterization represents each state transition with
#' standard-normal innovations. It is usually the better starting point when
#' observations are sparse or noisy, or when most adaptive local scales are
#' strongly shrunk. The centered parameterization samples `f` and `fprime`
#' directly under the same transition density. It can be more efficient when
#' the latent states and large local changes are strongly identified by the
#' likelihood. Neither parameterization dominates in every dataset. Fit both
#' when practical and compare divergences, R-hat, effective sample size per
#' second, leapfrog counts, and treedepth hits with [get_diagnostics()].
#'
#' With `adaptive_prior = "regularized_horseshoe"`, the effective local scale
#' is
#' \deqn{\tilde\xi_j^2 =
#' \frac{c^2\xi_j^2}{c^2+\gamma^2\xi_j^2},}
#' where \eqn{c = \code{slab_scale}\sqrt{c_{\mathrm{aux}}}} and
#' \eqn{c_{\mathrm{aux}}\sim\mathrm{InvGamma}(\nu/2,\nu/2)} with
#' \eqn{\nu=\code{slab_df}}. Thus the effective interval diffusion standard
#' deviation \eqn{\gamma\tilde\xi_j} approaches the finite slab scale for
#' extreme local scales. This can stabilize isolated, strongly identified
#' changes while preserving horseshoe shrinkage near zero.
#'
#' The adaptive transition over an interval of length \eqn{\Delta_j} is exactly
#' bivariate normal with mean
#' \eqn{(f'_j, f_j+\Delta_j f'_j)} and covariance
#' \eqn{\gamma^2\xi_j^2 Q(\Delta_j)}, where
#' \eqn{Q(\Delta)=((\Delta,\Delta^2/2),(\Delta^2/2,\Delta^3/3))}.
#'
#' @return An object of class `ibmfit`, or Stan code when `get_code = TRUE`.
#' @export
ibm <- function(t = NULL, y = NULL, infer_at = NULL, adaptive = FALSE,
                family = c("gaussian", "student_t", "bernoulli", "binomial",
                           "poisson", "negative_binomial", "lognormal",
                           "gamma", "beta", "exponential"),
                trials = NULL, exposure = 1, student_df = 4,
                parameterization = c("noncentered", "centered"),
                smoothing_prior = "tau ~ lognormal(-2, 0.5);",
                global_prior = c("half_cauchy", "half_normal"),
                global_upper = 1, global_alpha = 0.05,
                adaptive_prior = c("horseshoe", "regularized_horseshoe"),
                slab_scale = 1, slab_df = 4,
                regularized_retry = c("ask", "never"),
                log_sigma = list(mu = -1, sd = 1),
                log_phi = list(mu = 0, sd = 1),
                initial_sd = 5,
                iter = 2000, chains = 4,
                cores = getOption("mc.cores", chains),
                init = "random",
                max_treedepth = 12, adapt_delta = 0.9,
                get_code = FALSE, print_code = FALSE,
                stan_file = NULL, ...) {
  .ibm_fit(
    t = t, y = y, infer_at = infer_at, adaptive = adaptive,
    family = family, trials = trials, exposure = exposure,
    student_df = student_df,
    parameterization = parameterization,
    smoothing_prior = smoothing_prior, global_prior = global_prior,
    global_upper = global_upper, global_alpha = global_alpha,
    adaptive_prior = adaptive_prior, slab_scale = slab_scale,
    slab_df = slab_df, regularized_retry = regularized_retry,
    log_sigma = log_sigma, log_phi = log_phi, initial_sd = initial_sd,
    iter = iter, chains = chains, cores = cores, init = init,
    max_treedepth = max_treedepth, adapt_delta = adapt_delta,
    get_code = get_code, print_code = print_code, stan_file = stan_file,
    prior_only = FALSE, ...
  )
}
