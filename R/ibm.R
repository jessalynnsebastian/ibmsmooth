#' Fit an integrated Brownian motion smoother
#'
#' Fits either ordinary integrated Brownian motion (IBM) with one diffusion
#' scale or the locally adaptive IBM described by
#' \eqn{df'(t)=\gamma\xi(t)dW(t)}, \eqn{df(t)=f'(t)dt}.  In the adaptive model,
#' each interval has a half-Cauchy local scale \eqn{\xi_j}.
#'
#' @param t Numeric observation locations.
#' @param y Numeric observations. Replicates at the same location are allowed.
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
#'   prior through `Pr(gamma * sigma_ref > U) = global_alpha`.
#' @param global_alpha Tail probability used in the global-prior calibration.
#' @param adaptive_prior Adaptive local-scale prior. `"horseshoe"` is the
#'   default; `"regularized_horseshoe"` adds a finite Student-t slab.
#' @param slab_scale Positive slab scale for the regularized horseshoe, on the
#'   internally standardized derivative-diffusion scale.
#' @param slab_df Positive slab degrees of freedom. The default is 4.
#' @param regularized_retry Either `"ask"` (the default) or `"never"`.
#'   After a problematic interactive ordinary-horseshoe fit, `"ask"` offers
#'   to rerun the same model and HMC settings with a regularized horseshoe.
#'   Non-interactive sessions never start a refit automatically.
#' @param log_sigma List with `mu` and positive `sd` for the log observation
#'   standard deviation on the internally standardized response scale.
#' @param initial_sd Positive prior standard deviation for the initial function
#'   value and derivative on the internally standardized scale.
#' @param iter,chains,cores Stan sampling controls.
#' @param init Initial values passed to [rstan::sampling()].
#' @param max_treedepth,adapt_delta Stan HMC controls.
#' @param get_code If `TRUE`, return the selected Stan program rather than fit.
#' @param ... Additional arguments passed to [rstan::stan()].
#'
#' @details
#' Time is mapped to the unit interval and the response is centered and scaled
#' to unit standard deviation before fitting. Consequently, `smoothing_prior`
#' is expressed
#' on that stable internal scale. Natural-scale draws are returned by the
#' extraction helpers. Locations in `infer_at` are evaluated after fitting
#' using exact conditional IBM bridges, avoiding additional Stan parameters.
#'
#' For an adaptive fit, `sigma_ref` is the geometric mean of the marginal
#' standard deviations of the unit-diffusion IBM departure from its initial
#' linear trajectory. Conditioning makes the stochastic variance zero at the
#' first location, so the geometric mean is taken over subsequent locations.
#' The selected global prior is assigned a scale `s_gamma` such that
#' `Pr(gamma * sigma_ref > global_upper) = global_alpha`. The fitted object
#' retains these as `data$reference_sd` and `fit_spec$global_scale`.
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
                parameterization = c("noncentered", "centered"),
                smoothing_prior = "tau ~ lognormal(-2, 0.5);",
                global_prior = c("half_cauchy", "half_normal"),
                global_upper = 1, global_alpha = 0.05,
                adaptive_prior = c("horseshoe", "regularized_horseshoe"),
                slab_scale = 1, slab_df = 4,
                regularized_retry = c("ask", "never"),
                log_sigma = list(mu = -1, sd = 1),
                initial_sd = 5,
                iter = 2000, chains = 4,
                cores = getOption("mc.cores", chains),
                init = "random",
                max_treedepth = 12, adapt_delta = 0.9,
                get_code = FALSE, ...) {
  .ibm_fit(
    t = t, y = y, infer_at = infer_at, adaptive = adaptive,
    parameterization = parameterization,
    smoothing_prior = smoothing_prior, global_prior = global_prior,
    global_upper = global_upper, global_alpha = global_alpha,
    adaptive_prior = adaptive_prior, slab_scale = slab_scale,
    slab_df = slab_df, regularized_retry = regularized_retry,
    log_sigma = log_sigma, initial_sd = initial_sd,
    iter = iter, chains = chains, cores = cores, init = init,
    max_treedepth = max_treedepth, adapt_delta = adapt_delta,
    get_code = get_code, prior_only = FALSE, ...
  )
}
