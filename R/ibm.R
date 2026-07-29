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
#' @param global_scale Positive half-Cauchy scale for `gamma` in the adaptive
#'   global--local prior.
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
#' Time is shifted and divided by its mean grid spacing and the response is
#' standardized before fitting. Consequently, `smoothing_prior` is expressed
#' on that stable internal scale. Natural-scale draws are returned by the
#' extraction helpers. Locations in `infer_at` are evaluated after fitting
#' using exact conditional IBM bridges, avoiding additional Stan parameters.
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
                global_scale = 0.1,
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
    smoothing_prior = smoothing_prior, global_scale = global_scale,
    log_sigma = log_sigma, initial_sd = initial_sd,
    iter = iter, chains = chains, cores = cores, init = init,
    max_treedepth = max_treedepth, adapt_delta = adapt_delta,
    get_code = get_code, prior_only = FALSE, ...
  )
}
