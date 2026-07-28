# ibmsmooth

`ibmsmooth` fits integrated Brownian motion smoothers for Gaussian observations
using Stan. The package intentionally contains two models:

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
\(\gamma\sim\mathrm{C}^+(0,\texttt{global_scale})\).

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
  global_scale = 0.1
)

# Exact derivative-marginalized implementation of the same model
fit_adaptive_fast <- ibm(
  t, y,
  adaptive = TRUE,
  fast = TRUE,
  global_scale = 0.1
)

# Verify equality of the two conditional IBM priors
verify_ibm_equivalence(
  sort(unique(t)),
  lambda = rep(0.1, length(unique(t)) - 1)
)

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

The data are standardized internally, so custom smoothing priors are specified
on the standardized model scale. Extraction functions transform posterior
draws back to the original response and time scales. Locations supplied through
`infer_at` are not added to the Stan state vector; `predict_curve()` samples
them afterward using exact conditional IBM bridges.

Setting `fast = TRUE` selects separate Stan programs that analytically filter
out derivative states during HMC and draw their joint posterior trajectory
afterward. The observation-level latent function remains in the model, so this
representation can also support non-Gaussian likelihoods. The standard Stan
programs remain the default; the alternative programs are the files ending in
`_fast.stan`.

## License

This project is released under the MIT license. See [LICENSE](LICENSE).
