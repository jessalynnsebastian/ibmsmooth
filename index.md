# ibmsmooth

`ibmsmooth` fits integrated Brownian motion smoothers for Gaussian
observations using Stan. The package intentionally contains two models:

- ordinary IBM with one positive smoothing parameter and a user-supplied
  Stan prior for that parameter; and
- locally adaptive IBM with horseshoe global–local shrinkage on the
  derivative diffusion rate.

For locations (t_j\<t\_{j+1}), (*j=t*{j+1}-t_j), and state
(x_j=(f’(t_j),f(t_j))^T), both models use the exact transition

\[ x\_{j+1}x_j,\_j N_2. \]

Ordinary IBM sets (\_j=^2). The adaptive model sets (\_j=^(2_j)2), with
(\_j^+(0,1)) and (^+(0,)).

## Installation

``` r

remotes::install_github("jessalynnsebastian/ibmsmooth")
```

## Examples

``` r

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

plot(fit_adaptive)
plot_diffusion(fit_adaptive)
```

The data are standardized internally, so custom smoothing priors are
specified on the standardized model scale. Extraction functions
transform posterior draws back to the original response and time scales.

## License

This project is released under the MIT license. See
[LICENSE](https://jessalynnsebastian.github.io/ibmsmooth/LICENSE).
