# ibmsmooth

Note that this is a _work in progress_!

`ibmsmooth` fits global and locally adaptive integrated Brownian motion smoothers for Gaussian observations using either a Stan backend or an INLA backend.

The package is designed for flexible smoothing workflows with:

- global and locally adaptive smoothness models
- posterior summaries and plotting helpers
- support for both MCMC-based and latent Gaussian inference

## Installation

### From GitHub

```r
remotes::install_github("jessalynsebastian/ibmsmooth")
```

### Backend requirements

- Stan backend: `rstan` and a working C++ toolchain
- INLA backend: the `INLA` package

## Quick example

```r
library(ibmsmooth)

set.seed(1)
t <- sort(runif(40, 0, 10))
f_true <- sin(t)
y <- f_true + rnorm(length(t), 0, 0.25)

fit <- ibm(
  t = t,
  y = y,
  method = "stan",
  adaptive = FALSE,
  chains = 2,
  cores = 2,
  iter = 1000
)

plot(fit)
summary(fit)
```

## License

This project is released under the MIT license. See [LICENSE](LICENSE) for details.
