# ibmsmooth

`ibmsmooth` is an R package for Bayesian smoothing with integrated Brownian
motion models fitted in Stan. It supports ordinary and locally adaptive
smoothing, irregularly spaced observations, and inference for both a function
and its derivative.

The package includes several common observation models, exact continuous-time
prediction, posterior diagnostics, plotting, and simulation-performance
summaries.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("jessalynnsebastian/ibmsmooth")
```

The package uses `rstan`, so a working Stan toolchain is required. See the
[RStan getting started guide](https://mc-stan.org/rstan/articles/rstan.html) if
you need help configuring it.

## Documentation

Usage examples and model details are kept in the package vignettes:

- [Ordinary and locally adaptive IBM smoothing](https://jessalynnsebastian.github.io/ibmsmooth/articles/models.html)
- [Simulating integrated Brownian motion](https://jessalynnsebastian.github.io/ibmsmooth/articles/simulating-ibm.html)
- [Function reference](https://jessalynnsebastian.github.io/ibmsmooth/reference/)

## License

This project is released under the [MIT License](LICENSE).
