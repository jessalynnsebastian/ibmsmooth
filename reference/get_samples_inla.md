# Extract posterior samples from an INLA fit

Helper function to draw posterior samples for the latent function or its
derivative from an INLA fit. Samples are drawn using
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html).
The latent field is rescaled back to the original observation scale.

## Usage

``` r
get_samples_inla(ibmfit, param = c("f", "fprime"), n_samples = 1000)
```

## Arguments

- ibmfit:

  An object of class `ibmfit` produced by the INLA implementation.

- param:

  Character, either `"f"` for the latent function or `"fprime"` for its
  derivative.

- n_samples:

  Integer giving the number of posterior samples to draw.

## Value

A matrix of dimension `n_samples` by length of the latent time grid.
