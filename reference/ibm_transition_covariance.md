# Integrated Brownian motion transition matrices

Integrated Brownian motion transition matrices

## Usage

``` r
ibm_transition_covariance(delta, lambda = 1)

ibm_transition_cholesky(delta, lambda = 1)
```

## Arguments

- delta:

  Positive interval length.

- lambda:

  Positive derivative diffusion variance rate.

## Value

A 2 by 2 covariance matrix or its lower Cholesky factor, ordered as
`(fprime, f)`.
