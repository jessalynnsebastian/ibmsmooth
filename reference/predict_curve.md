# Predict an IBM curve after fitting

Draws function values and derivatives at arbitrary locations without
adding those locations to the Stan latent state. Within the observed
range, draws use the exact IBM bridge conditional on adjacent sampled
states. Outside the range they use forward or reverse IBM transitions.

## Usage

``` r
predict_curve(ibmfit, new_t = NULL, n_samples = 1000, seed = NULL)
```

## Arguments

- ibmfit:

  An `ibmfit`.

- new_t:

  Finite numeric prediction locations. The default uses the observed
  locations together with the fit's original `infer_at`.

- n_samples:

  Maximum number of posterior draws.

- seed:

  Optional random seed for bridge and extrapolation draws.

## Value

A list with `t`, `f`, and `fprime`. Draws are rows and locations are
columns.
