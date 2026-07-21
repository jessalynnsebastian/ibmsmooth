# Extract posterior hyperparameter samples

Extract posterior samples for model hyperparameters. By default,
returned quantities are transformed to the original response and time
scales.

## Usage

``` r
get_hyperparameter_samples(
  ibmfit,
  n_samples = 1000,
  format = c("long", "wide"),
  natural = TRUE,
  include_internal = FALSE
)
```

## Arguments

- ibmfit:

  An object of class `ibmfit`.

- n_samples:

  Number of samples to draw for INLA fits, or maximum number of Stan
  draws to keep. Use `NULL` to keep all Stan draws.

- format:

  Output format, either `"long"` or `"wide"`.

- natural:

  Logical. If `TRUE`, return quantities transformed to the original
  response/time scale whenever possible. If `FALSE`, return internal
  backend hyperparameters.

- include_internal:

  Logical. If `TRUE` and `natural = TRUE`, also include internal backend
  hyperparameters.

## Value

A data frame of hyperparameter samples.

## Details

For INLA fits, hyperparameters are sampled with `intern = TRUE`, so the
Gaussian likelihood hyperparameter is log precision and the rgeneric
field hyperparameters are log diffusion precisions on the scaled model
scale. On the natural scale this function returns the implied IBM
process standard deviation and precision. For adaptive INLA fits these
are returned over transition midpoints, and the global baseline process
standard deviation and precision are also returned as
`global_process_sd` and `global_process_precision`.

For Stan fits, recognized transformed parameters such as `sigma`, `tau`,
and `gamma` are converted from the scaled model units to the original
data units.
