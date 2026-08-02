data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T - 1] deltat;
  // OBSERVATION_DATA
  real<lower=0> global_scale;
  int<lower=1, upper=2> global_prior;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
  real<lower=0> initial_sd;
  int<lower=0, upper=1> prior_only;
}
transformed data {
  vector[T - 1] sqrt_deltat = sqrt(deltat);
  vector[T - 1] deltat_3_2 = deltat .* sqrt_deltat;
}
parameters {
  // OBSERVATION_PARAMETERS
  real<lower=0, upper=1> gamma_unif;
  vector<lower=0, upper=1>[T - 1] xi_unif;
  real<lower=0> slab_aux;
  vector[2] z_initial;
  array[T - 1] vector[2] z_transition;
}
transformed parameters {
  // OBSERVATION_TRANSFORMED_PARAMETERS
  real<lower=0> gamma = global_prior == 1
    ? global_scale * tan(0.5 * pi() * gamma_unif)
    : global_scale * inv_Phi(0.5 + 0.5 * gamma_unif);
  vector<lower=0>[T - 1] xi = tan(0.5 * pi() * xi_unif);
  real<lower=0> slab = slab_scale * sqrt(slab_aux);
  vector<lower=0>[T - 1] xi_regularized =
    sqrt(square(slab) * square(xi)
         ./ (square(slab) + square(gamma) * square(xi)));
  vector[T] fprime;
  vector[T] f;

  fprime[1] = initial_sd * z_initial[1];
  f[1] = initial_sd * z_initial[2];
  for (i in 2:T) {
    real h = deltat[i - 1];
    real tau_i = gamma * xi_regularized[i - 1];
    fprime[i] = fprime[i - 1]
      + tau_i * sqrt_deltat[i - 1] * z_transition[i - 1][1];
    f[i] = f[i - 1] + h * fprime[i - 1]
      + tau_i * deltat_3_2[i - 1]
        * (0.5 * z_transition[i - 1][1]
           + inv_sqrt(12.0) * z_transition[i - 1][2]);
  }
}
model {
  // OBSERVATION_PRIORS
  gamma_unif ~ uniform(0, 1);
  xi_unif ~ uniform(0, 1);
  slab_aux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  z_initial ~ std_normal();
  for (i in 1:(T - 1)) z_transition[i] ~ std_normal();
  if (!prior_only) {
    // OBSERVATION_LIKELIHOOD
  }
}
generated quantities {
  vector<lower=0>[T - 1] lambda_interval =
    square(gamma * xi_regularized);
}
