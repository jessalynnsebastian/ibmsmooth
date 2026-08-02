data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T - 1] deltat;
  // OBSERVATION_DATA
  real<lower=0> initial_sd;
  int<lower=0, upper=1> prior_only;
}
parameters {
  // OBSERVATION_PARAMETERS
  real<lower=0> tau;
  vector[2] z_initial;
  array[T - 1] vector[2] z_transition;
}
transformed parameters {
  // OBSERVATION_TRANSFORMED_PARAMETERS
  vector[T] fprime;
  vector[T] f;

  fprime[1] = initial_sd * z_initial[1];
  f[1] = initial_sd * z_initial[2];
  for (i in 2:T) {
    real h = deltat[i - 1];
    real h32 = h * sqrt(h);
    fprime[i] = fprime[i - 1]
      + tau * sqrt(h) * z_transition[i - 1][1];
    f[i] = f[i - 1] + h * fprime[i - 1]
      + tau * h32
        * (0.5 * z_transition[i - 1][1]
           + inv_sqrt(12.0) * z_transition[i - 1][2]);
  }
}
model {
  // OBSERVATION_PRIORS
  z_initial ~ std_normal();
  for (i in 1:(T - 1)) z_transition[i] ~ std_normal();
  // SMOOTHING_PRIOR
  if (!prior_only) {
    // OBSERVATION_LIKELIHOOD
  }
}
