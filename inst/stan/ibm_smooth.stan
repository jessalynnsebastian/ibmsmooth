data {
  int<lower=1> N_obs; // number of observations (including replicates if there are any)
  int<lower=1> T; // number of latent grid points (unique times + infer_at)
  array[N_obs] int<lower=1, upper=T> obs_time_idx;  // maps each obs to a grid index
  vector[T-1] deltat; // time diffs on the (scaled) grid
  vector[N_obs] y_obs; // observed y values (also scaled)
  int<lower=0, upper=1> regular;

  // gaussian lik - lognormal prior for sigma
  real log_sigma_mu;
  real<lower=0> log_sigma_sd;

  // ibm - lognormal prior for tau
  real log_tau_mu;
  real<lower=0> log_tau_sd;
}

transformed data {
  array[T-1] matrix[2,2] Lt;

  Lt[1] = to_matrix(
    [[sqrt(deltat[1]), 0],
     [(deltat[1]^1.5) * 0.5, (deltat[1]^1.5) * inv_sqrt(12)]]
  );

  for (i in 2:(T-1)) {
    if (regular == 0) {
      Lt[i] = to_matrix(
        [[sqrt(deltat[i]), 0],
         [(deltat[i]^1.5) * 0.5, (deltat[i]^1.5) * inv_sqrt(12)]]
      );
    } else {
      Lt[i] = Lt[1];
    }
  }
}

parameters {
  real log_sigma;
  real log_tau;
  matrix[2, T] z;
}

transformed parameters {
  vector[T] f;
  vector[T] fprime;

  real tau = exp(log_tau_mu + log_tau * log_tau_sd);
  real sigma = exp(log_sigma_mu + log_sigma * log_sigma_sd);

  // transform z to f, fprime
  f[1] = 5 * z[2, 1];
  fprime[1] = tau^2 * z[1, 1];

  for (i in 2:T) {
    vector[2] state =
      [ fprime[i - 1],
        f[i - 1] + deltat[i - 1] * fprime[i - 1] ]'
      + tau * Lt[i - 1] * z[, i];

    fprime[i] = state[1];
    f[i] = state[2];
  }
}

model {
  // hyperpriors
  log_sigma ~ std_normal();
  log_tau ~ std_normal();

  // ibm noncentered
  to_vector(z) ~ std_normal();

  // likelihood (replicates just repeat obs_time_idx)
  for (n in 1:N_obs) {
    y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
  }
}
