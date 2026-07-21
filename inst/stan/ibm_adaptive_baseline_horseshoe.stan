data {
  int<lower=1> N_obs;
  int<lower=2> T;
  array[N_obs] int<lower=1, upper=T> obs_time_idx;
  vector<lower=0>[T-1] deltat;
  vector[N_obs] y_obs;
  int<lower=0, upper=1> regular;
  real log_sigma_mu;
  real<lower=0> log_sigma_sd;
  real log_tau_mu;
  real<lower=0> log_tau_sd;
  real<lower=0> zeta;
  real<lower=0> initial_sd;
}

transformed data {
  array[T-1] matrix[2, 2] Lt;
  matrix[2, 2] L_initial = initial_sd * to_matrix(
    [[1.0, 0.0], [0.5, inv_sqrt(12.0)]]);

  Lt[1] = to_matrix(
    [[sqrt(deltat[1]), 0],
     [0.5 * deltat[1]^1.5, deltat[1]^1.5 * inv_sqrt(12.0)]]);
  if (T > 2) {
    for (i in 2:(T-1)) {
      if (regular == 0) {
        Lt[i] = to_matrix(
          [[sqrt(deltat[i]), 0],
           [0.5 * deltat[i]^1.5, deltat[i]^1.5 * inv_sqrt(12.0)]]);
      } else {
        Lt[i] = Lt[1];
      }
    }
  }
}

parameters {
  real log_sigma;
  real log_tau0;
  vector<lower=0, upper=1>[T-1] z_lambda;
  real<lower=0, upper=1> z_gamma;
  matrix[2, T] z;
}

transformed parameters {
  real sigma = exp(log_sigma_mu + log_sigma * log_sigma_sd);
  real tau0 = exp(log_tau_mu + log_tau0 * log_tau_sd);
  real gamma = zeta * tan(z_gamma * pi() / 2);
  vector[T-1] lambda;
  vector[T-1] psi;
  vector[T-1] local_variance;
  vector[T-1] tau;
  vector[T] f;
  vector[T] fprime;

  for (i in 1:(T-1)) {
    lambda[i] = tan(z_lambda[i] * pi() / 2);
    psi[i] = gamma * lambda[i];
    local_variance[i] = square(tau0) + square(psi[i]);
    tau[i] = sqrt(local_variance[i]);
  }

  {
    vector[2] initial_state = L_initial * z[, 1];
    fprime[1] = initial_state[1];
    f[1] = initial_state[2];
  }
  for (i in 2:T) {
    vector[2] state =
      [fprime[i - 1], f[i - 1] + deltat[i - 1] * fprime[i - 1]]'
      + tau[i - 1] * Lt[i - 1] * z[, i];
    fprime[i] = state[1];
    f[i] = state[2];
  }
}

model {
  log_sigma ~ std_normal();
  log_tau0 ~ std_normal();
  z_lambda ~ uniform(0, 1);
  z_gamma ~ uniform(0, 1);
  to_vector(z) ~ std_normal();

  for (n in 1:N_obs) {
    y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
  }
}
