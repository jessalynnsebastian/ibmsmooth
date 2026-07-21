data {
  int<lower=1> N_obs; // number of observations (including replicates if there are any)
  int<lower=1> T; // number of latent grid points (unique times + infer_at)
  array[N_obs] int<lower=1, upper=T> obs_time_idx;  // maps each obs to a grid index
  vector[T-1] deltat; // time diffs on the (scaled) grid
  vector[N_obs] y_obs; // observed y values (also scaled)
  int<lower=0, upper=1> regular;

  // gaussian likelihood - lognormal prior for sigma
  real log_sigma_mu;
  real<lower=0> log_sigma_sd;

  // ibm - lognormal prior for tau
  real log_tau_mu;
  real<lower=0> log_tau_sd;
}

transformed data{
  array[T-1] matrix[2,2] Lt;

  Lt[1] = to_matrix( // cholesky decomp of ibm covariance
    [[sqrt(deltat[1]), 0],
    [(deltat[1]^1.5) * 0.5, (deltat[1]^1.5) * inv_sqrt(12)]]);

  for(i in 2:(T-1)) {
    if (regular == 0) {
      Lt[i] = to_matrix( // cholesky decomp of ibm covariance
        [[sqrt(deltat[i]), 0],
        [(deltat[i]^1.5) * 0.5, (deltat[i]^1.5) * inv_sqrt(12)]]);
    } else {
      Lt[i] = Lt[1];
    }
  }
}

parameters {
  vector[T] log_tau;
  real log_sigma;
  vector[T] f;
  vector[T] fprime;
}

transformed parameters {
  real sigma = exp(log_sigma_mu + log_sigma * log_sigma_sd);
  vector[T-1] tau;
  tau[1] = exp(log_tau_mu + log_tau[1] * log_tau_sd);
  for(i in 2:(T-1)) {
    tau[i] = tau[i-1] * exp(log_tau[i] * log_tau_sd);
  }
}

model {
  // hyperpriors
  log_sigma ~ std_normal();
  log_tau ~ std_normal();

  // ibm prior
  fprime[1] ~ std_normal();
  f[1] ~ std_normal();
  for (i in 2:T) {
    [fprime[i], f[i]]' ~ multi_normal_cholesky(
      [fprime[i - 1], f[i - 1] + deltat[i - 1] * fprime[i - 1]]',
      tau[i - 1] * Lt[i - 1]);
  }

  // likelihood
  for (n in 1:N_obs) {
    y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
  }
}
