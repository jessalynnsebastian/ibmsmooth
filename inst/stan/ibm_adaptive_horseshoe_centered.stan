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

  // ibm - global scale for tau/gamma
  real zeta;
  real<lower=0> initial_sd;
}

transformed data{
  array[T-1] matrix[2,2] Lt;
  matrix[2,2] L_initial = initial_sd * to_matrix(
    [[1.0, 0.0], [0.5, inv_sqrt(12.0)]]);

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
  real log_sigma;
  vector[T] f;
  vector[T] fprime;
  vector<lower=0, upper=1>[T-1] z_tau;
  real<lower=0, upper=1> z_gamma;
}

transformed parameters {
  real sigma = exp(log_sigma_mu + log_sigma * log_sigma_sd);
  vector[T-1] tau;
  real gamma = zeta * tan(z_gamma * pi() / 2); // half-Cauchy(0, zeta)
  for(i in 1:(T-1)) {
    tau[i] = gamma * tan(z_tau[i] * pi() / 2);
  }
}

model {
  // local shrinkage parameters
  z_tau ~ uniform(0, 1);
  z_gamma ~ uniform(0, 1);

  // hyperpriors
  log_sigma ~ std_normal();

  // ibm prior (centered)
  to_vector({fprime[1], f[1]}) ~
    multi_normal_cholesky(rep_vector(0, 2), L_initial);
  for (i in 2:T) {
    vector[2] mu;
    mu[1] = fprime[i - 1];
    mu[2] = f[i - 1] + deltat[i - 1] * fprime[i - 1];
    to_vector({fprime[i], f[i]}) ~ multi_normal_cholesky(mu, tau[i - 1] * Lt[i - 1]);
  }

  // likelihood
  for (n in 1:N_obs) {
    y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
  }
}
