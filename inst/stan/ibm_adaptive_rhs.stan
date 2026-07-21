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
  // regularized horseshoe - slab scale
  real<lower=0> slab_scale;
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
  real log_sigma;
  vector<lower=0, upper=1>[T-1] z_tau;
  real<lower=0, upper=1> z_gamma;
  matrix[2, T] z;
}

transformed parameters {
  vector[T] f;
  vector[T] fprime;
  vector[T-1] tau;
  real gamma = zeta * tan(z_gamma * pi() / 2); // half-Cauchy(0, zeta)
  real sigma = exp(log_sigma_mu + log_sigma * log_sigma_sd);

  vector[T-1] tau_hc; // raw local half-Cauchy(0,1)
  vector[T-1] tau_tilde; // regularized (slab) local scales
  {
    real c2 = square(slab_scale);
    real gamma2 = square(gamma);

    for (i in 1:(T-1)) {
      // raw half-Cauchy(0,1)
      tau_hc[i] = tan(z_tau[i] * pi() / 2);

      // regularized horseshoe effective local scale
      tau_tilde[i] = sqrt(
        (c2 * square(tau_hc[i])) /
        (c2 + gamma2 * square(tau_hc[i]))
      );

      // final local diffusion per interval
      tau[i] = gamma * tau_tilde[i];
    }
  }

  // transform z to f, fprime
  f[1] = z[2, 1];
  fprime[1] = z[1, 1];
  for (i in 2:T) {
    vector[2] state = [fprime[i - 1], f[i - 1] + deltat[i - 1] * fprime[i - 1]]' + tau[i - 1] * Lt[i - 1] * z[, i];
    fprime[i] = state[1];
    f[i] = state[2];
  }
}

model {
  // hyperpriors
  log_sigma ~ std_normal();
  z_tau ~ uniform(0, 1);
  z_gamma ~ uniform(0, 1);

  // ibm prior
  z[1, ] ~ std_normal();
  z[2, ] ~ std_normal();

  // likelihood
  for (n in 1:N_obs) {
    y_obs[n] ~ normal(f[obs_time_idx[n]], sigma);
  }
}
