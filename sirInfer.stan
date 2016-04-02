functions {
  real[] sir(real t,
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {
    real dydt[3];
    dydt[1] <- - theta[1] * y[1] * y[2];
    dydt[2] <-   theta[1] * y[1] * y[2] - theta[2] * y[2];
    dydt[3] <-   theta[2] * y[2];
    return dydt;
  }
}


data {
  int<lower=1> T;
  real ts[T];
  real y0[3];
  real i[T];
}

transformed data {
  real t0;
  real x_r[0];
  int x_i[0];
  t0 <- 0;
}

parameters {
  real<lower=0> theta[2];
  vector<lower=0>[1] sigma;
}

model {
  real y_hat[T,3];
  theta[1] ~ normal(0,1);
  theta[2] ~ normal(0,1);
  sigma ~ cauchy(0,2.5);
  y_hat <- integrate_ode(pendulum, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T) {
    z[t] ~ normal(y_hat[t], sigma);
  }
}
