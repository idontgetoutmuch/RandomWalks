functions {
  real[] pendulum(real t,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i) {
    real dydt[2];
    dydt[1] <- y[2];
    dydt[2] <- - theta[1] * sin(y[1]);
    return dydt;
  }
}
data {
  int<lower=1> T;
  real y0[2];
  real z[T];
  real t0;
  real ts[T];
}
transformed data {
  real x_r[0];
  int x_i[0];
}
parameters {
  real theta[1];
  vector<lower=0>[1] sigma;
}
model {
  real y_hat[T,2];
  real z_hat[T];
  theta ~ normal(0,1);
  sigma ~ cauchy(0,2.5);
  y_hat <- integrate_ode(pendulum, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T) {
    z_hat[t] <- sin(y_hat[t,1]);
    z[t] ~ normal(z_hat[t], sigma);
  }
}
