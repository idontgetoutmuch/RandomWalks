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
  real y[T,2];
  real t0;
  real ts[T];
  vector<lower=0>[2] sigma;
}
transformed data {
  real x_r[0];
  int x_i[0];
}
parameters {
  real theta[1];
}
model {
  real y_hat[T,2];
  theta ~ normal(0,1);
  y_hat <- integrate_ode(pendulum, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T) {
    y[t] ~ normal(y_hat[t], sigma);
  }
}
