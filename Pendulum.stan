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
  real t0;
  real ts[T];
  real theta[1];
  real sigma[2];
}
transformed data {
  real x_r[0];
  int x_i[0];
}
model {
}
generated quantities {
  real y_hat[T,2];
  y_hat <- integrate_ode(pendulum, y0, t0, ts, theta, x_r, x_i);
  for (t in 1:T) {
    y_hat[t,1] <- y_hat[t,1] + normal_rng(0,sigma[1]);
    y_hat[t,2] <- y_hat[t,2] + normal_rng(0,sigma[2]);
  }
}
