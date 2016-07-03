functions {
    real[] lv(real t,real[] y,real[] theta,real[] x_r,int[] x_i) {
    real dydt[3];

    dydt[1] =  exp(y[3]) * y[1] * (1 - y[1] / x_r[1]) - x_r[2] * y[1] * y[2];
    dydt[2] = -x_r[3] * y[2] * (1 + y[2] / x_r[4]) + x_r[5] * y[1] * y[2];
    dydt[3] = -theta[1] * theta[1] * 0.5 - theta[1] * x_r[6] / x_r[7];

    return dydt;
    }
}

data {
  int<lower=1> T;
  int<lower=0> N;
  real t0;
  real ts[T];
  matrix[T,N] y;
  real k1;          // Hare carrying capacity
  real b;           // Hare death rate per lynx
  real d;           // Lynx death rate
  real k2;          // Lynx carrying capacity
  real c;           // Lynx birth rate per hare
  real p;           // Initial hares
  real l;           // Initial lynxes
  real h;           // Time step
}

transformed data {
  real x_r[7];
  int x_i[0];

  x_r[1] = k1;
  x_r[2] = b;
  x_r[3] = d;
  x_r[4] = k2;
  x_r[5] = c;
  x_r[6] = 1; // FIXME
  x_r[7] = h;
}

parameters {
      real<lower=0> theta[1];
      real<lower=0> sigma_y;
      real<lower=0> y0[3];
}

transformed parameters {
  real y_hat[T,3];
  y_hat = integrate_ode_rk45(lv, y0, t0, ts, theta, x_r, x_i);
}

model {
  sigma_y ~ cauchy(0,1);
  theta ~ normal(0,2);
  y0 ~ normal(5,2);
  for (i in 1:N){
    for (t in 1:T) {
          y[t,i] ~ normal(y_hat[t,1],sigma_y);
      }
  }
}
