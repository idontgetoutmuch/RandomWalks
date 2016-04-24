functions {
  real[] lv(real t,
            real[] y,
            real[] param,
            real[] x_r,
            int[] x_i) {
    real dydt[2];
    dydt[1] =  param[1] * y[1] - param[2] * y[1] * y[2];
    dydt[2] = -param[3] * y[2] + param[4] * y[1] * y[2];
    return dydt;
  }
}

data {
  int<lower=1> T;
  real t0;
  real ts[T];
  real y[T,2];
}

transformed data {
  real x_r[0];
  int x_i[0];
  real rel_tol;
  real abs_tol;
  real max_num_steps;
  rel_tol = 1e-10;
  abs_tol = 1e-10;
  max_num_steps = 1e8;
}

parameters {
  real<lower=0> y0[2];
  vector<lower=0>[2] sigma;
  real<lower=0> param[4];
}

model {
  real y_hat[T,2];
  sigma ~ cauchy(0, 2.5);

  // Note that the priors are *very* informed
  y0[1] ~ normal(35.0, 1);
  y0[2] ~ normal(3.9, 1e-1);

  param[1] ~ normal(0.4,4e-1);
  param[2] ~ normal(0.018,8e-2);
  param[3] ~ normal(0.8,4e-1);
  param[4] ~ normal(0.023,8e-2);

  y_hat = integrate_ode_cvode(lv, y0, t0, ts, param, x_r, x_i, rel_tol, abs_tol, max_num_steps);
  for (t in 1:T){
    y[t] ~ normal(y_hat[t], sigma);
  }
}
