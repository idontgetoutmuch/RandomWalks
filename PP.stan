// Infer growth rate for hares

// This is the model in mathematical notation.
//
// $$
// \begin{aligned}
// \frac{\mathrm{d}N_1}{\mathrm{d}t} & = & \rho_1 N_1 \bigg(1 - \frac{N_1}{K_1}\bigg) - c_1 N_1 N_2 \\
// \frac{\mathrm{d}N_2}{\mathrm{d}t} & = & -\rho_2 N_2 \bigg(1 + \frac{N_2}{K_2}\bigg) + c_2 N_1 N_2 \\
// \mathrm{d} \rho_1 & = & \rho_1 \sigma_{\rho_1} \mathrm{d}W_t
// \end{aligned}
//
// $$
// theta[1] = sigma
// theta[2] = k1
// theta[3] = b
// theta[4] = d
// theta[5] = k2
// theta[6] = c
// theta[7] = w
// theta[8] = h

functions {
  real[] pp(real   t,
            real[] y,
            real[] theta,
            real[] x_r,
            int[]  x_i) {
    real dydt[2];
    dydt[1] <-  exp(y[3]) * y[1] * (1 - y[1] / theta[2]) - theta[3] * y[1] * y[2];
    dydt[2] <- -theta[4] * y[2] * (1 + y[2] / theta[5]) + theta[6] * y[1] * y[2];
    dydt[3] <- -theta[1] * theta[1] * 0.5 - theta[1] * theta[7] / theta[8];
    return dydt;
  }
}

data {
  int<lower=1> N_t; // Number of time points
  real y[N_t, 1];   // Observations
  real t0;          // Start time
  real k1;          // Hare carrying capacity
  real b;           // Hare death rate per lynx
  real d;           // Lynx death rate
  real k2;          // Lynx carrying capacity
  real c;           // Lynx birth rate per hare
  real p;           // Initial hares
  real z;           // Initial lynxes
  real h;           // Time step
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
  real<lower=0> mu;    // \log{\mu} is the mean log birth rate of hares
  real<lower=0> sigma; // Standard deviation of hare birth rate
}

transformed parameters{
  theta[2] = k1;
  theta[3] = b;
  theta[4] = d;
  theta[5] = k2;
  theta[6] = c;
  theta[8] = h;
}

model {
  real w;
  real y1[3];
  real y_hat[T,3];

  // Priors
  y0[1] ~ lognormal(log(100.0), 0.2)
  y0[2] ~ lognormal(log(50.0), 0.1)
  mu    ~ uniform(0.0, 1.0)
  sigma ~ uniform(0.0, 0.5)
  y0[3] ~ normal(log(mu), sigma)

  for (t in 1:T){
    w ~ normal(0, sqrt(h));
    theta[1] = sigma;
    theta[7] = w;
    y1 = integrate_ode_cvode(pp, y0, t0, ts[t:t], theta, x_r, x_i, rel_tol, abs_tol, max_num_steps);
    y_hat[t] = y1;
    t0 = ts[t];
    y0 = y1;
    y[t] ~ lognormal(log(y_hat[t,1]), 0.1);
  }
}
