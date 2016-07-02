// $$
// \begin{aligned}
// \frac{\mathrm{d}N_1}{\mathrm{d}t} & = & \rho_1 N_1 \bigg(1 - \frac{N_1}{K_1}\bigg) - c_1 N_1 N_2 \\
// \frac{\mathrm{d}N_2}{\mathrm{d}t} & = & -\rho_2 N_2 \bigg(1 + \frac{N_2}{K_2}\bigg) + c_2 N_1 N_2 \\
// \mathrm{d} \rho_1 & = & \rho_1 \sigma_{\rho_1} \mathrm{d}W_t
// \end{aligned}
//
// $$
// theta[1] = ln_alpha
// theta[2] = k1
// theta[3] = b
// theta[4] = d
// theta[5] = k2
// theta[6] = c

functions {
  real[] pp(real   t,
            real[] y,
            real[] theta,
            real[] x_r,
            int[]  x_i) {
    real dydt[2];
    dydt[1] <-  exp(theta[1]) * y[1] * (1 - y[1] / theta[2]) - theta[3] * y[1] * y[2];
    dydt[2] <- -theta[4] * y[2] * (1 + y[2] / theta[5]) + theta[6] * y[1] * y[2];
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
  real<lower=0> mu; // \log{\mu} is the mean log birth rate of hares
  real<lower=0> sigma; // Standard deviation of hare birth rate
}
