model LV {

  const h         = 0.02;   // time step
  const delta_abs = 1.0e-2; // absolute error tolerance
  const delta_rel = 1.0e-5; // relative error tolerance
  const epsilon   = 2.0e-2; // observation error tolerance

  dim s(2); // number of species
  dim n(2); // number of parameters for each species process

  param mu[s,n];      // mean of birth and death rates
  param sigma[s,n];   // volatility of birth and death rates
  noise w[s,n];       // noise terms
  state ln_beta[s,n]; // birth and death rates
  state H, L;         // hares, lynxes
  obs y_H, y_L;       // observations

  // Reasonable parameters are 0.5, 0.02, 0.4 and 0.01.
  sub parameter {
    mu[0,0] ~ uniform(0.4, 0.6);
    sigma[0,0] ~ uniform(0.0, 0.1);
    mu[0,1] ~ uniform(0.01, 0.03);
    sigma[0,1] ~ uniform(0.0, 0.01);
    mu[1,0] ~ uniform(0.3, 0.5);
    sigma[1,0] ~ uniform(0.0, 0.1);
    mu[1,1] ~ uniform(0.009, 0.011);
    sigma[1,1] ~ uniform(0.0, 0.001);
  }

  sub proposal_parameter {
    mu[s,n] ~ truncated_gaussian(mu[s,n], 0.02, 0.0, 1.0);
    sigma[s,n] ~ truncated_gaussian(sigma[s,n], 0.01, 0.0, 0.5);
  }

  // Starting populations of 50 hares and 50 lynxes.
  sub initial {
    param beta[s,n];
    H ~ log_normal(log(50.0), 0.5);
    L ~ log_normal(log(50.0), 0.5);
    beta[s,n] ~ gaussian(mu[s,n], sigma[s,n]);
    ln_beta <- log(beta);
  }

  sub transition(delta = h) {
    w[s,n] ~ normal(0.0, sqrt(h));
    ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg = 'RK4(3)') {
      dH/dt = exp(ln_beta[0,0]) * H - exp(ln_beta[0,1]) * H * L;
      dL/dt = exp(ln_beta[1,1]) * H * L - exp(ln_beta[1,0]) * L;
      dln_beta[s,n]/dt = -sigma[s,n] * sigma[s,n] / 2 - sigma[s,n] * w[s,n] / h;
    }
  }

  sub observation {
    y_H ~ log_normal(log(H), 0.2);
    y_L ~ log_normal(log(L), 0.2);
  }
}
