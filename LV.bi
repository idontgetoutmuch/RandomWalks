model LV {

  const h         = 0.02;   // time step
  const delta_abs = 1.0e-2; // absolute error tolerance
  const delta_rel = 1.0e-5; // relative error tolerance
  const epsilon   = 2.0e-2; // observation error tolerance

  dim s(2); // number of species
  dim p(2); // number of parameters for each species process

  param mu[s,p];      // mean of birth and death rates
  param sigma[s,p];   // volatility of birth and death rates
  noise w[s,p];       // noise terms
  state ln_beta[s,p]; // birth and death rates
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
    mu[s,p] ~ truncated_gaussian(mu[s,p], 0.02, 0.0, 1.0);
    sigma[s,p] ~ truncated_gaussian(sigma[s,p], 0.01, 0.0, 0.5);
  }

  // Starting populations of 50 hares and 50 lynxes.
  sub initial {
    param beta[s,p];
    H ~ log_normal(log(50.0), 0.5);
    L ~ log_normal(log(50.0), 0.5);
    beta[s,p] ~ gaussian(mu[s,p], sigma[s,p]);
    ln_beta <- log(beta);
  }

  sub transition(delta = h) {
    w[s,p] ~ normal(0.0, sqrt(h));
    ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg = 'RK4(3)') {
      dH/dt = exp(ln_beta[0,0]) * H - exp(ln_beta[0,1]) * H * L;
      dL/dt = exp(ln_beta[1,1]) * H * L - exp(ln_beta[1,0]) * L;
      dln_beta[s,p]/dt = -sigma[s,p] * sigma[s,p] / 2 - sigma[s,p] * w[s,p] / h;
    }
  }

  sub observation {
    y_H ~ log_normal(log(H), 0.2);
    y_L ~ log_normal(log(L), 0.2);
  }
}
