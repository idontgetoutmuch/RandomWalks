// Infer growth rate for hares

model PP {
  const h         = 0.1;    // time step
  const delta_abs = 1.0e-3; // absolute error tolerance
  const delta_rel = 1.0e-6; // relative error tolerance

  const a  = 5.0e-1 // Hare growth rate - superfluous for inference
                    // but a reminder of what we should expect
  const k1 = 2.0e2  // Hare carrying capacity
  const b  = 2.0e-2 // Hare death rate per lynx
  const d  = 4.0e-1 // Lynx death rate
  const k2 = 2.0e1  // Lynx carrying capacity
  const c  = 4.0e-3 // Lynx birth rate per hare

  state P, Z       // Hares and lynxes
  state ln_alpha   // Hare growth rate - we express it in log form for
                   // consistency with the inference model
  obs P_obs        // Observations of hares
  param mu, sigma  // Mean and standard deviation of hare growth rate
  noise w          // Noise

  sub parameter {
    mu ~ uniform(0.0, 1.0)
    sigma ~ uniform(0.0, 0.5)
  }

  sub proposal_parameter {
     mu ~ truncated_gaussian(mu, 0.02, 0.0, 1.0);
     sigma ~ truncated_gaussian(sigma, 0.01, 0.0, 0.5);
   }

  sub initial {
    P ~ log_normal(log(100.0), 0.2)
    Z ~ log_normal(log(50.0), 0.1)
    ln_alpha ~ gaussian(log(mu), sigma)
  }

  sub transition(delta = h) {
    ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg = 'RK4(3)') {
      dP/dt =  exp(ln_alpha) * P * (1 - P / k1) - b * P * Z
      dZ/dt = -d * Z * (1 + Z / k2) + c * P * Z
      dln_alpha/dt = -sigma * sigma / 2 - sigma * w / h
    }
  }

  sub observation {
    P_obs ~ log_normal(log(P), 0.1)
  }
}
