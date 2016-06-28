/**
 * Lotka-Volterra-like phytoplankton-zooplankton (PZ) model.
 *
 * @author Lawrence Murray <lawrence.murray@csiro.au>
 * $Rev$
 * $Date$
 */
model PZ {
  const h         = 0.1;   // time step
  const delta_abs = 1.0e-3; // absolute error tolerance
  const delta_rel = 1.0e-6; // relative error tolerance
  const epsilon   = 2.0e-2; // observation error tolerance

  const c = 0.02 // 0.25 // 0.02   // zooplankton clearance rate
  const e = 0.30 // 3.0    // zooplankton growth efficiency
  const m_l = 0.1  // zooplankton linear mortality
  const m_q = 0.0 // 0.1  // zooplankton quadratic mortality

  const bb = 0.02
  const cc = 0.4
  const dd = 0.004

  param mu, sigma  // mean and standard deviation of phytoplankton growth
  state P, Z       // phytoplankton, zooplankton
  noise w          // noise term
  state ln_alpha   // stochastic phytoplankton growth rate
  obs P_obs        // observations of phytoplankton

  sub parameter {
    mu ~ uniform(0.0, 1.0)
    sigma ~ uniform(0.0, 0.5)
  }

  sub proposal_parameter {
    mu ~ truncated_gaussian(mu, 0.02, 0.0, 1.0);
    sigma ~ truncated_gaussian(sigma, 0.01, 0.0, 0.5);
  }

  sub initial {
    P ~ log_normal(log(20.0), 0.2)
    Z ~ log_normal(log(20.0), 0.1)
    ln_alpha ~ gaussian(log(0.5), 0.001)
    // ln_alpha ~ gaussian(1.454322, 0.001)
    // ln_alpha ~ gaussian(mu,sigma)
  }

  sub transition(delta = h) {
    w ~ normal(0.0, sqrt(h));
    ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg = 'RK4(3)') {
      // dP/dt = exp(ln_alpha) * P - c * P * Z
      dP/dt = exp(ln_alpha) * P - bb * P * Z
      // dZ/dt = e * c * P * Z - m_l * Z - m_q * Z * Z
      dZ/dt = dd * P * Z - cc * Z
      dln_alpha/dt = 0.0 // -sigma * sigma / 2 - sigma * w / h
    }
  }

  sub observation {
    P_obs ~ log_normal(log(P), 0.01)
  }
}
