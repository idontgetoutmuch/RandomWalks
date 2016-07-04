functions {
  real f1 (real a, real k1, real b, real p, real z) {
    real q;

    q = a * p * (1 - p / k1) - b * p * z;
    return q;
  }

  real f2 (real d, real k2, real c, real p, real z) {
    real q;

    q = -d * z * (1 + z / k2) + c * p * z;
    return q;
  }
}

data {
  int<lower=1> T;   // Number of observations
  real y[T];        // Observed hares
  real k1;          // Hare carrying capacity
  real b;           // Hare death rate per lynx
  real d;           // Lynx death rate
  real k2;          // Lynx carrying capacity
  real c;           // Lynx birth rate per hare
  real deltaT;      // Time step
}

parameters {
  real<lower=0> mu;
  real<lower=0> sigma;
  real<lower=0> p0;
  real<lower=0> z0;
  real<lower=0> rho0;
  real w[T];
}

transformed parameters {
  real<lower=0> p[T];
  real<lower=0> z[T];
  real<lower=0> rho[T];

  p[1] = p0;
  z[1] = z0;
  rho[1] = rho0;

  for (t in 1:(T-1)){
    p[t+1] = p[t] + deltaT * f1 (exp(rho[t]), k1, b, p[t], z[t]);
    z[t+1] = z[t] + deltaT * f2 (d, k2, c, p[t], z[t]);

    rho[t+1] = rho[t] * exp(sigma * sqrt(deltaT) * w[t] - 0.5 * sigma * sigma * deltaT);
  }
}

model {
  mu    ~ uniform(0.0,1.0);
  sigma ~ uniform(0.0, 0.5);
  p0    ~ lognormal(log(100.0), 0.2);
  z0    ~ lognormal(log(50.0), 0.1);
  rho0  ~ normal(log(mu), sigma);
  w     ~ normal(0.0,1.0);

  for (t in 1:T) {
    y[t] ~ lognormal(log(p[t]),0.1);
  }
}

