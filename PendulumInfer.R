library(rstan)

qc1 = 0.0001
deltaT = 0.01
nSamples = 100
m0 = c(1.07, 0)
g = 9.81
t0 = 0.0
ts = seq(deltaT,nSamples * deltaT,deltaT)

bigQ = matrix(c(qc1 * deltaT^3 / 3, qc1 * deltaT^2 / 2,
                qc1 * deltaT^2 / 2,       qc1 * deltaT
                ),
              nrow = 2,
              ncol = 2,
              byrow = TRUE
              )

bigR = matrix(c(0.1),
              nrow = 1,
              ncol = 1,
              byrow = TRUE)

## Generate the angle and angular velocity of pendulum with added noise
samples <- stan(file = 'Pendulum.stan',
                data = list (T  = nSamples,
                             y0 = m0,
                             t0 = t0,
                             ts = ts,
                             theta = array(g, dim = 1),
                             sigma = c(bigQ[1,1], bigQ[2,2])
                             ),
                algorithm="Fixed_param",
                seed = 42,
                chains = 1,
                iter =1
                )

## We can only observe the horizontal displacement
s <- extract(samples,permuted=FALSE)
zStan <- sin(s[1,1,1:nSamples])

estimates <- stan(file = 'PendulumInfer.stan',
                  data = list (T  = nSamples,
                               y0 = m0,
                               z  = zStan,
                               t0 = t0,
                               ts = ts
                               ),
                  seed = 42,
                  chains = 1,
                  iter = 1000,
                  warmup = 500
                )

set.seed(42)
etas <- rmvnorm(n=nSamples,mean=c(0.0,0.0),sigma=bigQ)

y <- matrix()
length(y) <- 2 * (nSamples + 1)
dim(y) <- c(nSamples + 1, 2)

z <- vector()
length(z) <- 1 * nSamples
dim(z) <- c(nSamples)

y[1,1] = m0[1]
y[1,2] = m0[2]

## y is the state and z is the observable
for (i in 1:nSamples) {
    y[i+1,1] <- y[i,1] + y[i,2] * deltaT
    y[i+1,2] <- y[i,2] - g * (sin(y[i,1])) * deltaT
    y[i+1,] <- y[i+1,] + etas[i,]
    z[i] <- sin(y[i+1,1])
    }

estimatesfromR <- stan(file = 'PendulumInfer.stan',
                       data = list (T  = nSamples,
                                    y0 = m0,
                                    z  = z,
                                    t0 = t0,
                                    ts = ts
                                    ),
                       seed = 42,
                       chains = 1,
                       iter = 1000,
                       warmup = 500
                       )
