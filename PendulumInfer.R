library(rstan)
library(mvtnorm)

qc1 = 0.0001
deltaT = 0.01
nSamples = 100
m0 = c(0.016, 0)
g = 9.81

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

etas <- rmvnorm(n=nSamples,mean=c(0.0,0.0),sigma=bigQ)
epsilons <- rmvnorm(n=nSamples,mean=c(0.0),sigma=bigR)

x <- matrix()
length(x) <- 2 * (nSamples + 1)
dim(x) <- c(nSamples + 1, 2)

y <- vector()
length(y) <- 1 * nSamples
dim(y) <- c(nSamples)

x[1,1] = m0[1]
x[1,2] = m0[2]

## x is the state and y is the observable
for (i in 1:nSamples) {
    x[i+1,1] <- x[i,1] + x[i,2] * deltaT
    x[i+1,2] <- x[i,2] - g * (sin(x[i,1])) * deltaT
    x[i+1,] <- x[i+1,] + etas[i,]
    y[i] <- sin(x[i+1,1]) ## + epsilons[i,1]
    }

samples <- stan(file = 'PendulumInfer.stan',
                data = list (T  = nSamples,
                             y0 = x[1,],
                             y  = y,
                             t0 = 0.0,
                             ts = seq(0.01,nSamples/100,0.01)
                             ),
                chains = 4,
                iter = 4000,
                warmup = 2000
                )

s <- extract(samples,permuted=FALSE)
print(dimnames(s))
plot(s[1,1,1:nSamples])

plot(x[,1])
plot(y)

