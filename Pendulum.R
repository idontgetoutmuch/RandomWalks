library(rstan)
library(mvtnorm)

qc1 = 0.01
deltaT = 0.01
nSamples = 500
m0 = c(1.6, 0)
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

y <- matrix()
length(y) <- 1 * nSamples
dim(y) <- c(nSamples,1)

x[1,1] = m0[1]
x[1,2] = m0[2]

for (i in 1:nSamples) {
    x[i+1,1] <- x[i,1] + x[i,2] * deltaT
    x[i+1,2] <- x[i,2] - g * (sin(x[i,1])) * deltaT
    x[i+1,] <- x[i+1,] + etas[i,]
    y[i] <- sin(x[i+1,1]) + epsilons[i,1]
    }

plot(x[,1])
plot(y[,1])

samples <- stan(file = 'Pendulum.stan',
                data = list (T  = nSamples,
                             y0 = c(1.6, 0.0),
                             t0 = 0.0,
                             ts = seq(0.01,5,0.01),
                             theta = array(9.81, dim = 1)
                             ),
                algorithm="Fixed_param",
                chains = 1,
                iter =1
                )

s <- extract(samples,permuted=FALSE)
print(dimnames(s))
plot(s[1,1,1:nSamples])
