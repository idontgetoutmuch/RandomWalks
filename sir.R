library(rstan)

nSamples = 500

samples <- stan(file = 'sir.stan',
                data = list (T  = nSamples,
                             y0 = c(762.0, 1.0, 0.0),
                             ts = seq(0.05,25,0.05),
                             theta = c(0.0026, 0.5)
                             ),
                algorithm="Fixed_param",
                chains = 1,
                iter =1
                )

s <- extract(samples,permuted=FALSE)
print(dimnames(s))
plot(s[1,1,1:nSamples])
plot(s[1,1,(nSamples + 1):(2 * nSamples)], type="l")
plot(s[1,1,(2 * nSamples + 1):(3 * nSamples)])
