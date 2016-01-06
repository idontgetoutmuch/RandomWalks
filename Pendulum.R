library(rstan)

samples <- stan(file = 'Pendulum.stan',
                data = list (T  = 500,
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
plot(s[1,1,1:500])
