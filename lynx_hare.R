library(deSolve)
library('optimx')

td = c(0:20)
hare = c(30,47.2,70.2,77.4,36.3,20.6,18.1,21.4,22,25.4,27.1,40.3,57,76.6,52.3,19.5,11.2,7.6,14.6,16.2,24.7)
lynx = c(4,6.1,9.8,35.2,59.4,41.7,19,13,8.3,9.1,7.4,8,12.3,19.5,45.7,51.1,29.7,15.8,9.7,10.1,8.6)
p = c(30,4,0.4,0.018,0.8,0.023)

## function dydt = lotvol(t,y,a1,a2,b1,b2)
##   ## Predator and Prey Model
##   tmp1 = a1*y(1) - a2*y(1)*y(2);
##   tmp2 = -b1*y(2) + b2*y(1)*y(2);
##   dydt = [tmp1; tmp2];
##   end

parametersLV = c(a1=0.4,
                 a2=0.018,
                 b1=0.8,
                 b2=0.023)

stateLV      = c(y1=30,
                 y2=4)

lotvol <- function(t, state, parameters) {
    ## Predator and Prey Model
    with(as.list(c(state, parameters)),{
        dy1 = a1*y1 - a2*y1*y2
        dy2 = -b1*y2 + b2*y1*y2
        list(c(dy1,dy2))
    })
}

timesLV <- seq(0, 20, by = 1.0)

outLV <- ode(y = stateLV, times = timesLV, func = lotvol, parms = parametersLV)

lvForNm <- function(parameters) {
    outLV <- ode(y = stateLV, times = timesLV, func = lotvol, parms = parameters)
    erry1 = outLV[,2] - hare
    erry2 = outLV[,3] - lynx
    J = erry1 %*% erry1 + erry2 %*% erry2
    J[1,1]
}

lvMin <- optimx(parametersLV, lvForNm, method = "Nelder-Mead", control = list(trace=6,maxit=1000))

parametersNM = c(a1=lvMin$a1, a2=lvMin$a2, b1=lvMin$b1, b2=lvMin$b2)

lvMin1 <- optimx(parametersNM, lvForNm, method = "Nelder-Mead", control = list(trace=6,maxit=1000))

parametersNM1 = c(a1=lvMin1$a1, a2=lvMin1$a2, b1=lvMin1$b1, b2=lvMin1$b2)

parametersHMC = c(a1=0.40,
                  a2=0.021,
                  b1=1.0,
                  b2 = 0.027)

parametersHMC2 = c(a1=0.38,
                   a2=0.020,
                   b1=1.1,
                   b2 = 0.029)

parametersHMC3 = c(a1=0.36,
                   a2=0.019,
                   b1=1.2,
                   b2 = 0.031)

parametersHMC3 = c(a1=0.36,
                   a2=0.019,
                   b1=1.2,
                   b2 = 0.030)

initParms <- c(-1.9,-0.1)

g <- function(v){
    x <- v[1]*v[1] + v[2]*v[2]
    return (x)
}

myMin <- optimx(initParms, g, method = "Nelder-Mead", control = list(trace=6,maxit=1000))

myPars <- myMin$p1
myVal  <- myMin$value


parameters <- c(a = -8/3,
                b = -10,
                c= 28)

state <- c(X = 1,
           Y = 1,
           Z = 1)

Lorenz <- function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
        ## rate of change
        dX <- a * X + Y * Z
        dY <- b * (Y - Z)
        dZ <- -X * Y + c * Y - Z
        ## return the rate of change
        list(c(dX, dY, dZ))
    })
}

times <- seq(0, 100, by = 0.01)

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)

par(oma = c(0, 0, 3, 0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "X"], out[, "Z"], pch = ".")
mtext(outer = TRUE, side = 3, "Lorenz model", cex = 1.5)

timesLV <- seq(0, 20, by = 1.0)

outLV <- ode(y = stateLV, times = timesLV, func = lotvol, parms = parametersLV)

par(oma = c(0, 0, 3, 0))
plot(outLV, xlab = "time", ylab = "-")
plot(outLV[, "y1"], outLV[, "y2"], pch = ".")
mtext(outer = TRUE, side = 3, "Lotka-Volterra model", cex = 1.5)


samples <- read.csv('/Users/dom/cmdstan/output.csv', comment.char='#')

png("diagrams/HaresLynxes.png")
plot(main="Hare / Lynx Populations",x=c(1900:1920),xlab="Year",hare,ylim=c(0,100),ylab="Pelts ('000s)",col=3,type="l")
lines(x=c(1900:1920),lynx,col=4)
