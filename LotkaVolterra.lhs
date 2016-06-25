% Modelling an Ecosystem
% Dominic Steinitz
% 23rd April 2016

---
bibliography: Stochastic.bib
---

Introduction
============

In the 1920s, @doi:10.1021/j150111a004 and @Volterra1926 developed a
model of a very simple predator-prey ecosystem.

$$
\begin{eqnarray}
\frac{\mathrm{d}N_1}{\mathrm{d}t} & = & \rho_1 N_1  - c_1 N_1 N_2 \label{eq2a} \\
\frac{\mathrm{d}N_2}{\mathrm{d}t} & = & c_2 N_1 N_2 - \rho_2 N2 \label{eq2b}
\end{eqnarray}
$$

Although simple, it turns out that the Canadian lynx and showshoe hare
are well represented by such a model. Furthermore, the Hudson Bay
Company kept records of how many pelts of each species were trapped
for almost a century, giving a good proxy of the population of each
species.

![](diagrams/HaresLynxes.png)

We can capture the fact that we do not have a complete model by
describing our state of ignorance about the parameters. In order to
keep this as simple as possible let us assume that log parameters
undergo Brownian motion. That is, we know the parameters will jiggle
around and the further into the future we look the less certain we are
about what values they will have taken. By making the log parameters
undergo Brownian motion, we can also capture our modelling assumption
that birth, death and predation rates are always positive. A similar
approach is taken in @Dureau2013 where the (log) parameters of an
epidemiological model are taken to be Ornstein-Uhlenbeck processes
(which is biologically more plausible although adds to the complexity
of the model, something we wish to avoid in an example such as this).

@Andrieu2010 propose a method to estimate the parameters of such
models and the domain specific probabilistic language
[LibBi](http://libbi.org/) (@Murray) can be used to apply this (and
other inference methods).

A Dynamical System Aside
========================

The above dynamical system is structurally unstable (more on this in a
future post), a possible indication that it should not be considered
as a good model of predator–prey interaction. Let us modify this to
include carrying capacities for the populations of both species.

$$
\begin{eqnarray}
\frac{\mathrm{d}N_1}{\mathrm{d}t} & = & \rho_1 N_1 \bigg(1 - \frac{N_1}{K_1}\bigg) - c_1 N_1 N_2 \\
\frac{\mathrm{d}N_2}{\mathrm{d}t} & = & -\rho_2 N_2 \bigg(1 + \frac{N_2}{K_2}\bigg) + c_2 N_1 N_2
\end{eqnarray}
$$

Data Generation with LibBi
==========================

Let's generate some data using LibBi.

~~~~{.CPP include="PP.bi"}
~~~~

![](diagrams/LVdata.png)

We can look at phase space starting with different populations and see
they all converge to the same fixed point.

![](diagrams/PPviaLibBi.png)


Data Generation with Haskell
============================

Since at some point in the future, I plan to produce Haskell versions
of the methods given in @Andrieu2010, let's generate the data using
Haskell.

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}

> module LotkaVolterra (
>     solLv
>   , solPp
>   , h0
>   , l0
>   , baz
>   , logBM
>   )where

> import Numeric.GSL.ODE
> import Numeric.LinearAlgebra

> import Data.Random.Source.PureMT
> import Data.Random hiding ( gamma )
> import Control.Monad.State

Here's the unstable model.

> lvOde :: Double ->
>          Double ->
>          Double ->
>          Double ->
>          Double ->
>          [Double] ->
>          [Double]
> lvOde rho1 c1 rho2 c2 _t [h, l] =
>   [
>     rho1 * h - c1 * h * l
>   , c2 * h * l - rho2 * l
>   ]
> lvOde _rho1 _c1 _rho2 _c2 _t vars =
>   error $ "lvOde called with: " ++ show (length vars) ++ " variable"

> rho1, c1, rho2, c2 :: Double
> rho1 = 0.5
> c1 = 0.02
> rho2 = 0.4
> c2 = 0.004

> deltaT :: Double
> deltaT = 0.1

> solLv :: Matrix Double
> solLv = odeSolve (lvOde rho1 c1 rho2 c2)
>                  [50.0, 50.0]
>                  (fromList [0.0, deltaT .. 50])

![](diagrams/LV.png)

And here's the stable model.

> ppOde :: Double ->
>          Double ->
>          Double ->
>          Double ->
>          Double ->
>          Double ->
>          Double ->
>          [Double] ->
>          [Double]
> ppOde a k1 b d k2 c _t [p, z] =
>   [
>     a * p * (1 - p / k1) - b * p * z
>   , -d * z * (1 + z / k2) + c * p * z
>   ]
> ppOde _a _k1 _b _d _k2 _c _t vars =
>   error $ "ppOde called with: " ++ show (length vars) ++ " variable"

> a, k1, b, d, k2, c :: Double
> a = 0.5
> k1 = 200.0
> b = 0.02
> d = 0.4
> k2 = 50.0
> c = 0.004

> solPp :: Double -> Double -> Matrix Double
> solPp x y = odeSolve (ppOde a k1 b d k2 c)
>                  [x, y]
>                  (fromList [0.0, deltaT .. 50])

> gamma, alpha, beta :: Double
> gamma = d / a
> alpha = a / (c * k1)
> beta  = d / (a * k2)

> fp :: (Double, Double)
> fp = ((gamma + beta) / (1 + alpha * beta), (1 - gamma * alpha) / (1 + alpha * beta))

> h0, l0 :: Double
> h0 = a * fst fp / c
> l0 = a * snd fp / b

> foo, bar :: Matrix R
> foo = matrix 2 [a / k1, b, c, negate d / k2]
> bar = matrix 1 [a, d]

> baz :: Maybe (Matrix R)
> baz = linearSolve foo bar

This gives a stable fixed point of

    [ghci]
    baz

Here's an example of convergence to that fixed point in phase space.

![](diagrams/PP.png)

The Stochastic Model
--------------------

Let us now assume that the Hare growth parameter undergoes Brownian
motion so that the further into the future we go, the less certain we
are about it. In order to ensure that this parameter remains positive,
let's model the log of it to be Brownian motion.

$$
\begin{eqnarray}
\frac{\mathrm{d}N_1}{\mathrm{d}t} & = & \rho_1 N_1 \bigg(1 - \frac{N_1}{K_1}\bigg) - c_1 N_1 N_2 \\
\frac{\mathrm{d}N_2}{\mathrm{d}t} & = & -\rho_2 N_2 \bigg(1 + \frac{N_2}{K_2}\bigg) + c_2 N_1 N_2 \\
\mathrm{d} \rho_1 & = & \rho_1 \sigma_{\rho_1} \mathrm{d}W_t
\end{eqnarray}
$$

where the final equation is a stochastic differential equation with
$W_t$ being a Wiener process.

By Itô we have

$$
\mathrm{d} (\log{\rho_1}) = - \frac{\sigma_{\rho_1}^2}{2} \mathrm{d} t + \sigma_{\rho_1} \mathrm{d}W_t
$$

We can use this to generate paths for $\rho_1$.

$$
\rho_1(t + \Delta t) = \rho_1(t)\exp{\bigg(- \frac{\sigma_{\rho_1}^2}{2} \Delta t + \sigma_{\rho_1} \sqrt{\Delta t} Z\bigg)}
$$

where $Z \sim {\mathcal{N}}(0,1)$.

> oneStepLogBM :: MonadRandom m => Double -> Double -> Double -> m Double
> oneStepLogBM deltaT sigma rhoPrev = do
>   x <- sample $ rvar StdNormal
>   return $ rhoPrev * exp(sigma * (sqrt deltaT) * x - 0.5 * sigma * sigma * deltaT)

> iterateM :: Monad m => (a -> m a) -> m a -> Int -> m [a]
> iterateM f mx n = sequence . take n . iterate (>>= f) $ mx

> logBMM :: MonadRandom m => Double -> Double -> Int -> Int -> m [Double]
> logBMM initRho sigma n m =
>   iterateM (oneStepLogBM (recip $ fromIntegral n) sigma) (return initRho) (n * m)

> logBM :: Double -> Double -> Int -> Int -> Int -> [Double]
> logBM initRho sigma n m seed =
>   evalState (logBMM initRho sigma n m) (pureMT $ fromIntegral seed)

![](diagrams/LogBrownianPaths.png)


Bibliography
============