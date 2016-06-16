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
\frac{\mathrm{d}x}{\mathrm{d}t} & = & \beta_{00} x  - \beta{01} xy \label{eq2a} \\
\frac{\mathrm{d}y}{\mathrm{d}t} & = & \beta_{11} xy - \beta{10} y \label{eq2b}
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
of the model something we wish to avoid in an example such as this).

@Andrieu2010 propose a method to estimate the parameters of such
models and the domain specific probabilistic language LibBi (@Murray)
can be used to apply this (and other inference methods).

Some Typical Data
=================

Rather than start with actual data for which we do not know the
parameters, let us start with some generated data and add some
noise. Since at some point in the future, I plan to produce Haskell
versions of the methods given in @Andrieu2010, let's generate the data
using Haskell.

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}

> module LotkaVolterra ( solLv )where

> import Numeric.GSL.ODE
> import Numeric.LinearAlgebra hiding ( R, vector, matrix, sym )

> lvOde :: Double ->
>          Double ->
>          Double ->
>          Double ->
>          Double ->
>          [Double] ->
>          [Double]
> lvOde beta00 beta01 beta10 beta11 _t [h, l] =
>   [
>     beta00 * h - beta01 * h * l
>   , beta11 * h * l - beta10 * l
>   ]
> lvOde _beta00 _beta01 _beta10 _beta11 _t vars =
>   error $ "lvOde called with: " ++ show (length vars) ++ " variable"

> beta00, beta01, beta10, beta11 :: Double

> beta00 = 0.5
> beta01 = 0.02
> beta10 = 0.4
> beta11 = 0.004

> deltaT :: Double
> deltaT = 0.1

> solLv :: Matrix Double
> solLv = odeSolve (lvOde beta00 beta01 beta10 beta11)
>                  [50.0, 50.0]
>                  (fromList [0.0, deltaT .. 50])

The Model Expanded
------------------

$$
\begin{eqnarray}
\frac{\mathrm{d}x}{\mathrm{d}t} & = & \beta_{11} x  - \beta_{12} xy \\
\frac{\mathrm{d}y}{\mathrm{d}t} & = & \beta_{22} xy - \beta_{21} y \\
\frac{\mathrm{d}\beta_{11}}{\mathrm{d}t} & = & \theta_{{11}}\mathrm{d}W_{{11}}(t) \\
\frac{\mathrm{d}\beta_{12}}{\mathrm{d}t} & = & \theta_{{12}}\mathrm{d}W_{{12}}(t) \\
\frac{\mathrm{d}\beta_{21}}{\mathrm{d}t} & = & \theta_{{21}}\mathrm{d}W_{{21}}(t) \\
\frac{\mathrm{d}\beta_{22}}{\mathrm{d}t} & = & \theta_{{22}}\mathrm{d}W_{{22}}(t)
\end{eqnarray}
$$

LibBi
-----

~~~~{.CPP include="LV.bi"}
~~~~


Bibliography
============