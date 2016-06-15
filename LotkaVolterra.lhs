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
\frac{\mathrm{d}x}{\mathrm{d}t} & = & a_1 x  - a_2 xy \label{eq2a} \\
\frac{\mathrm{d}y}{\mathrm{d}t} & = & b_2 xy - b_1 y \label{eq2b}
\end{eqnarray}
$$

Although simple, it turns out that the Canadian lynx and showshoe hare
are well represented by such a model. Furthermore, the Hudson Bay
Company kept records of how many pelts of each species were trapped
for almost a century, giving a good proxy of the population of each
species.

![](diagrams/HaresLynxes.png)

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}


> module LotkaVolterra where

> import Numeric.GSL.ODE
> import Numeric.LinearAlgebra hiding ( R, vector, matrix, sym )


> lvOde :: Double -> Double -> Double -> Double -> Double -> [Double] -> [Double]
> lvOde a1 a2 b1 b2 _t [h, l] =
>   [
>     a1 * h - a2 * h * l
>   , b2 * h * l - b1 * l
>   ]
> lvOde _a1 _a2 _b1 _b2 _t vars = error $ "lvOde called with: " ++ show (length vars) ++ " variable"

> sirOde :: Double -> Double -> Double -> [Double] -> [Double]
> sirOde delta gamma _t [s, i, _r] =
>   [
>     negate (delta * i * s)
>   , (delta * i * s) - (gamma * i)
>   , gamma * i
>   ]
> sirOde _b _g _t vars = error $ "sirOde called with: " ++ show (length vars) ++ " variable"

> delta, gamma :: Double
> delta = 0.0026
> gamma = 0.5

> initS, initI, initR :: Double
> initS = 762.0
> initI = 1.0
> initR = 0.01

> a1, a2, b1, b2 :: Double

a1 = 0.5
a2 = 0.02
b1 = 0.4
b2 = 0.004

> a1 = 0.7509811
> a2 = 0.2133682
> b1 = 0.6937935
> b2 = 0.6497548

> deltaT :: Double
> deltaT = 1.0

> sol :: Matrix Double
> sol = odeSolve (sirOde delta gamma) [initS, initI, initR] (fromList [0.0,deltaT..14.0])

> solLv :: Matrix Double
> solLv = odeSolve (lvOde a1 a2 b1 b2) [50.0, 50.0] (fromList [0.0,0.1..50])

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