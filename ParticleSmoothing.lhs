% Naive Particle Smoothing is Degenerate
% Dominic Steinitz
% 21st November 2015

---
bibliography: Stochastic.bib
---

Introduction
============

Let $\{X_t\}_{t \geq 1}$ be a (hidden) Markov process. By hidden, we
mean that we are not able to observe it.

$$
X_1 \sim \mu(\centerdot) \quad X_t \,|\, (X_{t-1} = x) \sim f(\centerdot \,|\, x)
$$

And let $\{Y_t\}_{t \geq 1}$ be an observable Markov process such that

$$
Y_t \,|\, (X_{t} = x) \sim g(\centerdot \,|\, x)
$$

That is the observations are conditionally independent given the state
of the hidden process.

As an example let us take

$$
\begin{aligned}
X_t &= aX_{t-1} + \sigma \epsilon_t \\
Y_t &= bX_t + \tau \eta_t
\end{aligned}
$$

where

$$
A =
\begin{bmatrix}
1 & 0 & \Delta t & 0        \\
0 & 1 & 0        & \Delta t \\
0 & 0 & 1        & 0        \\
0 & 0 & 0        & 1
\end{bmatrix}
, \quad
Q =
\begin{bmatrix}
\frac{q_1^c \Delta t^3}{3} & 0                          & \frac{q_1^c \Delta t^2}{2} & 0                          \\
0                          & \frac{q_2^c \Delta t^3}{3} & 0                          & \frac{q_2^c \Delta t^2}{2} \\
\frac{q_1^c \Delta t^2}{2} & 0                          & {q_1^c \Delta t}           & 0                          \\
0                          & \frac{q_2^c \Delta t^2}{2} & 0                          & {q_2^c \Delta t}
\end{bmatrix}
$$

and

$$
H =
\begin{bmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0
\end{bmatrix}
, \quad
R =
\begin{bmatrix}
\sigma_1^2 & 0 \\
0 & \sigma_2^2
\end{bmatrix}
$$

We wish to determine

$$
p(x_{1:n} \,|\, y_{1:n}) = \frac{p(x_{1:n}, y_{1:n})}{p(y_{1:n})}
$$

By definition we have

$$
{p(y_{1:n})} = \int {p(x_{1:n}, y_{1:n})} \,\mathrm{d}x_{1:n}
$$

And from the Markov and conditional independence we have

$$
{p(x_{1:n}, y_{1:n})} =
\underbrace{\mu(x_1)\prod_{k = 2}^n f(x_k \,|\, x_{k-1})}_{p(x_{1:n})}
\underbrace{\prod_{k = 1}^n g(y_k \,|\, x_k)}_{p(y_{1:n} \,|\, x_{1:n})}
$$

using the Markov property and Chapman-Kolmogorov for the first
factorisation and conditional independence for the second
factorisation.

For finite state models and linear Gaussian models (such the example
above), the posterior can be calculated exactly, see, for example, the
[Kalman
filter](https://idontgetoutmuch.wordpress.com/2014/08/06/fun-with-kalman-filters-part-ii/). In
other cases we need to find a numerical method.

If we could sample $X_{1:n}^{(i)} \sim p(x_{1:n} \,|\, y_{1:n})$ then
we could approximate the posterior as

$$
\hat{p}(x_{1:n} \,|\, y_{1:n}) = \frac{1}{N}\sum_{i=1}^N \delta_{X_{1:n}^{(i)}}(x_{1:n})
$$

If we wish to, we can create marginal estimates

$$
\hat{p}(x_k \,|\, y_{1:n}) = \frac{1}{N}\sum_{i=1}^N \delta_{X_{k}^{(i)}}(x_{k})
$$

When $k = N$, this is the filtering estimate.

Standard Bayesian Recursion
===========================

**Prediction**

$$
\begin{aligned}
p(x_n \,|\, y_{1:n-1}) &= \int p(x_{n-1:n} \,|\, y_{1:n-1}) \,\mathrm{d}x_{n-1} \\
 &= \int p(x_{n} \,|\, x_{n-1}, y_{1:n-1}) \, p(x_{n-1} \,|\, y_{1:n-1}) \,\mathrm{d}x_{n-1} \\
 &= \int f(x_{n} \,|\, x_{n-1}) \, p(x_{n-1} \,|\, y_{1:n-1}) \,\mathrm{d}x_{n-1} \\
\end{aligned}
$$

**Update**

$$
\begin{aligned}
p(x_n \,|\, y_{1:n}) &= \frac{p(y_n \,|\, x_n, y_{1:n-1}) \, p(x_n \,|\, y_{1:n-1})}
                             {p(y_n \,|\, y_{1:n-1})} \\
                     &= \frac{g(y_n \,|\, x_n) \, p(x_n \,|\, y_{1:n-1})}
                             {p(y_n \,|\, y_{1:n-1})}
\end{aligned}
$$

where by definition

$$
{p(y_n \,|\, y_{1:n-1})} = \int {g(y_n \,|\, x_n) \, p(x_n \,|\, y_{1:n-1})} \,\mathrm{d}x_n
$$

Path Space Recursion
--------------------

We have

$$
\begin{aligned}
p(x_{1:n} \,|\, y_{1:n}) &= \frac{x_{1:n}, y_{1:n}}{p(y_{1:n})} \\
&= \frac{p(x_n, y_n \,|\, x_{1:n-1}, y_{1:n-1})}{p(y_{1:n})} \, p(x_{1:n-1}, y_{1:n-1}) \\
&= \frac{p(y_n \,|\, x_{1:n}, y_{1:n-1}) \, p(x_n \,|\, x_{1:n-1}, y_{1:n-1}) }{p(y_{1:n})} \, p(x_{1:n-1}, y_{1:n-1}) \\
&= \frac{g(y_n \,|\, x_{n}) \, f(x_n \,|\, x_{n-1})}{p(y_n \,|\, y_{1:n-1})}
\,
\frac{p(x_{1:n-1}, y_{1:n-1})}{ \, p(y_{1:n-1})} \\
&= \frac{g(y_n \,|\, x_{n}) \, f(x_n \,|\, x_{n-1})}{p(y_n \,|\, y_{1:n-1})}
\,
{p(x_{1:n-1} \,|\,y_{1:n-1})} \\
&= \frac{g(y_n \,|\, x_{n}) \, \overbrace{f(x_n \,|\, x_{n-1}) \, {p(x_{1:n-1} \,|\,y_{1:n-1})}}^{\mathrm{predictive}p(x_{1:n} \,|\, y_{1:n-1})}}
{p(y_n \,|\, y_{1:n-1})} \\
\end{aligned}
$$

where by definition

$$
p(y_n \,|\, y_{1:n-1}) =
\int g(y_n \,|\, x_n) \, p(x_{1:n} \,|\, y_{1:n-1}) \,\mathrm{d}x_{1:n}
$$

**Prediction**

$$
p(x_{1:n} \,|\, y_{1:n-1}) = f(x_n \,|\, x_{n-1}) \, {p(x_{1:n-1} \,|\,y_{1:n-1})}
$$

**Update**

$$
p(x_{1:n} \,|\, y_{1:n}) = \frac{g(y_n \,|\, x_{n}) \, {p(x_{1:n} \,|\, y_{1:n-1})}}
{p(y_n \,|\, y_{1:n-1})}
$$

Algorithm
---------

The idea is to simulate paths using the recursion we derived above.


1. At time $n-1$ we have
$$
\hat{p}(x_{1:n-1} \,|\, y_{n-1}) = \frac{1}{N}\sum_{i=1}^N \delta_{X_{1:n-1}}^{(i)}(x_{1:n-1})
$$

2. Sample $\tilde{X}_n^{(i)} \sim f(\centerdot \,|\, X_{n-1}^{(i)})$
and set $\tilde{X}_{1:n}^{(i)} = (\tilde{X}_{1:n-1}^{(i)},
\tilde{X}_n^{(i)})$. We then have an approximation of the prediction
step
$$
\hat{p}(x_{1:n} \,|\, y_{n-1}) =
\frac{1}{N}\sum_{i=1}^N \delta_{\tilde{X}_{1:n}}^{(i)}(x_{1:n})
$$

Substituting

$$
\begin{aligned}
{\hat{p}(y_n \,|\, y_{1:n-1})} &=
\int {g(y_n \,|\, x_n) \, \hat{p}(x_n \,|\, y_{1:n-1})} \,\mathrm{d}x_n \\
&=
\int {g(y_n \,|\, x_n)}\frac{1}{N}\sum_{i=1}^N \delta_{\tilde{X}_{1:n-1}}^{(i)}(x_{1:n}) \,\mathrm{d}x_n \\
&=
\frac{1}{N}\sum_{i=1}^N {g(y_n \,|\, \tilde{X}_n^{(i)})}
\end{aligned}
$$

FIXME fill in the missing step

$$
\tilde{p}(x_{1:n} \,|\, y_{1:n}) =
\frac{g(y_n \,|\, x_{n}) \, {\hat{p}(x_{1:n} \,|\, y_{1:n-1})}}
     {\hat{p}(y_n \,|\, y_{1:n-1})}
=
\sum_{i=1}^N W_n^{(i)} \delta_{\tilde{X}_{1:n}^{(i)}} (x_{1:n})
$$

3. Now sample
$$
X_{1:n}^{(i)} \sim \tilde{p}(x_{1:n} \,|\, y_{1:n})
$$

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}


> {-# LANGUAGE FlexibleInstances            #-}
> {-# LANGUAGE MultiParamTypeClasses        #-}
> {-# LANGUAGE FlexibleContexts             #-}
> {-# LANGUAGE TypeFamilies                 #-}
> {-# LANGUAGE BangPatterns                 #-}

> module ParticleSmoothing where

> import Data.Random.Source.PureMT
> import Data.Random hiding ( StdNormal, Normal )
> import qualified Data.Random as R
> import Control.Monad.State
> import Data.Number.Erf
> import qualified Numeric.LinearAlgebra.HMatrix as H
> import Foreign.Storable ( Storable )
> import Data.Maybe ( fromJust )
> import Data.Bits ( shiftR )
> import qualified Data.Vector.Unboxed as V
> import Control.Monad.ST
> import System.Random.MWC

> deltaT, sigma1, sigma2, qc1, qc2 :: Double
> deltaT = 0.001
> sigma1 = 1/2
> sigma2 = 1/2
> qc1 = 1
> qc2 = 1

> bigA :: H.Matrix Double
> bigA = (4 H.>< 4) bigAl

> bigAl :: [Double]
> bigAl = [1, 0 , deltaT,      0,
>          0, 1,       0, deltaT,
>          0, 0,       1,      0,
>          0, 0,       0,      1]

> bigQ :: H.Matrix Double
> bigQ = (4 H.>< 4) bigQl

> bigQl :: [Double]
> bigQl = [qc1 * deltaT^3 / 3,                  0, qc1 * deltaT^2 / 2,                  0,
>                           0, qc2 * deltaT^3 / 3,                  0, qc2 * deltaT^2 / 2,
>          qc1 * deltaT^2 / 2,                  0,       qc1 * deltaT,                  0,
>                           0, qc2 * deltaT^2 / 2,                  0,       qc2 * deltaT]

> bigH :: H.Matrix Double
> bigH = (2 H.>< 4) [1, 0, 0, 0,
>                    0, 1, 0, 0]

> bigR :: H.Herm Double
> bigR = H.trustSym $ (2 H.>< 2) [sigma1^2,        0,
>                                        0, sigma2^2]

> normalMultivariate :: H.Vector Double -> H.Herm Double -> RVarT m (H.Vector Double)
> normalMultivariate mu bigSigma = do
>   z <- replicateM (H.size mu) (rvarT R.StdNormal)
>   return $ mu + bigA H.#> (H.fromList z)
>   where
>     (vals, bigU) = H.eigSH bigSigma
>     lSqrt = H.diag $ H.cmap sqrt vals
>     bigA = bigU H.<> lSqrt

> data family Normal k :: *

> data instance Normal (H.Vector Double) = Normal (H.Vector Double) (H.Herm Double)

> instance Distribution Normal (H.Vector Double) where
>   rvar (Normal m s) = normalMultivariate m s

> normalPdf :: (H.Numeric a, H.Field a, H.Indexable (H.Vector a) a, Num (H.Vector a)) =>
>              H.Vector a -> H.Herm a -> H.Vector a -> a
> normalPdf mu sigma x = exp $ normalLogPdf mu sigma x

> normalLogPdf :: (H.Numeric a, H.Field a, H.Indexable (H.Vector a) a, Num (H.Vector a)) =>
>                  H.Vector a -> H.Herm a -> H.Vector a -> a
> normalLogPdf mu bigSigma x = - H.sumElements (H.cmap log (diagonals dec))
>                               - 0.5 * (fromIntegral (H.size mu)) * log (2 * pi)
>                               - 0.5 * s
>   where
>     dec = fromJust $ H.mbChol bigSigma
>     t = fromJust $ H.linearSolve (H.tr dec) (H.asColumn $ x - mu)
>     u = H.cmap (\x -> x * x) t
>     s = H.sumElements u

> diagonals :: (Storable a, H.Element t, H.Indexable (H.Vector t) a) =>
>              H.Matrix t -> H.Vector a
> diagonals m = H.fromList (map (\i -> m H.! i H.! i) [0..n-1])
>   where
>     n = max (H.rows m) (H.cols m)

> instance PDF Normal (H.Vector Double) where
>   pdf (Normal m s) = normalPdf m s
>   logPdf (Normal m s) = normalLogPdf m s

> m0 :: H.Vector Double
> m0 = H.fromList [0, 0, 1, -1]

> bigP0 :: H.Herm Double
> bigP0 = H.trustSym $ H.ident 4

> bar :: IO [H.Vector Double]
> bar = replicateM n $ sample $ rvar (Normal m0 bigP0)

> n :: Int
> n = 500

> indices :: V.Vector Double -> V.Vector Double -> V.Vector Int
> indices bs xs = V.map (binarySearch bs) xs

> binarySearch :: (V.Unbox a, Ord a) =>
>                 V.Vector a -> a -> Int
> binarySearch vec x = loop 0 (V.length vec - 1)
>   where
>     loop !l !u
>       | u <= l    = l
>       | otherwise = let e = vec V.! k in if x <= e then loop l k else loop (k+1) u
>       where k = l + (u - l) `shiftR` 1

> y = undefined

> main :: IO ()
> main = do
>   xTilde1 <- bar
>   let weightUnnormalised = map (normalLogPdf y bigR) $
>                            map (bigH H.#>) xTilde1
>   let vs  :: V.Vector Double
>       vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))
>   return ()

> a = 2
> sigma = 0.1

> foo :: MonadRandom m => Double -> m Double
> foo xPrev = do
>   epsilon <- sample stdNormal
>   let xNew = a * xPrev + sigma * epsilon
>   return xNew

> samples :: (Foldable f, MonadRandom m) =>
>                     (Int -> RVar Double -> RVar (f Double)) ->
>                     Int ->
>                     m (f Double)
> samples repM n = sample $ repM n $ stdNormal

> biggerThan5 :: Int
> biggerThan5 = length (evalState xs (pureMT 42))
>   where
>     xs :: MonadRandom m => m [Double]
>     xs = liftM (filter (>= 5.0)) $ samples replicateM 100000

> biggerThan5' :: Double
> biggerThan5' = sum (evalState xs (pureMT 42)) / (fromIntegral n)
>   where
>     xs :: MonadRandom m => m [Double]
>     xs = liftM (map g) $
>          liftM (filter (>= 5.0)) $
>          liftM (map (+5)) $
>          samples replicateM n
>     g x = exp $ (5^2 / 2) - 5 * x
>     n = 100000

> epsilons :: (Foldable f, MonadRandom m) =>
>                     (Int -> RVar Double -> RVar (f Double)) ->
>                     Double ->
>                     Int ->
>                     m (f Double)
> epsilons repM deltaT n = sample $ repM n $ rvar (R.Normal 0.0 (sqrt deltaT))

> bM0to1 :: Foldable f =>
>           ((Double -> Double -> Double) -> Double -> f Double -> f Double)
>           -> (Int -> RVar Double -> RVar (f Double))
>           -> Int
>           -> Int
>           -> f Double
> bM0to1 scan repM seed n =
>   scan (+) 0.0 $
>   evalState (epsilons repM (recip $ fromIntegral n) n) (pureMT (fromIntegral seed))

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/BrownianPaths.png") 600 600 (translationX 0.0))
```

> p :: Double -> Double -> Double
> p a t = 2 * (1 - normcdf (a / sqrt t))

> m = 10000

> supAbove :: Double -> Double
> supAbove a = fromIntegral count / fromIntegral n
>   where
>     count = length $
>             filter (>= a) $
>             map (\seed -> maximum $ bM0to1 scanl replicateM seed m) [0..n - 1]

> bM0to1WithDrift seed mu n =
>   zipWith (\m x -> x + mu * m * deltaT) [0..] $
>   bM0to1 scanl replicateM seed n
>     where
>       deltaT = recip $ fromIntegral n

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/BrownianWithDriftPaths.png") 600 600 (translationX 0.0))
```

> supAbove' a = (sum $ zipWith (*) ns ws) / fromIntegral n
>   where
>     deltaT = recip $ fromIntegral m
>
>     uss = map (\seed -> bM0to1 scanl replicateM seed m) [0..n - 1]
>     ys = map last uss
>     ws = map (\x -> exp (-a * x - 0.5 * a^2)) ys
>
>     vss = map (zipWith (\m x -> x + a * m * deltaT) [0..]) uss
>     sups = map maximum vss
>     ns = map fromIntegral $ map fromEnum $ map (>=a) sups

Bibliography
============
