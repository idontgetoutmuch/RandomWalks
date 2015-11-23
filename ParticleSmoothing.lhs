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



> module ParticleSmoothing where

> import qualified Data.Vector as V
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
> import Data.Histogram.Fill
> import Data.Histogram.Generic ( Histogram )
> import Data.Number.Erf
> import Data.List ( transpose )

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
> epsilons repM deltaT n = sample $ repM n $ rvar (Normal 0.0 (sqrt deltaT))

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

> n = 500
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
