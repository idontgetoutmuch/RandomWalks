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
X_1 \sim \mu(\centerdot) \quad X_t \,|\, X_{t-1} \sim f(\centerdot \,|\, x)
$$

> module ParticleSmoothing where

> import qualified Data.Vector as V
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
> import Data.Histogram.Fill
> import Data.Histogram.Generic ( Histogram )
> import Data.Number.Erf
> import Data.List ( transpose )

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
