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
> {-# LANGUAGE GeneralizedNewtypeDeriving   #-}
> {-# LANGUAGE ScopedTypeVariables          #-}

> module ParticleSmoothing where

> import Data.Random.Source.PureMT
> import Data.Random hiding ( StdNormal, Normal )
> import qualified Data.Random as R
> import Control.Monad.State
> import Control.Monad.Writer hiding ( Any, All )
> import qualified Numeric.LinearAlgebra.HMatrix as H
> import Foreign.Storable ( Storable )
> import Data.Maybe ( fromJust )
> import Data.Bits ( shiftR )
> import qualified Data.Vector as V
> import Control.Monad.ST
> import System.Random.MWC

> import           Data.Array.Repa ( Z(..), (:.)(..), Any(..), computeP,
>                                    extent, DIM1, DIM2, slice, All(..)
>                                  )
> import qualified Data.Array.Repa as Repa

> import qualified Control.Monad.Loops as ML

> import Debug.Trace
> import PrettyPrint ()
> import Text.PrettyPrint
> import Text.PrettyPrint.HughesPJClass ( pPrint )

> deltaT, sigma1, sigma2, qc1, qc2 :: Double
> deltaT = 0.1
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

> bigQ :: H.Herm Double
> bigQ = H.trustSym $ (4 H.>< 4) bigQl

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

> m0 :: H.Vector Double
> m0 = H.fromList [0, 0, 1, -1]

> bigP0 :: H.Herm Double
> bigP0 = H.trustSym $ H.ident 4

> indices :: V.Vector Double -> V.Vector Double -> V.Vector Int
> indices bs xs = V.map (binarySearch bs) xs

> binarySearch :: Ord a =>
>                 V.Vector a -> a -> Int
> binarySearch vec x = loop 0 (V.length vec - 1)
>   where
>     loop !l !u
>       | u <= l    = l
>       | otherwise = let e = vec V.! k in if x <= e then loop l k else loop (k+1) u
>       where k = l + (u - l) `shiftR` 1

> type SystemState = (Double, Double, Double, Double)

> type ArraySmoothing = Repa.Array Repa.U DIM2

> singleStep :: ArraySmoothing SystemState -> H.Vector Double ->
>               WriterT [ArraySmoothing SystemState] (StateT PureMT IO) (ArraySmoothing SystemState)
> singleStep x y = do
>   let (Z :. ix :. jx) = extent x
>
>   xTildeNextH <- lift $ do
>     xHatR :: Repa.Array Repa.U DIM1 SystemState <- computeP $ Repa.slice x (Any :. jx - 1)
>     let xHatH = map (\(a, b, c, d) -> H.fromList [a, b, c, d]) $ Repa.toList xHatR
>     xTildeNextH <- mapM (\x -> sample $ rvar (Normal (bigA H.#> x) bigQ)) xHatH
>     return xTildeNextH
>
>   let systemState = map ((\[a,b,c,d] -> (a, b, c, d)) . H.toList) xTildeNextH
>   -- trace ("\nSystem state: " ++ show systemState) $ return ()
>   tell [x]
>
>   lift $ do
>   let xTildeNextR = Repa.fromListUnboxed (Z :. ix :. (1 :: Int)) $
>                     systemState
>       xTilde = Repa.append x xTildeNextR
>
>       weights = map (normalPdf y bigR) $
>                 map (bigH H.#>) xTildeNextH
>       vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))
>       cumSumWeights = V.scanl (+) 0 (V.fromList weights)
>       js = indices (V.tail cumSumWeights) vs
>   -- foo :: Repa.Array Repa.U DIM2 SystemState <- computeP xTilde
>   -- trace (show foo ++ " " ++ show js ++ " " ++ show ix ++ " " ++ show jx) $ return ()
>   let xNewV = V.map (\j -> Repa.transpose $
>                            Repa.reshape (Z :. (1 :: Int) :. jx + 1) $
>                            slice xTilde (Any :. j :. All)) js
>       xNewR = Repa.transpose $ V.foldr Repa.append (xNewV V.! 0) (V.tail xNewV)
>   computeP xNewR

> n :: Int
> n = 23

> bigA1 = 0.5
> bigQ1 = 0.1
> bigR1 = 0.1
> bigH1 = 1.0

> singleStep1 :: ArraySmoothing Double -> Double ->
>               WriterT [ArraySmoothing Double] (StateT PureMT IO) (ArraySmoothing Double)
> singleStep1 x y = do
>   let (Z :. ix :. jx) = extent x
>
>   xTildeNextH <- lift $ do
>     xHatR :: Repa.Array Repa.U DIM1 Double <- computeP $ Repa.slice x (Any :. jx - 1)
>     let xHatH = Repa.toList xHatR
>     xTildeNextH <- mapM (\x -> sample $ rvar (R.Normal (bigA1 * x) bigQ1)) xHatH
>     return xTildeNextH
>
>   let systemState = xTildeNextH
>   -- trace ("\nSystem state: " ++ show systemState) $ return ()
>   tell [x]
>
>   lift $ do
>   let xTildeNextR = Repa.fromListUnboxed (Z :. ix :. (1 :: Int)) $
>                     systemState
>       xTilde = Repa.append x xTildeNextR
>
>       weights = map (pdf (R.Normal y bigR1)) $
>                 map (bigH1 *) xTildeNextH
>       vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))
>       cumSumWeights = V.scanl (+) 0 (V.fromList weights)
>       totWeight = sum weights
>       js = indices (V.map (/ totWeight) $ V.tail cumSumWeights) vs
>   -- foo :: Repa.Array Repa.U DIM2 SystemState <- computeP xTilde
>   -- trace (show foo ++ " " ++ show js ++ " " ++ show ix ++ " " ++ show jx) $ return ()
>   let xNewV = V.map (\j -> Repa.transpose $
>                            Repa.reshape (Z :. (1 :: Int) :. jx + 1) $
>                            slice xTilde (Any :. j :. All)) js
>       xNewR = Repa.transpose $ V.foldr Repa.append (xNewV V.! 0) (V.tail xNewV)
>   -- trace ("\njs " ++ show js) $ return ()
>   computeP xNewR

> carSample :: MonadRandom m =>
>              H.Vector Double ->
>              m (Maybe ((H.Vector Double, H.Vector Double), H.Vector Double))
> carSample xPrev = do
>   xNew <- sample $ rvar (Normal (bigA H.#> xPrev) bigQ)
>   yNew <- sample $ rvar (Normal (bigH H.#> xNew) bigR)
>   return $ Just ((xNew, yNew), xNew)

> carSamples :: [(H.Vector Double, H.Vector Double)]
> carSamples = evalState (ML.unfoldrM carSample m0) (pureMT 17)

> carSample1 :: MonadRandom m =>
>               Double ->
>               m (Maybe ((Double, Double), Double))
> carSample1 xPrev = do
>   xNew <- sample $ rvar (R.Normal (bigA1 * xPrev) 0.1 {- bigQ1 -})
>   yNew <- sample $ rvar (R.Normal (bigH1 * xNew) bigR1)
>   return $ Just ((xNew, yNew), xNew)

> carSamples1 :: [(Double, Double)]
> carSamples1 = evalState (ML.unfoldrM carSample1 0.0) (pureMT 17)

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/CarPosition.png") 600 600 (translationX 0.0))
```

> y :: H.Vector Double
> y = snd $ head carSamples

> y1 :: Double
> y1 = snd $ head carSamples1

> initXHat :: MonadRandom m => m (ArraySmoothing SystemState)
> initXHat = do
>   xTilde1 <- replicateM n $ sample $ rvar (Normal m0 bigP0)
>   -- trace ("\nxTilde1 " ++ show xTilde1) $ return ()
>   let weights = map (normalPdf y bigR) $
>                 map (bigH H.#>) xTilde1
>       vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))
>       cumSumWeights = V.scanl (+) 0 (V.fromList weights)
>       js = indices (V.tail cumSumWeights) vs
>       xHat1 = Repa.fromListUnboxed (Z :. n :. (1 :: Int)) $
>               map ((\[a,b,c,d] -> (a, b, c, d)) . H.toList) $
>               V.toList $
>               V.map ((V.fromList xTilde1) V.!) js
>   -- trace ("\nWeights1 " ++ show weights) $ return ()
>   return xHat1

> initXHat1 :: MonadRandom m => m (ArraySmoothing Double)
> initXHat1 = do
>   xTilde1 <- replicateM n $ sample $ rvar $ R.Normal y1 bigR1
>   -- trace ("\nxTilde1 " ++ show xTilde1) $ return ()
>   let weights = map (pdf (R.Normal y1 bigR1)) $
>                 map (bigH1 *) xTilde1
>       totWeight = sum weights
>       vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))
>       cumSumWeights = V.scanl (+) 0 (V.fromList $ map (/ totWeight) weights)
>       js = indices (V.tail cumSumWeights) vs
>       xHat1 = Repa.fromListUnboxed (Z :. n :. (1 :: Int)) $
>               V.toList $
>               V.map ((V.fromList xTilde1) V.!) js
>   -- trace ("cumSumWeights " ++ show cumSumWeights) $ return ()
>   -- trace ("\nWeights1 " ++ show weights) $ return ()
>   -- trace ("\njs1 " ++ show js) $ return ()
>   return xHat1

> test :: IO ()
> test = do
>   states <- snd <$> evalStateT smoother (pureMT 24)
>   putStrLn "States"
>   mapM_ (putStrLn . render . pPrint) states

> test1 :: IO [[Double]]
> test1 = do
>   states <- snd <$> evalStateT smoother1 (pureMT 24)
>   let foo :: Int -> IO (Repa.Array Repa.U DIM1 Double)
>       foo i = computeP $ Repa.slice (last states) (Any :. i :. All)
>   bar <- mapM foo [0..22]
>   let baz :: [[Double]]
>       baz = map Repa.toList bar
>   return baz
>   -- error (show bar)
>   -- error (show $ extent $ (last states))
>   -- mapM_ (putStrLn . render . pPrint) states

> smoother :: StateT PureMT IO (ArraySmoothing SystemState, [ArraySmoothing SystemState])
> smoother = runWriterT $ do
>   xHat1 <- lift initXHat
>   foldM singleStep xHat1 (take 3 $ map snd $ tail carSamples)

> smoother1 :: StateT PureMT IO (ArraySmoothing Double, [ArraySmoothing Double])
> smoother1 = runWriterT $ do
>   xHat1 <- lift initXHat1
>   foldM singleStep1 xHat1 (take 20 $ map snd $ tail carSamples1)

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/Smooth.png") 600 600 (translationX 0.0))
```

Notes
=====

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

Bibliography
============
