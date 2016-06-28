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

As an example let us take the one given in @Srkk:2013:BFS:2534502
where the movement of a car is given by Newton's laws of motion and
the acceleration is modelled as white noise.

$$
\begin{aligned}
X_t &= AX_{t-1} + Q \epsilon_t \\
Y_t &= HX_t + R \eta_t
\end{aligned}
$$

Although we do not do so here, $A, Q, H$ and $R$ can be derived from
the dynamics. For the purpose of this blog post, we note that they are
given by

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

We wish to determine the position and velocity of the car given noisy
observations of the position. In general we need the distribution of
the hidden path given the observable path. We use the notation
$x_{m:n}$ to mean the path of $x$ starting a $m$ and finishing at $n$.

$$
p(x_{1:n} \,|\, y_{1:n}) = \frac{p(x_{1:n}, y_{1:n})}{p(y_{1:n})}
$$

Haskell Preamble
----------------

{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
{-# OPTIONS_GHC -fno-warn-missing-methods #-}
{-# OPTIONS_GHC -fno-warn-orphans         #-}


> {-# LANGUAGE FlexibleInstances            #-}
> {-# LANGUAGE MultiParamTypeClasses        #-}
> {-# LANGUAGE FlexibleContexts             #-}
> {-# LANGUAGE TypeFamilies                 #-}
> {-# LANGUAGE BangPatterns                 #-}
> {-# LANGUAGE GeneralizedNewtypeDeriving   #-}
> {-# LANGUAGE ScopedTypeVariables          #-}
> {-# LANGUAGE TemplateHaskell              #-}
> {-# LANGUAGE DataKinds                    #-}

> module ParticleSmoothingII
>   ( simpleSamples
>   , carSamples
>   , testCar
>   , nSimples
>   , testSimple
>   , testSimple1
>   , testSimple2
>   , smoothBackwards
>   , testSmoothSimple
>   , testSmoothCar
>   ) where

> import Data.Random.Source.PureMT
> import Data.Random hiding ( StdNormal, Normal )
> import qualified Data.Random as R
> import Control.Monad.State
> import Control.Monad.Writer hiding ( Any, All )
> import qualified Numeric.LinearAlgebra.HMatrix as H
> import Numeric.LinearAlgebra.HMatrix ( (#>) )
> import qualified Numeric.LinearAlgebra.Static as S
> import Foreign.Storable ( Storable )
> import Data.Maybe ( fromJust )
> import Data.Bits ( shiftR )
> import qualified Data.Vector as V
> import qualified Data.Vector.Unboxed as U
> import Control.Monad.ST
> import System.Random.MWC

> import Data.Array.Repa ( Z(..), (:.)(..), Any(..), computeP,
>                          extent, DIM1, DIM2, slice, All(..)
>                        )
> import qualified Data.Array.Repa as Repa

> import qualified Control.Monad.Loops as ML

> import PrettyPrint ()
> import Text.PrettyPrint ( render )
> import Text.PrettyPrint.HughesPJClass ( Pretty, pPrint )

> import Data.Vector.Unboxed.Deriving
> import Foreign.Storable.Record as Store
> import Foreign.Storable ( Storable (..) )

> import Data.List ( scanl' )

> import Debug.Trace


Some Theory
===========

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
---------------------------

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
p(x_{1:n} \,|\, y_{1:n}) &= \frac{p(x_{1:n}, y_{1:n})}{p(y_{1:n})} \\
&= \frac{p(x_n, y_n \,|\, x_{1:n-1}, y_{1:n-1})}{p(y_{1:n})} \, p(x_{1:n-1}, y_{1:n-1}) \\
&= \frac{p(y_n \,|\, x_{1:n}, y_{1:n-1}) \, p(x_n \,|\, x_{1:n-1}, y_{1:n-1}) }{p(y_{1:n})} \, p(x_{1:n-1}, y_{1:n-1}) \\
&= \frac{g(y_n \,|\, x_{n}) \, f(x_n \,|\, x_{n-1})}{p(y_n \,|\, y_{1:n-1})}
\,
\frac{p(x_{1:n-1}, y_{1:n-1})}{ \, p(y_{1:n-1})} \\
&= \frac{g(y_n \,|\, x_{n}) \, f(x_n \,|\, x_{n-1})}{p(y_n \,|\, y_{1:n-1})}
\,
{p(x_{1:n-1} \,|\,y_{1:n-1})} \\
&= \frac{g(y_n \,|\, x_{n}) \, \overbrace{f(x_n \,|\, x_{n-1}) \, {p(x_{1:n-1} \,|\,y_{1:n-1})}}^{\mathrm{predictive}\,p(x_{1:n} \,|\, y_{1:n-1})}}
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


At time $n-1$ we have an approximating distribution

$$
\hat{p}(x_{1:n-1} \,|\, y_{1:n-1}) = \frac{1}{N}\sum_{i=1}^N \delta_{X_{1:n-1}}^{(i)}(x_{1:n-1})
$$

Sample $\tilde{X}_n^{(i)} \sim f(\centerdot \,|\, X_{n-1}^{(i)})$
and set $\tilde{X}_{1:n}^{(i)} = (\tilde{X}_{1:n-1}^{(i)},
\tilde{X}_n^{(i)})$. We then have an approximation of the prediction
step.

$$
\hat{p}(x_{1:n} \,|\, y_{1:n-1}) =
\frac{1}{N}\sum_{i=1}^N \delta_{\tilde{X}_{1:n}}^{(i)}(x_{1:n})
$$

Substituting

$$
\begin{aligned}
{\hat{p}(y_n \,|\, y_{1:n-1})} &=
\int {g(y_n \,|\, x_n) \, \hat{p}(x_{1:n} \,|\, y_{1:n-1})} \,\mathrm{d}x_n \\
&=
\int {g(y_n \,|\, x_n)}\frac{1}{N}\sum_{i=1}^N \delta_{\tilde{X}_{1:n-1}}^{(i)}(x_{1:n}) \,\mathrm{d}x_n \\
&=
\frac{1}{N}\sum_{i=1}^N {g(y_n \,|\, \tilde{X}_n^{(i)})}
\end{aligned}
$$

and again

$$
\begin{aligned}
\tilde{p}(x_{1:n} \,|\, y_{1:n}) &=
\frac{g(y_n \,|\, x_{n}) \, {\hat{p}(x_{1:n} \,|\, y_{1:n-1})}}
     {\hat{p}(y_n \,|\, y_{1:n-1})} \\
&=
\frac{g(y_n \,|\, x_{n}) \, \frac{1}{N}\sum_{i=1}^N \delta_{\tilde{X}_{1:n}}^{(i)}(x_{1:n})}
     {\frac{1}{N}\sum_{i=1}^N {g(y_n \,|\, \tilde{X}_n^{(i)})}} \\
&=
\frac{ \sum_{i=1}^N g(y_n \,|\, \tilde{X}_n^{(i)}) \, \delta_{\tilde{X}_{1:n}}^{(i)}(x_{1:n})}
     {\sum_{i=1}^N {g(y_n \,|\, \tilde{X}_n^{(i)})}} \\
&=
\sum_{i=1}^N W_n^{(i)} \delta_{\tilde{X}_{1:n}^{(i)}} (x_{1:n})
\end{aligned}
$$

where $W_n^{(i)} \propto g(y_n \,|\, \tilde{X}_n^{(i)})$ and $\sum_{i=1}^N W_n^{(i)} = 1$.

Now sample

$$
X_{1:n}^{(i)} \sim \tilde{p}(x_{1:n} \,|\, y_{1:n})
$$

A Haskell Implementation
========================

Let's specify some values for the example of the car moving in two
dimensions.

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

> n :: Int
> n = 5

With these we generate hidden and observable sample path.

> carSample :: MonadRandom m =>
>              H.Vector Double ->
>              m (Maybe ((H.Vector Double, H.Vector Double), H.Vector Double))
> carSample xPrev = do
>   xNew <- sample $ rvar (Normal (bigA H.#> xPrev) bigQ)
>   yNew <- sample $ rvar (Normal (bigH H.#> xNew) bigR)
>   return $ Just ((xNew, yNew), xNew)

> carSamples :: [(H.Vector Double, H.Vector Double)]
> carSamples = evalState (ML.unfoldrM carSample m0) (pureMT 17)

We can plot an example trajectory for the car and the noisy
observations that are available to the smoother / filter.

``` {.dia height='600'}
dia = image (DImage (ImageRef "diagrams/CarPosition.png") 600 600 (translationX 0.0))
```

Sadly there is no equivalent to [numpy](http://www.numpy.org/) in
Haskell. Random number packages generate
[vectors](https://hackage.haskell.org/package/vector), for multi-rank
arrays there is [repa](https://hackage.haskell.org/package/repa) and
for fast matrix manipulation there is
[hmtatrix](https://hackage.haskell.org/package/hmatrix). Thus for our
single step path update function, we have to pass in functions to
perform type conversion. Clearly with all the copying inherent in this
approach, performance is not going to be great.

The type synonym *ArraySmoothing* is used to denote the cloud of
particles at each time step.

> type ArraySmoothing = Repa.Array Repa.U DIM2

> singleStep :: forall a . U.Unbox a =>
>               (a -> H.Vector Double) ->
>               (H.Vector Double -> a) ->
>               H.Matrix Double ->
>               H.Herm Double ->
>               H.Matrix Double ->
>               H.Herm Double ->
>               ArraySmoothing a -> H.Vector Double ->
>               WriterT [(ArraySmoothing a, V.Vector Double)] (StateT PureMT IO) (ArraySmoothing a)
> singleStep f g bigA bigQ bigH bigR x y = do
>   let (Z :. ix :. jx) = extent x
>
>   xHatR <- lift $ computeP $ Repa.slice x (Any :. jx - 1)
>   let xHatH = map f $ Repa.toList (xHatR  :: Repa.Array Repa.U DIM1 a)
>   xTildeNextH <- lift $ mapM (\x -> sample $ rvar (Normal (bigA H.#> x) bigQ)) xHatH
>
>   let xTildeNextR = Repa.fromListUnboxed (Z :. ix :. (1 :: Int)) $
>                     map g xTildeNextH
>       xTilde = Repa.append x xTildeNextR
>
>       weights = map (normalPdf y bigR) $
>                 map (bigH H.#>) xTildeNextH
>       vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))
>       totWeight = sum weights
>       normWeights = V.map (/ totWeight) (V.fromList weights)
>       cumSumWeights = V.tail $ V.scanl (+) 0 normWeights
>       js = indices cumSumWeights vs
>       xNewV = V.map (\j -> Repa.transpose $
>                            Repa.reshape (Z :. (1 :: Int) :. jx + 1) $
>                            slice xTilde (Any :. j :. All)) js
>       xNewR = Repa.transpose $ V.foldr Repa.append (xNewV V.! 0) (V.tail xNewV)
>   tell [(x, normWeights)]
>   -- trace ("\nNormalised weights: " ++ show normWeights) $ return ()
>   computeP xNewR


The state for the car is a 4-tuple.

> data SystemState a = SystemState { xPos  :: a
>                                  , yPos  :: a
>                                  , _xSpd :: a
>                                  , _ySpd :: a
>                                  }
>   deriving Show

We initialize the smoother from some prior distribution.

> initCar :: StateT PureMT IO (ArraySmoothing (SystemState Double))
> initCar = do
>   xTilde1 <- replicateM n $ sample $ rvar (Normal m0 bigP0)
>   let weights = map (normalPdf (snd $ head carSamples) bigR) $
>                 map (bigH H.#>) xTilde1
>       vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))
>       cumSumWeights = V.scanl (+) 0 (V.fromList weights)
>       js = indices (V.tail cumSumWeights) vs
>       xHat1 = Repa.fromListUnboxed (Z :. n :. (1 :: Int)) $
>               map ((\[a,b,c,d] -> SystemState a b c d) . H.toList) $
>               V.toList $
>               V.map ((V.fromList xTilde1) V.!) js
>   return xHat1

Now we can run the smoother.

> nCarPositions :: Int
> nCarPositions = 4

> smootherCar :: StateT PureMT IO
>             (ArraySmoothing (SystemState Double)
>             , [(ArraySmoothing (SystemState Double), V.Vector Double)])
> smootherCar = runWriterT $ do
>   xHat1 <- lift initCar
>   foldM (singleStep f h bigA bigQ bigH bigR) xHat1 (take nCarPositions $ map snd $ tail carSamples)

> f :: SystemState Double -> H.Vector Double
> f (SystemState a b c d) = H.fromList [a, b, c, d]

> h :: H.Vector Double -> SystemState Double
> h = (\[a,b,c,d] -> (SystemState a b c d)) . H.toList

And create inferred positions for the car which we then plot.

> testCar :: IO ([Double], [Double])
> testCar = do
>   states <- (map fst . snd) <$> evalStateT smootherCar (pureMT 24)
>   let xs :: [Repa.Array Repa.D DIM2 Double]
>       xs = map (Repa.map xPos) states
>   sumXs :: [Repa.Array Repa.U DIM1 Double] <- mapM Repa.sumP (map Repa.transpose xs)
>   let ixs = map extent sumXs
>       sumLastXs = map (* (recip $ fromIntegral n)) $
>                   zipWith (Repa.!) sumXs (map (\(Z :. x) -> Z :. (x - 1)) ixs)
>   let ys :: [Repa.Array Repa.D DIM2 Double]
>       ys = map (Repa.map yPos) states
>   sumYs :: [Repa.Array Repa.U DIM1 Double] <- mapM Repa.sumP (map Repa.transpose ys)
>   let ixsY = map extent sumYs
>       sumLastYs = map (* (recip $ fromIntegral n)) $
>                   zipWith (Repa.!) sumYs (map (\(Z :. x) -> Z :. (x - 1)) ixsY)
>   return (sumLastXs, sumLastYs)

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/CarPositionInf.png") 600 600 (translationX 0.0))
```

So it seems our smoother does quite well at inferring the state at the
**latest** observation, that is, when it is working as a filter. But
what about estimates for earlier times? We should do better as we have
observations in the past and in the future. Let's try with a simpler
example and a smaller number of particles.

First we create some samples for our simple 1 dimensional linear
Gaussian model.

> bigA1, bigQ1, bigR1, bigH1 :: Double
> bigA1 = 0.5
> bigQ1 = 0.01
> bigR1 = 0.1
> bigH1 = 1.0

> simpleSample :: MonadRandom m =>
>               Double ->
>               m (Maybe ((Double, Double), Double))
> simpleSample xPrev = do
>   xNew <- sample $ rvar (R.Normal (bigA1 * xPrev) bigQ1)
>   yNew <- sample $ rvar (R.Normal (bigH1 * xNew) bigR1)
>   return $ Just ((xNew, yNew), xNew)

> simpleSamples :: [(Double, Double)]
> simpleSamples = evalState (ML.unfoldrM simpleSample 0.0) (pureMT 17)

Again create a prior.

> initSimple :: MonadRandom m => m (ArraySmoothing Double)
> initSimple = do
>   let y = snd $ head simpleSamples
>   xTilde1 <- replicateM n $ sample $ rvar $ R.Normal y bigR1
>   let weights = map (pdf (R.Normal y bigR1)) $
>                 map (bigH1 *) xTilde1
>       totWeight = sum weights
>       vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))
>       cumSumWeights = V.scanl (+) 0 (V.fromList $ map (/ totWeight) weights)
>       js = indices (V.tail cumSumWeights) vs
>       xHat1 = Repa.fromListUnboxed (Z :. n :. (1 :: Int)) $
>               V.toList $
>               V.map ((V.fromList xTilde1) V.!) js
>   return xHat1

Now we can run the smoother.

> nSimples :: Int
> nSimples = 7

> smootherSimple :: StateT PureMT IO (ArraySmoothing Double, [(ArraySmoothing Double, V.Vector Double)])
> smootherSimple = runWriterT $ do
>   xHat1 <- lift initSimple
>   foldM (singleStep f1 g1 ((1 H.>< 1) [bigA1]) (H.trustSym $ (1 H.>< 1) [bigQ1^2])
>                           ((1 H.>< 1) [bigH1]) (H.trustSym $ (1 H.>< 1) [bigR1^2]))
>         xHat1
>         (take nSimples $ map H.fromList $ map return . map snd $ tail simpleSamples)

> f1 :: Double -> H.Vector Double
> f1 a = H.fromList [a]

> g1 :: H.Vector Double -> Double
> g1 = (\[a] -> a) . H.toList

And finally we can look at the paths not just the means of the
marginal distributions at the latest observation time.

> testSimple :: IO [[[Double]]]
> testSimple = do
>   states <- (map fst . snd) <$> evalStateT smootherSimple (pureMT 24)
>   let path :: Int -> Int -> IO (Repa.Array Repa.U DIM1 Double)
>       path i j = computeP $ Repa.slice (states!!j) (Any :. i :. All)
>   pathss <- mapM (\i -> mapM (flip path i) [0..n - 1]) [0..4]
>   return $ map (\paths ->map Repa.toList paths) pathss

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/Smooth.png") 600 600 (translationX 0.0))
```

We can see that at some point in the past all the current particles
have one ancestor. The marginals of the smoothing distribution (at
some point in the past) have collapsed on to one particle.

> testSimple1 :: IO [[Double]]
> testSimple1 = do
>   states <- (map fst . snd) <$> evalStateT smootherSimple (pureMT 24)
>   let ijxs = map extent states
>       f ijx = let Z :. i :. j = ijx in (i,j)
>       ijs = map f ijxs
>   let prePss :: [Repa.Array Repa.D DIM1 Double]
>       prePss :: [Repa.Array Repa.D DIM1 Double] =
>         map (\(_i, j) -> Repa.slice (states!!(j - 1)) (Any :. j - 1)) ijs
>   pss :: [Repa.Array Repa.U DIM1 Double] <- mapM computeP prePss
>   return $ map Repa.toList pss

> testSimple2 = do
>   weightParticles <- snd <$> evalStateT smootherSimple (pureMT 24)
>   let weights :: [V.Vector Double]
>       weights = map snd weightParticles
>       histories = map fst weightParticles
>   let ijxs = map extent histories
>       ijs = map (\ijx -> let Z :. i :. j = ijx in (i,j)) ijxs
>       prePss = map (\(_i, j) -> Repa.slice (histories!!(j - 1)) (Any :. j - 1)) ijs
>   pss :: [Repa.Array Repa.U DIM1 Double] <- mapM computeP prePss
>   return (map Repa.toList pss, map V.toList weights)


> (.+), (.*), (.-) :: Num a => [a] -> [a] -> [a]
> (.+) = zipWith (+)
> (.*) = zipWith (*)
> (.-) = zipWith (-)

> (./) :: Fractional a => [a] -> [a] -> [a]
> (./) = zipWith (/)

> smoothOneStep :: H.Matrix Double ->
>                  H.Herm Double ->
>                  [Double] ->
>                  ([Double], [H.Vector Double], [H.Vector Double]) ->
>                  [Double]
> smoothOneStep bigA bigQ wSuccs (filteringWtCurrs, particleCurrs, particleSuccs) =
>   trace ("\nnormWeights = " ++ show normWeightsList ++
>          "\nlog normWeights = " ++ show logNormWeightsList ++
>          "\nnormTerms = " ++ show normTerms ++
>          "\nIndividual weights = " ++ show indWeightsList) $
>   --        "\nCurrent un-normalised weight = " ++ show wCurrs' ++
>   --        "\nCurrent weight = " ++ show wCurrs) $
>   -- trace ("\nNext weights = " ++ show wSuccs ++
>   --        "\nFiltering weights = " ++ show filteringWtCurrs ++
>   --        "\nCurrent weight = " ++ show wCurrs) $
>   wCurrs
>   where
>     normWeights :: H.Vector Double -> [Double]
>     normWeights particleSucc = map (normalPdf particleSucc bigQ) $
>                                map (bigA #>) particleCurrs
>     logNormWeights particleSucc = map (logPdf (Normal particleSucc bigQ)) $
>                                   map (bigA #>) particleCurrs
>     normWeightsList = map normWeights particleSuccs
>     logNormWeightsList = map logNormWeights particleSuccs
>     normTerm :: H.Vector Double -> Double
>     normTerm particleSucc = sum $ filteringWtCurrs .* (normWeights particleSucc)
>     normTerms :: [Double]
>     normTerms = map normTerm particleSuccs
>     indWeights :: H.Vector Double -> [Double]
>     indWeights particleCurr = map (\x -> normalPdf x bigQ (bigA #> particleCurr))
>                                   particleSuccs
>     indWeightsList = map indWeights particleCurrs
>     unSmoothWeight :: Double -> H.Vector Double -> Double
>     unSmoothWeight filteringWtCurr particleCurr =
>       filteringWtCurr * (sum (wSuccs .* (indWeights particleCurr)))
>     smoothWeight :: Double -> H.Vector Double -> Double
>     smoothWeight filteringWtCurr particleCurr =
>       filteringWtCurr * (sum (wSuccs .* (indWeights particleCurr) ./ normTerms))
>     wCurrs :: [Double]
>     wCurrs = zipWith smoothWeight filteringWtCurrs particleCurrs
>     wCurrs' :: [Double]
>     wCurrs' = zipWith unSmoothWeight filteringWtCurrs particleCurrs

Let's run the filter with

    [ghci]
    n

> foo = do
>   weightParticles <- snd <$> evalStateT smootherSimple (pureMT 24)
>   let weights :: [V.Vector Double]
>       weights = map snd weightParticles
>       particleHistories = map fst weightParticles
>   return (particleHistories, weights)

> baz = do
>   (histories, pss, weights) <- bar
>   mapM_ putStrLn $ map show weights
>   mapM_ putStrLn $ map render $ map pPrint histories
>   mapM_ putStrLn $ map render $ map pPrint pss

> bar = do
>   (histories, weights) <- foo
>   let ijxs = map extent histories
>       ijs = map (\ijx -> let Z :. i :. j = ijx in (i,j)) ijxs
>       prePss = map (\(_i, j) -> Repa.slice (histories!!(j - 1)) (Any :. j - 1)) ijs
>   pss :: [Repa.Array Repa.U DIM1 Double] <- mapM computeP prePss
>   return (histories, pss, weights)

    [ghci]
    baz

> urk = do
>   (_, pss, weights) <- bar
>   let revWeights = reverse weights
>       initSmoothingWeights = V.toList $ head revWeights
>       filteringWeights1 = V.toList $ head $ tail revWeights
>       filteringWeights2 = V.toList $ head $ tail $ tail revWeights
>       zs :: [[H.Vector Double]]
>       zs = map (map (\x -> H.fromList [x]) . Repa.toList) $ reverse pss
>       eek1 = smoothOneStep ((1 H.>< 1) [bigA1]) (H.trustSym $ (1 H.>< 1) [bigQ1]) initSmoothingWeights (filteringWeights1, head $ tail zs, head zs)
>       eek2 = smoothOneStep ((1 H.>< 1) [bigA1]) (H.trustSym $ (1 H.>< 1) [bigQ1]) eek1 (filteringWeights2, head $ tail $ tail zs, head $ tail zs)
>   return [eek1, eek2]

    [ghci]
    urk >>= mapM_ (putStrLn . render. pPrint)

> smoothBackwards ::
>   forall a . (U.Unbox a, Show a) =>
>   (a -> H.Vector Double) ->
>   H.Matrix Double ->
>   H.Herm Double ->
>   (StateT PureMT IO (ArraySmoothing a, [(Repa.Array Repa.U DIM2 a, V.Vector Double)])) ->
>   IO ([[a]], [[Double]])
> smoothBackwards f bigA bigQ filter = do
>   weightParticles <- snd <$> evalStateT filter (pureMT 24)
>   let weights :: [V.Vector Double]
>       weights = map snd weightParticles
>       particleHistories = map fst weightParticles
>   let ijxs = map extent particleHistories
>       ijs = map (\ijx -> let Z :. i :. j = ijx in (i,j)) ijxs
>   -- trace (show ijs) $ return ()
>   let prePss :: [Repa.Array Repa.D DIM1 a]
>       prePss = map (\(_i, j) -> Repa.slice (particleHistories!!(j - 1)) (Any :. j - 1)) ijs
>   pss :: [Repa.Array Repa.U DIM1 a] <- mapM computeP prePss
>   let revWeights = reverse weights
>       initSmoothingWeights = V.toList $ head revWeights
>       tailFilteringWeights = map V.toList $ tail revWeights
>       zs :: [[H.Vector Double]]
>       zs = map (map f . Repa.toList) $ reverse pss
>       xs :: [([Double], [H.Vector Double], [H.Vector Double])]
>       xs = zip3 (tailFilteringWeights) zs (tail zs)
>       smoothingWeights = scanl' (smoothOneStep bigA bigQ) initSmoothingWeights xs
>   -- trace (show weights) $ return ()
>   -- trace (show initSmoothingWeights) $ return ()
>   return $ (map Repa.toList pss, smoothingWeights)

> testSmoothSimple :: IO ([[Double]], [[Double]])
> testSmoothSimple = do
>     smoothBackwards (\x -> H.fromList [x]) ((1 H.>< 1) [bigA1]) (H.trustSym $ (1 H.>< 1) [bigQ1]) smootherSimple

> testSmoothCar :: IO ([[ Double]], [[Double]], [[Double]])
> testSmoothCar = do
>     (xyss, wss) <- smoothBackwards f bigA bigQ smootherCar
>     -- trace (show xyss) $ return ()
>     -- trace (show wss) $ return ()
>     let xss = map (map xPos) xyss
>     let yss = map (map yPos) xyss
>     return (xss, yss, wss)

> fooCar = do
>   weightParticles <- snd <$> evalStateT smootherCar (pureMT 24)
>   let weights :: [V.Vector Double]
>       weights = map snd weightParticles
>       particleHistories = map fst weightParticles
>   return (particleHistories, weights)

> barCar = do
>   (histories, weights) <- fooCar
>   let ijxs = map extent histories
>       ijs = map (\ijx -> let Z :. i :. j = ijx in (i,j)) ijxs
>       prePss = map (\(_i, j) -> Repa.slice (histories!!(j - 1)) (Any :. j - 1)) ijs
>   pss :: [Repa.Array Repa.U DIM1 (SystemState Double)] <- mapM computeP prePss
>   return (histories, pss, weights)

> urkCar = do
>   (_, pss, weights) <- barCar
>   let revWeights = reverse weights
>       initSmoothingWeights = V.toList $ head revWeights
>       filteringWeights1 = V.toList $ head $ tail revWeights
>       filteringWeights2 = V.toList $ head $ tail $ tail revWeights
>       filteringWeights3 = V.toList $ head $ tail $ tail $ tail revWeights
>       zs :: [[H.Vector Double]]
>       zs = map (map f . Repa.toList) $ reverse pss
>       eek1 = smoothOneStep bigA bigQ initSmoothingWeights (filteringWeights1, head $ tail zs, head zs)
>       eek2 = smoothOneStep bigA bigQ eek1 (filteringWeights2, head $ tail $ tail zs, head $ tail zs)
>       eek3 = smoothOneStep bigA bigQ eek2 (filteringWeights3, head $ tail $ tail $ tail zs, head $ tail $ tail zs)
>   return [eek1, eek2, eek3]

Pendulum
========

$$
\frac{\mathrm{d}^2\alpha}{\mathrm{d}t^2} = -g\sin\alpha + w(t)
$$

$$
\begin{bmatrix}
x_{1,i} \\
x_{2,i}
\end{bmatrix}
=
\begin{bmatrix}
x_{1,i-1} + x_{2,i-1}\Delta t \\
x_{2,i-1} - g\sin x_{1,i-1}\Delta t
\end{bmatrix}
+
\mathbf{q}_{i-1}
$$

where $q_i \sim {\mathcal{N}}(0,Q)$

$$
Q
=
Q =
\begin{bmatrix}
\frac{q^c \Delta t^3}{3} & \frac{q^c \Delta t^2}{2} \\
\frac{q^c \Delta t^2}{2} & {q^c \Delta t}
\end{bmatrix}
$$

Assume that we can only measure the horizontal position of the
pendulum so that

$$
y_i = \sin x_i + r_k
$$

where $r_i \sim {\mathcal{N}}(0,R)$.

> deltaT' = 0.01
> g  = 9.81
> bigQ' :: H.Herm Double
> bigQ' = H.trustSym $ (2 H.>< 2) bigQl'
> qc1' = 0.01

> bigQl' :: [Double]
> bigQl' = [qc1' * deltaT'^3 / 3, qc1' * deltaT'^2 / 2,
>           qc1' * deltaT'^2 / 2,       qc1' * deltaT']

> bigR' :: H.Matrix Double
> bigR'  = (1 H.>< 1) [0.1]

> type PendulumState = S.R 2

> pendulumStateUpdate :: MonadRandom m => PendulumState -> m PendulumState
> pendulumStateUpdate xPrev = do
>   let x1Prev = fst $ S.headTail xPrev
>       x2Prev = fst $ S.headTail $ snd $ S.headTail xPrev
>   eta <- sample $ rvar (Normal 0.0 bigQ')
>   let x1New = x1Prev + x2Prev * deltaT'
>       x2New = x2Prev - g * sin (x1Prev * deltaT')
>   undefined



m0 = [1.6;0]; % Slightly off
P0 = 0.1*eye(2);

> pendulumSample :: MonadRandom m =>
>              H.Vector Double ->
>              m (Maybe ((H.Vector Double, H.Vector Double), H.Vector Double))
> pendulumSample xPrev = do
>   xNew <- sample $ rvar (Normal (bigA H.#> xPrev) bigQ)
>   yNew <- sample $ rvar (Normal (bigH H.#> xNew) bigR)
>   return $ Just ((xNew, yNew), xNew)

> pendulumSamples :: [(H.Vector Double, H.Vector Double)]
> pendulumSamples = evalState (ML.unfoldrM pendulumSample m0) (pureMT 17)


Notes
=====

Helpers for the Inverse CDF
---------------------------

That these are helpers for the inverse CDF is delayed to another blog
post.

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

Multivariate Normal
-------------------

The *random-fu* package does not contain a sampler or pdf for a
multivariate normal so we create our own. This should be added to
*random-fu-multivariate* package or something similar.

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

Misellaneous
------------

> derivingUnbox "SystemState"
>     [t| forall a . (U.Unbox a) => SystemState a -> (a, a, a, a) |]
>     [| \ (SystemState x y xdot ydot) -> (x, y, xdot, ydot) |]
>     [| \ (x, y, xdot, ydot) -> SystemState x y xdot ydot |]

> instance Pretty a => Pretty (SystemState a) where
>   pPrint (SystemState x y xdot ydot ) = pPrint (x, y, xdot, ydot)

> liftA4 :: Applicative f => (a -> b -> c -> d -> e) -> f a -> f b -> f c -> f d -> f e
> liftA4 f a b c d = fmap f a <*> b <*> c <*> d

> store :: Storable a => Store.Dictionary (SystemState a)
> store =
>   Store.run $
>   liftA4 SystemState
>   (Store.element xPos)
>   (Store.element yPos)
>   (Store.element _xSpd)
>   (Store.element _ySpd)

> instance (Storable a) => Storable (SystemState a) where
>   sizeOf = Store.sizeOf store
>   alignment = Store.alignment store
>   peek = Store.peek store
>   poke = Store.poke store

Bibliography
============
