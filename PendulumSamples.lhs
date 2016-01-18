% Particle Smoothing
% Dominic Steinitz
% 18th January 2016

---
bibliography: Stochastic.bib
---

Introduction
============

The [equation of
motion](https://en.wikipedia.org/wiki/Pendulum_(mathematics)#Simple_gravity_pendulum)
for a pendulum of unit length subject to [Gaussian white
noise](https://en.wikipedia.org/wiki/White_noise#Mathematical_definitions)
is

$$
\frac{\mathrm{d}^2\alpha}{\mathrm{d}t^2} = -g\sin\alpha + w(t)
$$

We can discretize this via the usual [Euler method](https://en.wikipedia.org/wiki/Euler_method)

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

where $q_i \sim {\mathcal{N}}(0,Q)$ and

$$
Q
=
\begin{bmatrix}
\frac{q^c \Delta t^3}{3} & \frac{q^c \Delta t^2}{2} \\
\frac{q^c \Delta t^2}{2} & {q^c \Delta t}
\end{bmatrix}
$$

The explanation of the precise form of the covariance matrix will be
the subject of another blog post; for the purpose of exposition of
forward filtering / backward smoothing, this detail is not important.

Assume that we can only measure the horizontal position of the
pendulum and further that this measurement is subject to error so that

$$
y_i = \sin x_i + r_k
$$

where $r_i \sim {\mathcal{N}}(0,R)$.

Particle Filtering can give us an estimate of where the pendulum is
and its velocity using all the observations up to that point in
time. But now suppose we have observed the pendulum for a fixed period
of time then at times earlier than the time at which we stop our
observations we now have observations in the future as well as in the
past. If we can somehow take account of these future observations then
we should be able to improve our estimate of where the pendulum was at
any given point in time (as well as its velocity). Forward Filtering /
Backward Smoothing is a technique for doing this.

Haskell Preamble
----------------

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}

> {-# LANGUAGE MultiParamTypeClasses        #-}
> {-# LANGUAGE TypeFamilies                 #-}
> {-# LANGUAGE ScopedTypeVariables          #-}
> {-# LANGUAGE DataKinds                    #-}


> {-# LANGUAGE FlexibleInstances            #-}
> {-# LANGUAGE MultiParamTypeClasses        #-}
> {-# LANGUAGE FlexibleContexts             #-}
> {-# LANGUAGE TypeFamilies                 #-}
> {-# LANGUAGE BangPatterns                 #-}
> {-# LANGUAGE GeneralizedNewtypeDeriving   #-}
> {-# LANGUAGE ScopedTypeVariables          #-}
> {-# LANGUAGE TemplateHaskell              #-}
> {-# LANGUAGE DataKinds                    #-}

> module PendulumSamples ( pendulumSamples
>                        , pendulumSamples'
>                        , testFiltering
>                        , testSmoothing
>                        , testFilteringG
>                        , testSmoothingG
>                        ) where

> import           Data.Random hiding ( StdNormal, Normal )
> import           Data.Random.Source.PureMT ( pureMT )
> import           Control.Monad.State ( evalState, replicateM )
> import qualified Control.Monad.Loops as ML
> import           Control.Monad.Writer ( tell, WriterT, lift,
>                                         runWriterT
>                                       )
> import           Numeric.LinearAlgebra.Static
>                  ( R, vector, Sym,
>                    headTail, matrix, sym,
>                    diag
>                  )
> import           GHC.TypeLits

> import           MultivariateNormal

> import qualified Data.Vector as V

> import           Data.Bits ( shiftR )

> import           Data.List ( transpose )

Simulation
==========

Let's first plot some typical trajectories of the pendulum.

> deltaT, g :: Double
> deltaT = 0.01
> g  = 9.81

> type PendulumState = R 2
> type PendulumObs = R 1

> pendulumSample :: MonadRandom m =>
>                   Sym 2 ->
>                   Sym 1 ->
>                   PendulumState ->
>                   m (Maybe ((PendulumState, PendulumObs), PendulumState))
> pendulumSample bigQ bigR xPrev = do
>   let x1Prev = fst $ headTail xPrev
>       x2Prev = fst $ headTail $ snd $ headTail xPrev
>   eta <- sample $ rvar (MultivariateNormal 0.0 bigQ)
>   let x1= x1Prev + x2Prev * deltaT
>       x2 = x2Prev - g * (sin x1Prev) * deltaT
>       xNew = vector [x1, x2] + eta
>       x1New = fst $ headTail xNew
>   epsilon <-  sample $ rvar (MultivariateNormal 0.0 bigR)
>   let yNew = vector [sin x1New] + epsilon
>   return $ Just ((xNew, yNew), xNew)

Let's try plotting some samples when we are in the linear region with
which we are familiar from school $\sin\alpha \approx \alpha$.

$$
\frac{\mathrm{d}^2\alpha}{\mathrm{d}t^2} = -g\alpha + w(t)
$$

In this case we expect the horizontal displacement to be approximately
equal to the angle of displacement and thus the observations to be
symmetric about the actuals.

> bigQ :: Sym 2
> bigQ = sym $ matrix bigQl

> qc1 :: Double
> qc1 = 0.0001

> bigQl :: [Double]
> bigQl = [ qc1 * deltaT^3 / 3, qc1 * deltaT^2 / 2,
>           qc1 * deltaT^2 / 2,       qc1 * deltaT
>         ]

> bigR :: Sym 1
> bigR  = sym $ matrix [0.0001]

> m0 :: PendulumState
> m0 = vector [0.01, 0]

> pendulumSamples :: [(PendulumState, PendulumObs)]
> pendulumSamples = evalState (ML.unfoldrM (pendulumSample bigQ bigR) m0) (pureMT 17)

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/PendulumObs.png") 600 600 (translationX 0.0))
```

But if we work in a region in which linearity breaks down then the
observations are no longer symmetrical about the actuals.

> bigQ' :: Sym 2
> bigQ' = sym $ matrix bigQl'

> qc1' :: Double
> qc1' = 0.01

> bigQl' :: [Double]
> bigQl' = [ qc1' * deltaT^3 / 3, qc1' * deltaT^2 / 2,
>            qc1' * deltaT^2 / 2,       qc1' * deltaT
>          ]

> bigR' :: Sym 1
> bigR'  = sym $ matrix [0.1]

> m0' :: PendulumState
> m0' = vector [1.6, 0]

> pendulumSamples' :: [(PendulumState, PendulumObs)]
> pendulumSamples' = evalState (ML.unfoldrM (pendulumSample bigQ' bigR') m0') (pureMT 17)

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/PendulumObs1.png") 600 600 (translationX 0.0))
```

Filtering
=========

> nParticles :: Int
> nParticles = 1000

The usual Bayesian update step.

> type Particles a = V.Vector a

> oneFilteringStep ::
>   MonadRandom m =>
>   (Particles a -> m (Particles a)) ->
>   (Particles a -> Particles b) ->
>   (b -> b -> Double) ->
>   Particles a ->
>   b ->
>   WriterT [Particles a] m (Particles a)
> oneFilteringStep stateUpdate obsUpdate weight statePrevs obs = do
>   tell [statePrevs]
>   stateNews <- lift $ stateUpdate statePrevs
>   let obsNews = obsUpdate stateNews
>   let weights       = V.map (weight obs) obsNews
>       cumSumWeights = V.tail $ V.scanl (+) 0 weights
>       totWeight     = V.last cumSumWeights
>   vs <- lift $ V.replicateM nParticles (sample $ uniform 0.0 totWeight)
>   let js = indices cumSumWeights vs
>       stateTildes = V.map (stateNews V.!) js
>   return stateTildes

The system state and observable.

> data SystemState a = SystemState { x1  :: a, x2  :: a }
>   deriving Show

> newtype SystemObs a = SystemObs { y1  :: a }
>   deriving Show

To make the system state update a bit more readable, let's introduce
some lifted arithmetic operators.

> (.+), (.*), (.-) :: (Num a) => V.Vector a -> V.Vector a -> V.Vector a
> (.+) = V.zipWith (+)
> (.*) = V.zipWith (*)
> (.-) = V.zipWith (-)

The state update itself.

> stateUpdate :: MonadRandom m =>
>                Sym 2 ->
>                Particles (SystemState Double) ->
>                m (Particles (SystemState Double))
> stateUpdate bigQ xPrevs = do
>   let ix = V.length xPrevs
>
>   let x1Prevs = V.map x1 xPrevs
>       x2Prevs = V.map x2 xPrevs
>
>   etas <- replicateM ix $ sample $ rvar (MultivariateNormal 0.0 bigQ)
>   let eta1s, eta2s :: V.Vector Double
>       eta1s = V.fromList $ map (fst . headTail) etas
>       eta2s = V.fromList $ map (fst . headTail . snd . headTail) etas
>
>   let deltaTs = V.replicate ix deltaT
>       gs = V.replicate ix g
>       x1s = x1Prevs .+ (x2Prevs .* deltaTs) .+ eta1s
>       x2s = x2Prevs .- (gs .* (V.map sin x1Prevs) .* deltaTs) .+ eta2s
>
>   return (V.zipWith SystemState x1s x2s)

The function which maps the state to the observable.

> obsUpdate :: Particles (SystemState Double) ->
>              Particles (SystemObs Double)
> obsUpdate xs = V.map (SystemObs . sin . x1) xs

And finally a function to calculate the weight of each particle given
an observation.

> weight :: forall a n . KnownNat n =>
>           (a -> R n) ->
>           Sym n ->
>           a -> a -> Double
> weight f bigR obs obsNew = pdf (MultivariateNormal (f obsNew) bigR) (f obs)

The variance of the prior.

> bigP :: Sym 2
> bigP = sym $ diag 0.1

And we can generate our ensemble of particles chosen from the prior.

> initParticles :: MonadRandom m =>
>                  m (Particles (SystemState Double))
> initParticles = V.replicateM nParticles $ do
>   r <- sample $ rvar (MultivariateNormal m0' bigP)
>   let x1 = fst $ headTail r
>       x2 = fst $ headTail $ snd $ headTail r
>   return $ SystemState { x1 = x1, x2 = x2}

Now we can run the particle filter

> test :: Int -> [Particles (SystemState Double)]
> test nTimeSteps = snd $ evalState action (pureMT 19)
>   where
>     action = runWriterT $ do
>       xs <- lift $ initParticles
>       V.foldM (oneFilteringStep (stateUpdate bigQ') obsUpdate (weight f bigR'))
>             xs
>             (V.fromList $ map (SystemObs . fst . headTail . snd) (take nTimeSteps pendulumSamples'))

and extract the estimated position of the filter.

> testFiltering :: Int -> [Double]
> testFiltering nTimeSteps = map ((/ (fromIntegral nParticles)). sum . V.map x1) (test nTimeSteps)

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/PendulumFitted.png") 600 600 (translationX 0.0))
```

Smoothing
=========

If we could calculate the marginal smoothing distributions $\{p(x_t
\,|\, y_{1:T})\}_{i=1}^T$ then we might be able to sample from
them. Using the Markov assumption of our model that $x_i$ is
independent of $y_{i+1:N}$ given $x_{i+1}$, we have

$$
\begin{aligned}
\overbrace{p(x_i \,|\, y_{1:N})}^{\mathrm{smoother}\,\mathrm{at}\, i} &=
\int p(x_i, x_{i+1} \,|\, y_{1:N}) \,\mathrm{d}x_{i+1} & \text{marginal distribution} \\
&=
\int p(x_{i+1} \,|\, y_{1:N}) \,p(x_{i} \,|\, y_{1:N}, x_{i+1}) \,\mathrm{d}x_{i+1} & \text{conditional density} \\
&=
\int p(x_{i+1} \,|\, y_{1:N}) \,p(x_{i} \,|\, y_{1:i}, x_{i+1}) \,\mathrm{d}x_{i+1} & \text{Markov model} \\
&=
\int p(x_{i+1} \,|\, y_{1:N}) \,
\frac{p(x_{i}, x_{i+1} \,|\, y_{1:i})}{p(x_{i+1} \,|\, y_{1:i})}
\,\mathrm{d}x_{i+1}
& \text{conditional density} \\
&=
\int p(x_{i+1} \,|\, y_{1:N}) \,
\frac{p(x_{i+1} \,|\, x_{i}, y_{1:i})
\,p(x_{i} \,|\, y_{1:i})}{p(x_{i+1} \,|\, y_{1:i})}
\,\mathrm{d}x_{i+1}
& \text{conditional density} \\
&=
\int \overbrace{p(x_{i+1} \,|\, y_{1:N})}^{\text{smoother at }i+1} \,
\underbrace{
\overbrace{p(x_{i} \,|\, y_{1:i})}^{\text{filter at }i}
\frac{p(x_{i+1} \,|\, x_{i})}
     {p(x_{i+1} \,|\, y_{1:i})}
}
_{\text{backward transition }p(x_{i} \,|\, y_{1:i},\,x_{i+1})}
\,\mathrm{d}x_{i+1}
& \text{Markov model}
\end{aligned}
$$

We observe that this is a (continuous state space) Markov process with
a non-homogeneous transition function albeit one which goes backwards
in time. Apparently for conditionally Gaussian linear state-space
models, this is known as RTS, or Rauch-Tung-Striebel smoothing
(@RTS-1965).

According to @Cappe:IntroToSMC:Online,

* It appears to be efficient and stable in the long term (although no
proof was available at the time the slides were presented).

* It is not sequential (in particular, one needs to store all particle
positions and weights).

* It has numerical complexity proportional $O(n^2)$ where $N$ is the
number of particles.

We can use this to sample paths starting at time $i = N$ and working
backwards. From above we have

$$
p(x_i \,|\, X_{i+1}, Y_{1:N}) =
\frac{p(X_{i+1} \,|\, x_{i})
\,p(x_{i} \,|\, Y_{1:i})}{p(X_{i+1} \,|\, Y_{1:i})}
=
Z
\,p(X_{i+1} \,|\, x_{i})
\,p(x_{i} \,|\, Y_{1:i})
$$

where $Z$ is some normalisation constant (Z for "Zustandssumme" - sum
over states).

From particle filtering we know that

$$
{p}(x_i \,|\, y_{1:i}) \approx \hat{p}(x_i \,|\, y_{1:i}) \triangleq \sum_{j=1}^M w_i^{(j)}\delta(x_i - x_i^{(j)})
$$

Thus

$$
\hat{p}(x_i \,|\, X_{i+1}, Y_{1:i})
=
Z
\,p(X_{i+1} \,|\, x_{i})
\,\hat{p}(x_{i} \,|\, Y_{1:i})
=
\sum_{j=1}^M w_i^{(j)}\delta(x_i - x_i^{(j)})
\,p(X_{i+1} \,|\, x_{i})
$$

and we can sample $x_i$ from $\{x_i^{(j)}\}_{1 \leq j \leq M}$ with
probability

$$
\frac{w_k^{(i)}
\,p(X_{i+1} \,|\, x_{i})}
{\sum_{i=1}^N w_k^{(i)}
\,p(X_{i+1} \,|\, x_{i})}
$$

Recalling that we have re-sampled the particles in the filtering
algorithm so that their weights are all $1/M$ and abstracting the
state update and state density function, we can encode this update
step in Haskell as

> oneSmoothingStep :: MonadRandom m =>
>          (Particles a -> V.Vector a) ->
>          (a -> a -> Double) ->
>          a ->
>          Particles a ->
>          WriterT (Particles a) m a
> oneSmoothingStep stateUpdate
>                  stateDensity
>                  smoothingSampleAtiPlus1
>                  filterSamplesAti = do it
>   where
>     it = do
>       let mus = stateUpdate filterSamplesAti
>           weights =  V.map (stateDensity smoothingSampleAtiPlus1) mus
>           cumSumWeights = V.tail $ V.scanl (+) 0 weights
>           totWeight     = V.last cumSumWeights
>       v <- lift $ sample $ uniform 0.0 totWeight
>       let ix = binarySearch cumSumWeights v
>           xnNew = filterSamplesAti V.! ix
>       tell $ V.singleton xnNew
>       return $ xnNew

To sample a complete path we start with a sample from the filtering
distribution at at time $i = N$ (which is also the smoothing
distribution)

> oneSmoothingPath :: MonadRandom m =>
>              (Int -> V.Vector (Particles a)) ->
>              (a -> Particles a -> WriterT (Particles a) m a) ->
>              Int -> m (a, V.Vector a)
> oneSmoothingPath filterEstss oneSmoothingStep nTimeSteps = do
>   let ys = filterEstss nTimeSteps
>   ix <- sample $ uniform 0 (nParticles - 1)
>   let xn = (V.head ys) V.! ix
>   runWriterT $ V.foldM oneSmoothingStep xn ys

Of course we need to run through the filtering distributions starting
at $i = N$

> filterEstss :: Int -> V.Vector (Particles (SystemState Double))
> filterEstss n = V.reverse $ V.fromList $ test n



> stateUpdate' :: Particles (SystemState Double) ->
>                 Particles (SystemState Double)
> stateUpdate' xPrevs = V.zipWith SystemState x1s x2s
>   where
>     ix = V.length xPrevs
>
>     x1Prevs = V.map x1 xPrevs
>     x2Prevs = V.map x2 xPrevs
>
>     deltaTs = V.replicate ix deltaT
>     gs = V.replicate ix g
>     x1s = x1Prevs .+ (x2Prevs .* deltaTs)
>     x2s = x2Prevs .- (gs .* (V.map sin x1Prevs) .* deltaTs)

> testSmoothing :: Int -> Int -> [Double]
> testSmoothing m n = V.toList $ evalState action (pureMT 23)
>   where
>     action = do
>       xss <- V.replicateM n $
>              snd <$> (oneSmoothingPath filterEstss (oneSmoothingStep stateUpdate' (weight h bigQ')) m)
>       let yss = V.fromList $ map V.fromList $
>                 transpose $
>                 V.toList $ V.map (V.toList) $
>                 xss
>       return $ V.map (/ (fromIntegral n)) $ V.map V.sum $ V.map (V.map x1) yss

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/PendulumSmoothed20.png") 600 600 (translationX 0.0))
```

Unknown Gravity
===============

Let us continue with the same example but now assume that $g$ is
unknown and that we wish to estimate it. Let us also assume that our
apparatus is not subject to noise.

Again we have

$$
\frac{\mathrm{d}^2\alpha}{\mathrm{d}t^2} = -g\sin\alpha + w(t)
$$

But now when we discretize it we include a third variable

$$
\begin{bmatrix}
x_{1,i} \\
x_{2,i} \\
x_{3,i}
\end{bmatrix}
=
\begin{bmatrix}
x_{1,i-1} + x_{2,i-1}\Delta t \\
x_{2,i-1} - x_{3,i-1}\sin x_{1,i-1}\Delta t \\
x_{3,i-1}
\end{bmatrix}
+
\mathbf{q}_{i-1}
$$

where $q_i \sim {\mathcal{N}}(0,Q)$

$$
Q
=
\begin{bmatrix}
0 & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & \sigma^2_g
\end{bmatrix}
$$

Again we assume that we can only measure the horizontal position of
the pendulum so that

$$
y_i = \sin x_i + r_k
$$

where $r_i \sim {\mathcal{N}}(0,R)$.


> data SystemStateG a = SystemStateG { gx1  :: a, gx2  :: a, gx3 :: a }
>   deriving Show

> stateUpdateG :: MonadRandom m =>
>                Sym 3 ->
>                Particles (SystemStateG Double) ->
>                m (Particles (SystemStateG Double))
> stateUpdateG bigQ xPrevs = do
>   let ix = V.length xPrevs
>
>   let x1Prevs = V.map gx1 xPrevs
>       x2Prevs = V.map gx2 xPrevs
>       x3Prevs = V.map gx3 xPrevs
>
>   etas <- replicateM ix $ sample $ rvar (MultivariateNormal 0.0 bigQ)
>   let eta1s, eta2s, eta3s :: V.Vector Double
>       eta1s = V.fromList $ map (fst . headTail) etas
>       eta2s = V.fromList $ map (fst . headTail . snd . headTail) etas
>       eta3s = V.fromList $ map (fst . headTail . snd . headTail . snd . headTail) etas
>
>   let deltaTs = V.replicate ix deltaT
>       x1s = x1Prevs .+ (x2Prevs .* deltaTs) .+ eta1s
>       x2s = x2Prevs .- (x3Prevs .* (V.map sin x1Prevs) .* deltaTs) .+ eta2s
>       x3s = x3Prevs .+ eta3s
>
>   return (V.zipWith3 SystemStateG x1s x2s x3s)

> stateUpdateG' :: Particles (SystemStateG Double) ->
>                  Particles (SystemStateG Double)
> stateUpdateG' xPrevs = V.zipWith3 SystemStateG x1s x2s x3s
>   where
>     ix = V.length xPrevs
>
>     x1Prevs = V.map gx1 xPrevs
>     x2Prevs = V.map gx2 xPrevs
>     x3Prevs = V.map gx3 xPrevs
>
>     deltaTs = V.replicate ix deltaT
>     x1s = x1Prevs .+ (x2Prevs .* deltaTs)
>     x2s = x2Prevs .- (x3Prevs .* (V.map sin x1Prevs) .* deltaTs)
>     x3s = x3Prevs

> mG :: R 3
> mG = vector [1.6, 0.0, 8.00]

> bigPg :: Sym 3
> bigPg = sym $ matrix [
>     1e-6, 0.0, 0.0
>   , 0.0, 1e-6, 0.0
>   , 0.0, 0.0, 1e-2
>   ]

> bigQg :: Sym 3
> bigQg = sym $ matrix bigQgl

> qc1G :: Double
> qc1G = 0.01

> sigmaG :: Double
> sigmaG = 5.0e-3

> bigQgl :: [Double]
> bigQgl = [ qc1G * deltaT^3 / 3, qc1G * deltaT^2 / 2, 0.0,
>            qc1G * deltaT^2 / 2,       qc1G * deltaT, 0.0,
>                            0.0,                 0.0, sigmaG
>          ]

> bigRg :: Sym 1
> bigRg  = sym $ matrix [0.1]

> initParticlesG :: MonadRandom m =>
>                  m (Particles (SystemStateG Double))
> initParticlesG = V.replicateM nParticles $ do
>   r <- sample $ rvar (MultivariateNormal mG bigPg)
>   let x1 = fst $ headTail r
>       x2 = fst $ headTail $ snd $ headTail r
>       x3 = fst $ headTail $ snd $ headTail $ snd $ headTail r
>   return $ SystemStateG { gx1 = x1, gx2 = x2, gx3 = x3}



> obsUpdateG :: Particles (SystemStateG Double) ->
>              Particles (SystemObs Double)
> obsUpdateG xs = V.map (SystemObs . sin . gx1) xs

> type LogG = [Particles (SystemStateG Double)]

> testG :: Int -> LogG
> testG n = snd $ evalState action (pureMT 19)
>   where
>     action = runWriterT $ do
>       xs <- lift $ initParticlesG
>       V.foldM (oneFilteringStep (stateUpdateG bigQg) obsUpdateG (weight f bigRg))
>             xs
>             (V.fromList $ map (SystemObs . fst . headTail . snd) (take n pendulumSamples'))

> testFilteringG :: Int -> [Double]
> testFilteringG n = map ((/ (fromIntegral nParticles)). sum . V.map gx3) (testG n)

> filterGEstss :: Int -> V.Vector (Particles (SystemStateG Double))
> filterGEstss n = V.reverse $ V.fromList $ testG n

> hG :: SystemStateG Double -> R 3
> hG u = vector [gx1 u , gx2 u, gx3 u]

> testSmoothingG :: Int -> Int -> ([Double], [Double], [Double])
> testSmoothingG m n = (\(x, y, z) -> (V.toList x, V.toList y, V.toList z))  $
>                      evalState action (pureMT 23)
>   where
>     action = do
>       xss <- V.replicateM n $
>              snd <$> (oneSmoothingPath filterGEstss (oneSmoothingStep stateUpdateG' (weight hG bigQg)) m)
>       let yss = V.fromList $ map V.fromList $
>                 transpose $
>                 V.toList $ V.map (V.toList) $
>                 xss
>       return ( V.map (/ (fromIntegral n)) $ V.map V.sum $ V.map (V.map gx1) yss
>              , V.map (/ (fromIntegral n)) $ V.map V.sum $ V.map (V.map gx2) yss
>              , V.map (/ (fromIntegral n)) $ V.map V.sum $ V.map (V.map gx3) yss
>              )

Notes
=====

Helpers for Converting Types
----------------------------

> f :: SystemObs Double -> R 1
> f = vector . pure . y1

> h :: SystemState Double -> R 2
> h u = vector [x1 u , x2 u]

Helpers for the Inverse CDF
---------------------------

That these are helpers for the inverse CDF is delayed to another blog
post.

> indices :: V.Vector Double -> V.Vector Double -> V.Vector Int
> indices bs xs = V.map (binarySearch bs) xs

> binarySearch :: (Ord a) =>
>                 V.Vector a -> a -> Int
> binarySearch vec x = loop 0 (V.length vec - 1)
>   where
>     loop !l !u
>       | u <= l    = l
>       | otherwise = let e = vec V.! k in if x <= e then loop l k else loop (k+1) u
>       where k = l + (u - l) `shiftR` 1

Bibliography
============