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
>                        , test
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

> import qualified Data.Vector.Unboxed as V
> import           Data.Vector.Unboxed.Deriving

> import           Data.Bits ( shiftR )


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
>       x2 = x2Prev - g * sin (x1Prev * deltaT)
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

> nParticles :: Int
> nParticles = 10

> data SystemState a = SystemState { x1  :: a, x2  :: a }
>   deriving Show

> newtype SystemObs a = SystemObs { y1  :: a }
>   deriving Show

> type Particles a = V.Vector a

> (.+), (.*), (.-) :: (V.Unbox a, Num a) => V.Vector a -> V.Vector a -> V.Vector a
> (.+) = V.zipWith (+)
> (.*) = V.zipWith (*)
> (.-) = V.zipWith (-)

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
>       x1 = x1Prevs .+ (x2Prevs .* deltaTs) .+ eta1s
>       x2 = x2Prevs .- (gs .* V.map sin (x1Prevs .* deltaTs)) .+ eta2s
>
>   return (V.zipWith SystemState x1 x2)

> obsUpdate :: Particles (SystemState Double) ->
>              Particles (SystemObs Double)
> obsUpdate xs = V.map (SystemObs . sin . x1) xs

> f :: SystemObs Double -> R 1
> f = vector . pure . y1

> weight :: forall n . KnownNat n =>
>           (SystemObs Double -> R n) ->
>           Sym n ->
>           SystemObs Double -> SystemObs Double -> Double
> weight f bigR obs obsNew = pdf (MultivariateNormal (f obsNew) bigR) (f obs)

> oneFilteringStep ::
>   MonadRandom m =>
>   (Particles (SystemState Double) -> m (Particles (SystemState Double))) ->
>   (Particles (SystemState Double) -> (Particles (SystemObs Double))) ->
>   (SystemObs Double -> SystemObs Double -> Double) ->
>   Particles (SystemState Double) ->
>   SystemObs Double ->
>   WriterT Log m (Particles (SystemState Double))
> oneFilteringStep stateUpdate obsUpdate weight statePrevs obs = do
>   tell [(V.sum (V.map x1 statePrevs) / fromIntegral nParticles,
>          V.sum (V.map x2 statePrevs) / fromIntegral nParticles)]
>   stateNews <- lift $ stateUpdate statePrevs
>   let obsNews = obsUpdate stateNews
>   let weights       = V.map (weight obs) obsNews
>       cumSumWeights = V.tail $ V.scanl (+) 0 weights
>       totWeight     = V.last cumSumWeights
>   vs <- lift $ V.replicateM nParticles (sample $ uniform 0.0 totWeight)
>   let js = indices cumSumWeights vs
>       stateTildes = V.map (stateNews V.!) js
>   return stateTildes

> bigP :: Sym 2
> bigP = sym $ diag 0.1

> initParticles :: MonadRandom m =>
>                  m (Particles (SystemState Double))
> initParticles = V.replicateM nParticles $ do
>   r <- sample $ rvar (MultivariateNormal 0.0 bigP)
>   let x1 = fst $ headTail r
>       x2 = fst $ headTail $ snd $ headTail r
>   return $ SystemState { x1 = x1, x2 = x2}

> type Log = [(Double, Double)]

> test :: Int -> Log
> test n = snd $ evalState action (pureMT 19)
>   where
>     action = runWriterT $ do
>       xs <- lift $ initParticles
>       V.foldM (oneFilteringStep (stateUpdate bigQ') obsUpdate (weight f bigR'))
>             xs
>             (V.fromList $ map (SystemObs . fst . headTail . snd) (take n pendulumSamples'))

Notes
=====

Helpers for the Inverse CDF
---------------------------

That these are helpers for the inverse CDF is delayed to another blog
post.

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

Other
-----

> derivingUnbox "SystemState"
>     [t| forall a . V.Unbox a => SystemState a -> (a, a) |]
>     [| \ x -> (x1 x, x2 x) |]
>     [| \ (u, v) -> SystemState u v |]

> derivingUnbox "SystemObs"
>     [t| forall a . V.Unbox a => SystemObs a -> a |]
>     [| \ x -> y1 x |]
>     [| \ u -> SystemObs u |]
