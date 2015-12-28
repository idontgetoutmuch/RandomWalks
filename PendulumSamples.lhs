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
> import           Numeric.LinearAlgebra.Static
>                  ( R, vector, Sym,
>                    headTail, matrix, sym
>                  )
> import           GHC.TypeLits

> import           MultivariateNormal

> import qualified Data.Vector.Unboxed as U
> import           Data.Vector.Unboxed.Deriving


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

> data SystemState a = SystemState { x1  :: a, x2  :: a }
>   deriving Show

> newtype SystemObs a = SystemObs { y1  :: a }
>   deriving Show

> type Particles a = U.Vector a

> stateUpdate :: MonadRandom m =>
>           Particles (SystemState Double) ->
>           m (Particles (SystemState Double))
> stateUpdate xPrevs = do
>   let ix = U.length xPrevs
>
>   let x1Prevs = U.map x1 xPrevs
>       x2Prevs = U.map x2 xPrevs
>
>   eta <- replicateM ix $ sample $ rvar (MultivariateNormal 0.0 bigQ')
>
>   let x1 = U.zipWith (\x1Prev x2Prev -> x1Prev + x2Prev * deltaT) x1Prevs x2Prevs
>       x2 = U.zipWith (\x1Prev x2Prev -> x2Prev - g * sin (x1Prev * deltaT)) x1Prevs x2Prevs
>
>   return (U.zipWith SystemState x1 x2)

> obsUpdate :: Particles (SystemState Double) ->
>              Particles (SystemObs Double)
> obsUpdate xs = U.map (SystemObs . sin . x1) xs

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
>   m (Particles (SystemState Double), Particles (SystemObs Double))
> oneFilteringStep stateUpdate obsUpdate weight statePrevs obs = do
>   stateNews <- stateUpdate statePrevs
>   let obsNews = obsUpdate stateNews
>   let weights :: U.Vector Double
>       weights = U.map (weight obs) obsNews
>       totWeight = U.sum weights
>   let normalisedWeights = U.map (/totWeight) weights
>   return (stateNews, obsNews)

>       -- cumSumWeights = V.tail $ V.scanl (+) 0 normWeights

> test :: MonadRandom m =>
>         Particles (SystemState Double) ->
>         SystemObs Double ->
>         m (Particles (SystemState Double), Particles (SystemObs Double))
> test = oneFilteringStep stateUpdate obsUpdate (weight undefined bigR')

> derivingUnbox "SystemState"
>     [t| forall a . U.Unbox a => SystemState a -> (a, a) |]
>     [| \ x -> (x1 x, x2 x) |]
>     [| \ (u, v) -> SystemState u v |]

> derivingUnbox "SystemObs"
>     [t| forall a . U.Unbox a => SystemObs a -> a |]
>     [| \ x -> y1 x |]
>     [| \ u -> SystemObs u |]
