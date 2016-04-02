% Modelling an Epidemic
% Dominic Steinitz
% 28th March 2016

---
bibliography: Stochastic.bib
---

Introduction
============

@RSPSA1927:115 give a simple model of the spread of an infectious
disease. Individuals move from being susceptible ($S$) to infected
($I$) to recovered ($R$).

$$
\begin{eqnarray}
\frac{dS}{dt} & = & - \beta S(t) I(t) \label{eq2a} \\
\frac{dI}{dt} & = & \beta S(t) I(t) - k I(t) \label{eq2b} \\
\frac{dR}{dt} & = & k I(t) . \label{eq2c}
\end{eqnarray}
$$

![](diagrams/Sir.png)

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}

> {-# LANGUAGE BangPatterns                 #-}

> module OdeExample where

> import Numeric.GSL.ODE
> import Numeric.LinearAlgebra

> import           Data.Random hiding ( StdNormal, Normal, gamma )
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
> import           GHC.TypeLits ( KnownNat )
> import           MultivariateNormal ( MultivariateNormal(..) )
> import qualified Data.Vector as V
> import qualified Data.Vector.Storable as VS
> import           Data.Bits ( shiftR )
> import           Data.List ( transpose )
> import           Control.Parallel.Strategies
> import           GHC.Generics (Generic)

> beta, gamma :: Double
> beta = 0.0026
> gamma = 0.5

> nParticles :: Int
> nParticles = 500

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

> data SystemState a = SystemState { stateS     :: a
>                                  , stateI     :: a
>                                  , stateR     :: a
>                                  , stateBeta  :: a
>                                  , stateGamma :: a
>                                  }
>   deriving Show

> type Time = Double

> deltaT :: Double
> deltaT = 1.0

> stateUpdate :: (Time, Particles (SystemState Double)) ->
>                (Time, Particles (SystemState Double))
> stateUpdate txPrevs =
>   (newT, V.zipWith5 SystemState ss is rs betas gammas)
>   where
>     absT   = fst txPrevs
>     xPrevs = snd txPrevs
>
>     sPrevs     = V.map stateS     xPrevs
>     iPrevs     = V.map stateI     xPrevs
>     rPrevs     = V.map stateR     xPrevs
>     betaPrevs  = V.map stateBeta  xPrevs
>     gammaPrevs = V.map stateGamma xPrevs
>
>     newT = absT + deltaT
>
>     f xs = odeSolve sirOde xs
>            (absT `VS.cons` (newT `VS.cons` VS.empty))
>
>     ms :: V.Vector (Matrix Double)
>     ms = V.map f (V.zipWith3 (\s i r -> [s, i, r]) sPrevs iPrevs rPrevs)
>
>     ns :: V.Vector (VS.Vector Double)
>     ns = V.map (\m -> (toRows m)!!1) ms
>
>     ss     = V.map (VS.! 0) ns
>     is     = V.map (VS.! 1) ns
>     rs     = V.map (VS.! 2) ns
>     betas  = betaPrevs
>     gammas = gammaPrevs


> sirOde :: Double -> [Double] -> [Double]
> sirOde _t [s, i, _r] =
>   [
>     negate (beta * i * s)
>   , (beta * i * s) - (gamma * i)
>   , gamma * i
>   ]
> sirOde _t vars = error $ "sirOde called with: " ++ show (length vars)


Notes
=====

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

Vector Helpers
--------------

> chunksOf :: Int -> V.Vector a -> V.Vector (V.Vector a)
> chunksOf n xs = ys
>   where
>     l = V.length xs
>     m  = 1 + (l - 1) `div` n
>     ys = V.unfoldrN m (\us -> Just (V.take n us, V.drop n us)) xs


Bibliography
============