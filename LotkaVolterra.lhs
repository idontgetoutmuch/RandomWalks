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
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}

> {-# LANGUAGE BangPatterns                 #-}
> {-# LANGUAGE DataKinds                    #-}
> {-# LANGUAGE ExplicitForAll               #-}

> module LotkaVolterra where

> import Numeric.GSL.ODE
> import Numeric.LinearAlgebra hiding ( R, vector, matrix, sym )

> import           Data.Random hiding ( StdNormal, Normal, gamma )
> import           Data.Random.Source.PureMT ( pureMT )
> import           Control.Monad.State ( evalState )
> import           Control.Monad.Writer ( tell, WriterT, lift,
>                                         runWriterT
>                                       )
> import           Numeric.LinearAlgebra.Static
>                  ( R, vector, Sym,
>                    headTail, matrix, sym
>                  )
> import           GHC.TypeLits ( KnownNat )
> import           MultivariateNormal ( MultivariateNormal(..) )
> import qualified Data.Vector as V
> import qualified Data.Vector.Storable as VS
> import           Data.Bits ( shiftR )
> import           Control.Parallel.Strategies


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

> sol :: Matrix Double
> sol = odeSolve (sirOde delta gamma) [initS, initI, initR] (fromList [0.0,deltaT..14.0])

> solLv :: Matrix Double
> solLv = odeSolve (lvOde a1 a2 b1 b2) [50.0, 50.0] (fromList [0.0,0.1..50])

> nParticles :: Int
> nParticles = 10000

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
>                                  , stateDelta  :: a
>                                  , stateGamma :: a
>                                  }
>   deriving Show

> deltaT :: Double
> deltaT = 1.0

> stateUpdate :: Particles (SystemState Double) ->
>                Particles (SystemState Double)
> stateUpdate xPrevs =
>   V.zipWith5 SystemState ss is rs deltas gammas
>   where
>     sPrevs     = V.map stateS     xPrevs
>     iPrevs     = V.map stateI     xPrevs
>     rPrevs     = V.map stateR     xPrevs
>     deltaPrevs  = V.map stateDelta  xPrevs
>     gammaPrevs = V.map stateGamma xPrevs
>
>     f b g xs = odeSolve (sirOde b g) xs
>                (0.0 `VS.cons` (deltaT `VS.cons` VS.empty))
>
>     ms :: V.Vector (Matrix Double)
>     ms = V.zipWith3 f deltaPrevs gammaPrevs
>                       (V.zipWith3 (\s i r -> [s, i, r]) sPrevs iPrevs rPrevs)
>
>     ns :: V.Vector (VS.Vector Double)
>     ns = V.map (\m -> (toRows m)!!1) ms
>
>     ss     = V.map (VS.! 0) ns
>     is     = V.map (VS.! 1) ns
>     rs     = V.map (VS.! 2) ns
>     deltas  = deltaPrevs
>     gammas = gammaPrevs



> (.+) :: (Num a) => V.Vector a -> V.Vector a -> V.Vector a
> (.+) = V.zipWith (+)

> stateUpdateNoisy :: MonadRandom m =>
>                     Sym 5 ->
>                     Particles (SystemState Double) ->
>                     m (Particles (SystemState Double))
> stateUpdateNoisy bigQ xPrevs = do
>   let xs = stateUpdate xPrevs
>
>       sPrevs :: V.Vector Double
>       sPrevs     = V.map stateS     xs
>       iPrevs     = V.map stateI     xs
>       rPrevs     = V.map stateR     xs
>       deltaPrevs = V.map stateDelta xs
>       gammaPrevs = V.map stateGamma xs
>
>   let mus :: V.Vector (R 5)
>       mus = V.zipWith5 (\a b c d e -> vector (map log [a, b, c, d, e]))
>                        sPrevs iPrevs rPrevs deltaPrevs gammaPrevs
>
>   nus <- mapM (\mu -> fmap exp $ sample $ rvar (MultivariateNormal mu bigQ)) mus
>
>   let nu1s, nu2s, nu3s, nu4s, nu5s :: V.Vector Double
>       nu1s = V.map (fst . headTail) nus
>       nu2s = V.map (fst . headTail . snd . headTail) nus
>       nu3s = V.map (fst . headTail . snd . headTail . snd . headTail) nus
>       nu4s = V.map (fst . headTail . snd . headTail . snd . headTail .
>                     snd . headTail) nus
>       nu5s = V.map (fst . headTail . snd . headTail . snd . headTail .
>                     snd . headTail . snd . headTail) nus
>
>   return (V.zipWith5 SystemState nu1s nu2s nu3s nu4s nu5s)

> newtype SystemObs a = SystemObs { obsI  :: a }
>   deriving Show

> obsUpdate :: Particles (SystemState Double) ->
>              Particles (SystemObs Double)
> obsUpdate xs = V.map (SystemObs . stateI) xs

> priorMu :: R 5
> priorMu = vector [log initS, log initI, log initR, log {- 0.005 -} delta, log {- 0.4 -} gamma]

> bigP :: Sym 5
> bigP = sym $ matrix [
>     1e-6, 0.0, 0.0, 0.0, 0.0
>   , 0.0, 1e-6, 0.0, 0.0, 0.0
>   , 0.0, 0.0, 1e-6, 0.0, 0.0
>   , 0.0, 0.0, 0.0, 5e-3, 0.0
>   , 0.0, 0.0, 0.0, 0.0, 5e-3
>   ]

> initParticles :: MonadRandom m =>
>                  m (Particles (SystemState Double))
> initParticles = V.replicateM nParticles $ do
>   r <- sample $ rvar (MultivariateNormal priorMu bigP)
>   let x1 = exp $ fst $ headTail r
>       x2 = exp $ fst $ headTail $ snd $ headTail r
>       x3 = exp $ fst $ headTail $ snd $ headTail $ snd $ headTail r
>       x4 = exp $ fst $ headTail $ snd $ headTail $ snd $ headTail $
>            snd $ headTail r
>       x5 = exp $ fst $ headTail $ snd $ headTail $ snd $ headTail $
>            snd $ headTail $ snd $ headTail r
>   return $ SystemState { stateS = x1, stateI = x2, stateR = x3,
>                          stateDelta = x4, stateGamma = x5}

> gens :: [[Double]]
> gens = map toList $ toRows $ tr sol

> obs :: [Double]
> -- obs = [3, 8, 28, 75, 221, 291, 255, 235, 190, 125, 70, 28, 12, 5]
> obs  = gens!!1

> bigR :: Sym 1
> bigR  = sym $ matrix [2.0]

> bigQ :: Sym 5
> bigQ = sym $ matrix
>        [ 1e-4, 0.0, 0.0, 0.0, 0.0
>        , 0.0, 1e-4, 0.0, 0.0, 0.0
>        , 0.0, 0.0, 1e-4, 0.0, 0.0
>        , 0.0, 0.0, 0.0, 1e-3, 0.0
>        , 0.0, 0.0, 0.0, 0.0, 1e-2
>        ]

> weight :: forall a n . KnownNat n =>
>           (a -> R n) ->
>           Sym n ->
>           a -> a -> Double
> weight f bigR obs obsNew = pdf (MultivariateNormal (f obsNew) bigR) (f obs)

> runFilter :: [Particles (SystemState Double)]
> runFilter = snd $ evalState action (pureMT 19)
>   where
>     action = runWriterT $ do
>       xs <- lift $ initParticles
>       V.foldM
>         (oneFilteringStep (stateUpdateNoisy bigQ) obsUpdate (weight f bigR))
>         xs
>         (V.fromList $ map SystemObs obs)

> testFiltering :: [Double]
> testFiltering = map ((/ (fromIntegral nParticles)). sum . V.map stateGamma) runFilter

> testFilteringF :: (([Double], [Double]), [Double])
> testFilteringF = ((s, i), r)
>   where
>     ps = runFilter
>     s = map ((/ (fromIntegral nParticles)). sum . V.map stateS) ps
>     i = map ((/ (fromIntegral nParticles)). sum . V.map stateI) ps
>     r = map ((/ (fromIntegral nParticles)). sum . V.map stateR) ps

> type ParticleStream = [[Double]]

> testFilteringS ::
>   (ParticleStream, (ParticleStream, (ParticleStream, (ParticleStream, ParticleStream))))
> testFilteringS = (s, (i, (r, (d, g))))
>   where
>     ps = runFilter
>     s = map (take 200 . V.toList . V.map stateS) ps
>     i = map (take 200 . V.toList . V.map stateI) ps
>     r = map (take 200 . V.toList . V.map stateR) ps
>     d = map (take 200 . V.toList . V.map stateDelta) ps
>     g = map (take 200 . V.toList . V.map stateGamma) ps

Notes
=====

Helpers for Converting Types
----------------------------

> f :: SystemObs Double -> R 1
> f = vector . pure . obsI


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