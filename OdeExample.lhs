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
> {-# LANGUAGE DataKinds                    #-}
> {-# LANGUAGE ExplicitForAll               #-}

> module OdeExample where

> import Numeric.GSL.ODE
> import Numeric.LinearAlgebra hiding ( R, vector, matrix, sym )

> import           Data.Random hiding ( StdNormal, Normal, gamma )
> import           Data.Random.Source.PureMT ( pureMT )
> import           Control.Monad.State ( evalState, replicateM )
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


> sirOde :: Double -> Double -> Double -> [Double] -> [Double]
> sirOde beta gamma _t [s, i, _r] =
>   [
>     negate (beta * i * s)
>   , (beta * i * s) - (gamma * i)
>   , gamma * i
>   ]
> sirOde _b _g _t vars = error $ "sirOde called with: " ++ show (length vars) ++ " variable"

> beta, gamma :: Double
> beta = 0.0026
> gamma = 0.5

> initS, initI, initR :: Double
> initS = 762.0
> initI = 1.0
> initR = 0.0

> sol :: Matrix Double
> sol = odeSolve (sirOde beta gamma) [initS, initI, initR] (fromList [0.0,deltaT..14.0])

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

> data SystemState a = SystemState { stateS     :: a
>                                  , stateI     :: a
>                                  , stateR     :: a
>                                  , stateBeta  :: a
>                                  , stateGamma :: a
>                                  }
>   deriving Show

> deltaT :: Double
> deltaT = 0.2

> stateUpdate :: Particles (SystemState Double) ->
>                Particles (SystemState Double)
> stateUpdate xPrevs =
>   V.zipWith5 SystemState ss is rs betas gammas
>   where
>     sPrevs     = V.map stateS     xPrevs
>     iPrevs     = V.map stateI     xPrevs
>     rPrevs     = V.map stateR     xPrevs
>     betaPrevs  = V.map stateBeta  xPrevs
>     gammaPrevs = V.map stateGamma xPrevs
>
>     f b g xs = odeSolve (sirOde b g) xs
>                (0.0 `VS.cons` (deltaT `VS.cons` VS.empty))
>
>     ms :: V.Vector (Matrix Double)
>     ms = V.zipWith3 f betaPrevs gammaPrevs
>                       (V.zipWith3 (\s i r -> [s, i, r]) sPrevs iPrevs rPrevs)
>
>     ns :: V.Vector (VS.Vector Double)
>     ns = V.map (\m -> (toRows m)!!1) ms
>
>     ss     = V.map (VS.! 0) ns
>     is     = V.map (VS.! 1) ns
>     rs     = V.map (VS.! 2) ns
>     betas  = betaPrevs
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
>       sPrevs     = V.map stateS     xs
>       iPrevs     = V.map stateI     xs
>       rPrevs     = V.map stateR     xs
>       betaPrevs  = V.map stateBeta  xs
>       gammaPrevs = V.map stateGamma xs
>
>   let ix     = V.length xPrevs
>   etas <- replicateM ix $ sample $ rvar (MultivariateNormal 0.0 bigQ)
>   let eta1s, eta2s, eta3s, eta4s, eta5s :: V.Vector Double
>       eta1s = V.fromList $ map (fst . headTail) etas
>       eta2s = V.fromList $ map (fst . headTail . snd . headTail) etas
>       eta3s = V.fromList $ map (fst . headTail . snd . headTail . snd . headTail) etas
>       eta4s = V.fromList $ map (fst . headTail . snd . headTail . snd . headTail .
>                                 snd . headTail) etas
>       eta5s = V.fromList $ map (fst . headTail . snd . headTail . snd . headTail .
>                                 snd . headTail . snd . headTail) etas
>
>   return (V.zipWith5 SystemState (sPrevs .+ eta1s) (iPrevs .+ eta2s) (rPrevs .+ eta3s)
>                                  (betaPrevs .+ eta4s) (gammaPrevs .+ eta5s))

> newtype SystemObs a = SystemObs { obsI  :: a }
>   deriving Show

> obsUpdate :: Particles (SystemState Double) ->
>              Particles (SystemObs Double)
> obsUpdate xs = V.map (SystemObs . stateI) xs

> priorMu :: R 5
> priorMu = vector [initS, initI, initR, log beta, log gamma]

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
>   let x1 = fst $ headTail r
>       x2 = fst $ headTail $ snd $ headTail r
>       x3 = fst $ headTail $ snd $ headTail $ snd $ headTail r
>       x4 = exp $ fst $ headTail $ snd $ headTail $ snd $ headTail $
>            snd $ headTail r
>       x5 = exp $ fst $ headTail $ snd $ headTail $ snd $ headTail $
>            snd $ headTail $ snd $ headTail r
>   return $ SystemState { stateS = x1, stateI = x2, stateR = x3,
>                          stateBeta = x4, stateGamma = x5}

> gens :: [[Double]]
> gens = map toList $ toRows $ tr sol

> obs :: [Double]
> -- obs = [3, 8, 28, 75, 221, 291, 255, 235, 190, 125, 70, 28, 12, 5]
> obs  = gens!!1

> bigR :: Sym 1
> bigR  = sym $ matrix [0.05]

> bigQ :: Sym 5
> bigQ = sym $ matrix
>        [ 1e-6, 0.0, 0.0, 0.0, 0.0
>        , 0.0, 1e-6, 0.0, 0.0, 0.0
>        , 0.0, 0.0, 1e-6, 0.0, 0.0
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

> sol1 :: Matrix Double
> sol1 = odeSolve (sirOde 0.005 2.5) [initS, initI, initR] (fromList [0.0,deltaT..14.0])

> testFilteringF :: (([Double], [Double]), [Double])
> testFilteringF = ((s, i), r)
>   where
>     ps = runFilter
>     s = map ((/ (fromIntegral nParticles)). sum . V.map stateS) ps
>     i = map ((/ (fromIntegral nParticles)). sum . V.map stateI) ps
>     r = map ((/ (fromIntegral nParticles)). sum . V.map stateR) ps

> testFilteringS :: (([[Double]], [[Double]]), [[Double]])
> testFilteringS = ((s, i), r)
>   where
>     ps = runFilter
>     s = map (take 200 . V.toList . V.map stateS) ps
>     i = map (take 200 . V.toList . V.map stateI) ps
>     r = map (take 200 . V.toList . V.map stateR) ps

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

Bibliography
============