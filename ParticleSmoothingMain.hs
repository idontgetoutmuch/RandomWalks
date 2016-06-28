{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module Main (
    main
    ) where

import Diagrams.Prelude
import Diagrams.Backend.CmdLine
import Diagrams.Backend.Cairo.CmdLine

import ParticleSmoothingII
import ParticleSmoothingChart
import qualified Numeric.LinearAlgebra.HMatrix as H

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

mean :: (Fractional a, Foldable t) => t a -> a
mean xs = sum xs / (fromIntegral $ length xs)

main :: IO ()
main = do
  -- let x1s = map ((H.! 0) . fst) $ take 100 carSamples
  -- let x2s = map ((H.! 1) . fst) $ take 100 carSamples
  -- let y1s = map ((H.! 0) . snd) $ take 100 carSamples
  -- let y2s = map ((H.! 1) . snd) $ take 100 carSamples
  -- displayHeader "diagrams/CarPosition.png"
  --               (diag "Car Path and Measurements"
  --                (zip x1s x2s)
  --                (zip y1s y2s))
  -- uss <- testCar
  -- displayHeader "diagrams/CarPositionInf.png"
  --               (diagFitted "Car Path and Measurements"
  --                (zip x1s x2s)
  --                (zip y1s y2s)
  --                (zip (fst uss) (snd uss)))
  -- (xss, yss, wss) <- testSmoothCar
  -- let xs, ys :: [Double]
  --     xs = zipWith (\xs ys -> sum $ zipWith (*) xs ys) xss wss
  --     ys = zipWith (\xs ys -> sum $ zipWith (*) xs ys) yss wss
  -- displayHeader "diagrams/CarPositionSmoothInf.png"
  --               (diagFitted "Smooth Car Path and Measurements"
  --                (zip x1s x2s)
  --                (zip y1s y2s)
  --                (zip xs ys))

  let xs = map fst $ take nSimples simpleSamples
      ys = map snd $ take nSimples simpleSamples
  -- uss <- testSimple1
  -- let zs = map mean uss
  -- displayHeader "diagrams/SimpleFilter.png"
  --               (diagParticles "Simple Filter"
  --                (zip [0.0,1.0..fromIntegral nSimples - 1] xs)
  --                (zip [0.0,1.0..fromIntegral nSimples - 1] ys)
  --                [zip [0.0,1.0..fromIntegral nSimples - 1] zs])
  (pss, vss) <- testSimple2
  let zs = zipWith (\xs ys -> sum $ zipWith (*) xs ys) pss vss
  displayHeader "diagrams/SimpleFilter.png"
                (diagParticles "Simple Filter"
                 (zip [0.0,1.0..fromIntegral nSimples - 1] xs)
                 (zip [0.0,1.0..fromIntegral nSimples - 1] ys)
                 [zip [0.0,1.0..fromIntegral nSimples - 1] zs])
  (xss, wss) <- testSmoothSimple
  let l = length xss
  let xs = zipWith (\xs ys -> sum $ zipWith (*) xs ys) xss (reverse wss)
  displayHeader "diagrams/SimpleSmoothInf.png"
                (diagSmoothed "Path and Measurements"
                 (zip [0.0,1.0..] (take l $ map fst simpleSamples))
                 (zip [0.0,1.0..] (take l $ map snd simpleSamples))
                 (zip [0.0,1.0..] zs)
                 (zip [0.0,1.0..] xs))
  -- xsss <- testSimple
  -- let l = length $ head $ last xsss
  -- displayHeader "diagrams/Smooth.png"
  --               (diagSmooth "Smoothing Degeneracy"
  --                (zip [0.0,1.0..] (take l $ map fst simpleSamples))
  --                (concatMap (\xss -> (map (\xs -> zip [0.0,1.0..] xs) xss)) xsss))
  putStrLn "hello"
