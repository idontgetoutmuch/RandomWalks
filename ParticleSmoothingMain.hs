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

import ParticleSmoothing
import ParticleSmoothingChart
import qualified Numeric.LinearAlgebra.HMatrix as H

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  let x1s = map ((H.! 0) . fst) $ take 100 carSamples
  let x2s = map ((H.! 1) . fst) $ take 100 carSamples
  let y1s = map ((H.! 0) . snd) $ take 100 carSamples
  let y2s = map ((H.! 1) . snd) $ take 100 carSamples
  displayHeader "diagrams/CarPosition.png"
                (diag "Car Path and Measurements"
                 (zip x1s x2s)
                 (zip y1s y2s))
  xss <- test1
  let l = length $ head xss
  displayHeader "diagrams/Smooth.png"
                (diagSmooth "Smoothing Degeneracy"
                 (zip [0.0,1.0..] (take l $ map fst carSamples1))
                 (map (\xs -> zip [0.0,1.0..] xs) xss))
  putStrLn "hello"
