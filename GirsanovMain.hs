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

import Control.Monad

import Girsanov
import GirsanovChart

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  let sims = map (\seed -> bM0to1 scanl replicateM seed 1000) [0..99]
  displayHeader "diagrams/BrownianPaths.png"
                (diag "Sample Paths for Brownian Motion"
                 (let ls = [0..1000 - 1] in zip ls (map (\x -> x*2e-3 + 2) ls) )
                 (map (zip [0..]) sims))
  putStrLn "hello"
