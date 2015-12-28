{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable )
import Diagrams.Backend.CmdLine

import System.IO.Unsafe

import PendulumSamples
import Numeric.LinearAlgebra.Static



denv :: DEnv Double
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

diagFitted :: String ->
              Double ->
             [(Double, Double)] ->
             [(Double, Double)] ->
             Diagram Cairo
diagFitted t r l xs =
  fst $ runBackend denv (render (chartFitted t r l xs) (600, 500))

chartFitted :: String ->
               Double ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              Renderable ()
chartFitted title range acts obs = toRenderable layout
  where

    actuals = plot_lines_values .~ [acts]
            $ plot_lines_style  . line_color .~ opaque red
            $ plot_lines_title .~ "Actual Trajectory"
            $ plot_lines_style  . line_width .~ 1.0
            $ def

    measurements = plot_points_values .~ obs
                 $ plot_points_style  . point_color .~ opaque blue
                 $ plot_points_title .~ "Measurements"
                 $ def

    layout = layoutlr_title .~ title
           $ layoutlr_left_axis . laxis_override .~ axisGridHide
           $ layoutlr_left_axis . laxis_title .~ "Actual Angle"
           $ layoutlr_left_axis . laxis_generate .~ scaledAxis def (-range, range)
           $ layoutlr_right_axis . laxis_override .~ axisGridHide
           $ layoutlr_right_axis . laxis_title .~ "Measured Horizontal Displacement"
           $ layoutlr_right_axis . laxis_generate .~ scaledAxis def (-range, range)
           $ layoutlr_x_axis . laxis_override .~ axisGridHide
           $ layoutlr_plots .~ [Left (toPlot actuals),
                                Right (toPlot measurements)]
           $ layoutlr_grid_last .~ False
           $ def


displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  let xs = take 500 $ pendulumSamples
      states = map fst $ map headTail $ map fst xs
      obs    = map fst $ map headTail $ map snd xs
  displayHeader "diagrams/PendulumObs.png"
                (diagFitted "Pendulum" 0.04 (zip [0,1..] states) (zip [0,1..] obs))
  let xs = take 500 $ pendulumSamples'
      states = map fst $ map headTail $ map fst xs
      obs    = map fst $ map headTail $ map snd xs
  displayHeader "diagrams/PendulumObs1.png"
                (diagFitted "Pendulum" 2.5 (zip [0,1..] states) (zip [0,1..] obs))
  putStrLn "Hello"