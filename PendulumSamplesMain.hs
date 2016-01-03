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

diagEstimated :: String ->
                 [(Double, Double)] ->
                 [(Double, Double)] ->
                 [(Double, Double)] ->
                 Diagram Cairo
diagEstimated t l xs es =
  fst $ runBackend denv (render (chartEstimated t l xs es) (600, 500))

chartEstimated :: String ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              Renderable ()
chartEstimated title acts obs ests = toRenderable layout
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

    estimas = plot_lines_values .~ [ests]
            $ plot_lines_style  . line_color .~ opaque black
            $ plot_lines_title .~ "Inferred Trajectory"
            $ plot_lines_style  . line_width .~ 1.0
            $ def

    layout = layout_title .~ title
           $ layout_plots .~ [toPlot actuals, toPlot measurements, toPlot estimas]
           $ layout_y_axis . laxis_title .~ "Angle / Horizontal Displacement"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def


displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  -- let xs = take 500 $ pendulumSamples
  --     states = map fst $ map headTail $ map fst xs
  --     obs    = map fst $ map headTail $ map snd xs
  -- displayHeader "diagrams/PendulumObs.png"
  --               (diagFitted "Pendulum" 0.04 (zip [0,1..] states) (zip [0,1..] obs))
  -- let xs = take 500 $ pendulumSamples'
  --     states = map fst $ map headTail $ map fst xs
  --     obs    = map fst $ map headTail $ map snd xs
  -- displayHeader "diagrams/PendulumObs1.png"
  --               (diagFitted "Pendulum" 2.5 (zip [0,1..] states) (zip [0,1..] obs))
  -- let x1s = testFiltering 500
  -- displayHeader "diagrams/PendulumFitted.png"
  --               (diagEstimated "Fitted Pendulum"
  --                              (zip [0,1..] states)
  --                              (zip [0,1..] obs)
  --                              (zip [0,1..] x1s))
  -- let x1ss = reverse $ testSmoothing 500 20
  -- displayHeader "diagrams/PendulumSmoothed20.png"
  --               (diagEstimated "Smoothed Pendulum"
  --                              (zip [0,1..] states)
  --                              (zip [0,1..] obs)
  --                              (zip [0,1..] x1ss))
  let x3s = testFilteringG 1000
  displayHeader "diagrams/PendulumG.png"
                (diagFitted "Gravity" 12.0 (zip [0,1..] (replicate 1000 9.81)) (zip [0,1..] x3s))
  putStrLn "Hello"