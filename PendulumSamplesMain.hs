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

import qualified Diagrams.Backend.Rasterific as R


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


diagFittedG :: String ->
               Double ->
               [(Double, Double)] ->
               [(Double, Double)] ->
               Diagram Cairo
diagFittedG t r l xs =
  fst $ runBackend denv (render (chartFittedG t r l xs) (600, 500))

chartFittedG :: String ->
                Double ->
               [(Double, Double)] ->
               [(Double, Double)] ->
               Renderable ()
chartFittedG title range acts obs = toRenderable layout
  where

    actuals = plot_lines_values .~ [acts]
            $ plot_lines_style  . line_color .~ opaque red
            $ plot_lines_title .~ "Actual Gravity"
            $ plot_lines_style  . line_width .~ 1.0
            $ def

    measurements = plot_lines_values .~ [obs]
                 $ plot_lines_style  . line_color .~ opaque blue
                 $ plot_lines_title .~ "Inferred Gravity"
                 $ plot_lines_style  . line_width .~ 1.0
                 $ def

    layout = layout_title .~ title
           $ layout_plots .~ [toPlot actuals, toPlot measurements]
           $ layout_y_axis . laxis_title .~ "Acceleration"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_y_axis . laxis_generate .~ scaledAxis def (0.0, range)
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ layout_grid_last .~ False
           $ def

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

-- dia :: Diagram B
-- dia :: QDiagram Cairo V2 Double Any
dia :: QDiagram R.Rasterific V2 Double Any
dia = image (DImage (ImageRef "diagrams/PendulumObs.B.png")
                    600 600 (translationX (0.0 :: Double)))

main :: IO ()
main = do
  let xs = take 500 $ pendulumSamples
      states = map fst $ map headTail $ map fst xs
      obs    = map fst $ map headTail $ map snd xs
  displayHeader "diagrams/PendulumObs.B.png"
                (diagFitted "Pendulum" 0.04 (zip [0,1..] states) (zip [0,1..] obs))

  let xs = take 500 $ pendulumSamples'
      states = map fst $ map headTail $ map fst xs
      obs    = map fst $ map headTail $ map snd xs
  displayHeader "diagrams/PendulumObs1.B.png"
                (diagFitted "Pendulum" 2.5 (zip [0,1..] states) (zip [0,1..] obs))


  let x1s = testFiltering 500
  putStrLn $ show $ last x1s
  putStrLn $ show $ sum $ zipWith (\x y -> (x - y)^2) states x1s
  displayHeader "diagrams/PendulumFitted.B.png"
                (diagEstimated "Fitted Pendulum"
                               (zip [0,1..] states)
                               (zip [0,1..] obs)
                               (zip [0,1..] x1s))

  let x1ss = reverse $ testSmoothing 500 20
  putStrLn $ show $ last x1ss
  putStrLn $ show $ sum $ zipWith (\x y -> (x - y)^2) states x1ss
  displayHeader "diagrams/PendulumSmoothed20.B.png"
                (diagEstimated "Smoothed Pendulum"
                               (zip [0,1..] states)
                               (zip [0,1..] obs)
                               (zip [0,1..] x1ss))

  let x3s = testFilteringG 500
  putStrLn $ show $ last x3s
  displayHeader "diagrams/PendulumG.B.png"
                (diagFittedG "Gravity" 12.0 (zip [0,1..] (replicate 500 9.81)) (zip [0,1..] x3s))

  let us = testSmoothingG 500 20
      x1s = reverse $ (\(v, _, _) -> v) us
      x3s = reverse $ (\(_, _, v) -> v) us
  putStrLn $ show $ last x3s
  displayHeader "diagrams/PendulumSmoothedG.B.png"
                (diagFittedG "Smoothed Gravity" 12.0 (zip [0,1..] (replicate 500 9.81)) (zip [0,1..] x3s))

  displayHeader "diagrams/PendulumSmoothedG1X20.B.png"
                (diagEstimated "Smoothed Pendulum"
                               (zip [0,1..] states)
                               (zip [0,1..] obs)
                               (zip [0,1..] x1s))
  putStrLn "Hello"