{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module ParticleSmoothingChart (
    diag
  , diagSmooth
  , diagFitted
  , diagSmoothed
  , diagParticles
  ) where

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable )

import System.IO.Unsafe


denv :: DEnv Double
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

diag :: String ->
        [(Double, Double)] ->
        [(Double, Double)] ->
        Diagram Cairo
diag t l xs =
  fst $ runBackend denv (render (chart t l xs) (600, 500))

chart :: String ->
         [(Double, Double)] ->
         [(Double, Double)] ->
         Renderable ()
chart t l obs = toRenderable layout
  where

    boundry = plot_lines_values .~ [l]
              $ plot_lines_style  . line_color .~ opaque red
              $ plot_lines_title .~ "Actual Trajectory"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    actuals = plot_points_values .~ obs
              $ plot_points_style  . point_color .~ opaque blue
              $ plot_points_title .~ "Measurements"
              $ def

    layout = layout_title .~ t
           $ layout_plots .~ [toPlot actuals, toPlot boundry]
           $ layout_y_axis . laxis_title .~ "y co-ordinate"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "x co-ordindate"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

diagFitted :: String ->
             [(Double, Double)] ->
             [(Double, Double)] ->
             [(Double, Double)] ->
             Diagram Cairo
diagFitted t l xs es =
  fst $ runBackend denv (render (chartFitted t l xs es) (600, 500))

chartFitted :: String ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              Renderable ()
chartFitted t l obs ests = toRenderable layout
  where

    boundry = plot_lines_values .~ [l]
              $ plot_lines_style  . line_color .~ opaque red
              $ plot_lines_title .~ "Actual Trajectory"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    estimas = plot_lines_values .~ [ests]
              $ plot_lines_style  . line_color .~ opaque green
              $ plot_lines_title .~ "Inferred Trajectory"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    actuals = plot_points_values .~ obs
              $ plot_points_style  . point_color .~ opaque blue
              $ plot_points_title .~ "Measurements"
              $ def

    layout = layout_title .~ t
           $ layout_plots .~ [toPlot actuals, toPlot boundry, toPlot estimas]
           $ layout_y_axis . laxis_title .~ "y co-ordinate"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "x co-ordindate"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

diagSmoothed :: String ->
                [(Double, Double)] ->
                [(Double, Double)] ->
                [(Double, Double)] ->
                [(Double, Double)] ->
                Diagram Cairo
diagSmoothed t l xs es em =
  fst $ runBackend denv (render (chartSmoothed t l xs es em) (600, 500))

chartSmoothed :: String ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              Renderable ()
chartSmoothed t l obs ests sms = toRenderable layout
  where

    boundry = plot_lines_values .~ [l]
              $ plot_lines_style  . line_color .~ opaque red
              $ plot_lines_title .~ "Actual Trajectory"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    estimas = plot_lines_values .~ [ests]
              $ plot_lines_style  . line_color .~ opaque green
              $ plot_lines_title .~ "Inferred Trajectory"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    estimsm = plot_lines_values .~ [sms]
              $ plot_lines_style  . line_color .~ opaque yellow
              $ plot_lines_title .~ "Smoothed Trajectory"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    actuals = plot_points_values .~ obs
              $ plot_points_style  . point_color .~ opaque blue
              $ plot_points_title .~ "Measurements"
              $ def

    layout = layout_title .~ t
           $ layout_plots .~ [toPlot actuals, toPlot boundry, toPlot estimas, toPlot estimsm]
           $ layout_y_axis . laxis_title .~ "y co-ordinate"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "x co-ordindate"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

chartSmooth :: String ->
              [(Double, Double)] ->
              [[(Double, Double)]] ->
              Renderable ()
chartSmooth t l obss = toRenderable layout
  where

    boundry = plot_points_values .~ l
              $ plot_points_title .~ "Actual"
              $ plot_points_style .~ filledCircles 4 (red `withOpacity` 0.25)
              $ def

    actuals obs = plot_points_values .~ obs
                  $ plot_points_style  .~  filledCircles 2 (blue `withOpacity` 0.25)
                  $ def

    paths = plot_lines_values .~ obss
            $ plot_lines_style  . line_color .~ (green `withOpacity` 0.1)
            $ plot_lines_title .~ "Particle Trajectories"
            -- $ plot_lines_style  . line_width .~ 1.0
            $ def

    layout = layout_title .~ t
           $ layout_plots .~ (toPlot paths) : (toPlot boundry) : (map (toPlot . actuals) obss)
           $ layout_y_axis . laxis_title .~ "Position"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

diagSmooth :: String ->
             [(Double, Double)] ->
             [[(Double, Double)]] ->
             Diagram Cairo
diagSmooth t l xss =
  fst $ runBackend denv (render (chartSmooth t l xss) (600, 500))

chartParticles :: String ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              [[(Double, Double)]] ->
              Renderable ()
chartParticles t xs ys zss = toRenderable layout
  where

    actual  = plot_lines_values .~ [xs]
              $ plot_lines_style  . line_color .~ (blue `withOpacity` 0.7)
              $ plot_lines_title .~ "Actual"
              $ def

    observd = plot_points_values .~ ys
              $ plot_points_title .~ "Observed"
              $ plot_points_style .~ filledCircles 4 (red `withOpacity` 0.25)
              $ def

    paths = plot_lines_values .~ zss
            $ plot_lines_style  . line_color .~ (green `withOpacity` 0.7)
            $ plot_lines_title .~ "Particle Mean"
            $ def

    layout = layout_title .~ t
           $ layout_plots .~ (toPlot actual) : (toPlot observd) : (toPlot paths) : []
           $ layout_y_axis . laxis_title .~ "Position"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

diagParticles :: String ->
             [(Double, Double)] ->
             [(Double, Double)] ->
             [[(Double, Double)]] ->
             Diagram Cairo
diagParticles t xs ys zss =
  fst $ runBackend denv (render (chartParticles t xs ys zss) (600, 500))
