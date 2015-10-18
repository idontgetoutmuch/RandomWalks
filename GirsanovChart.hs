{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module GirsanovChart (
    diag
  ) where

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable )

import System.IO.Unsafe


denv :: DEnv Double
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

diag :: String ->
        [[(Double, Double)]] ->
        Diagram Cairo
diag t xss =
  fst $ runBackend denv (render (chart t xss) (600, 500))

chart :: String ->
         [[(Double, Double)]] ->
         Renderable ()
chart t obss = toRenderable layout
  where

    actuals = plot_lines_values .~ obss
              $ plot_lines_style  . line_color .~ opaque blue
              $ plot_lines_title .~ "Path"
              $ plot_lines_style  . line_width .~ 1.0
              $ def

    layout = layout_title .~ t
           $ layout_plots .~ [toPlot actuals]
           $ layout_y_axis . laxis_title .~ "Value"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def
