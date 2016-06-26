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

import LotkaVolterra
import Numeric.LinearAlgebra hiding ( diag )


denv :: DEnv Double
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

diagSirGen :: String ->
                 [(Double, Double)] ->
                 [(Double, Double)] ->
                 Diagram Cairo
diagSirGen t l xs =
  fst $ runBackend denv (render (chartSirGen t l xs) (600, 500))

chartSirGen :: String ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              Renderable ()
chartSirGen title acts obs = toRenderable layout
  where

    actuals = plot_lines_values .~ [acts]
            $ plot_lines_style  . line_color .~ opaque red
            $ plot_lines_title .~ "Hares=100, Lynxes=50"
            $ plot_lines_style  . line_width .~ 2.0
            $ def

    measurements = plot_lines_values .~ [obs]
                 $ plot_lines_style  . line_color .~ opaque blue
            $ plot_lines_style  . line_width .~ 2.0
                 $ plot_lines_title .~ "Hares=100, Lynxes=25"
                 $ def

    layout = layout_title .~ title
           $ layout_plots .~ [toPlot actuals, toPlot measurements]
           $ layout_y_axis . laxis_title .~ "Lynxes"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Hares"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

sims :: [[Double]]
sims = map (logBM 1.0 0.05 100 1) [1..10]

sims' :: [[(Double, Double, Double)]]
sims' = map (eulerEx (log 0.5, 100.0, 50.0) 0.05 100 30) [3]

rho1ss, pss, zss :: [[Double]]
rho1ss = map (map (\(x, _, _) -> x)) sims'
pss    = map (map (\(_, y, _) -> y)) sims'
zss    = map (map (\(_, _, z) -> z)) sims'

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

diag :: String ->
        [[(Double, Double)]] ->
        Diagram Cairo
diag t xss =
  fst $ runBackend denv (render (chart t xss) (600, 500))

main :: IO ()
main = do
  -- let xs = map toList $ toRows $ tr $ solPp h0 l0 -- 100.0 50.0
  --     ps = xs!!0
  --     zs = xs!!1
  -- let ys = map toList $ toRows $ tr $ solPp h0 l0 -- 100.0 25.0
  --     ps' = ys!!0
  --     zs' = ys!!1
  -- displayHeader "diagrams/PP.png"
  --               (diagSirGen "Hares and Lynxes"
  --                              (zip ps zs)
  --                              (zip ps' zs'))

  -- displayHeader "diagrams/LogBrownianPaths.png"
  --               (diag "Sample Paths for Log Brownian Motion"
  --                (map (zip [0..]) sims))

  displayHeader "diagrams/StochPathsA.png"
                (diag "Sample Paths for Stochastic LV"
                 [zip [0..] (pss!!0), zip [0..] (zss!!0)])

  putStrLn "Hello"