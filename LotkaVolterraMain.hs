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
import Numeric.LinearAlgebra


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

main :: IO ()
main = do
  -- let xs = map toList $ toRows $ tr solLv
  --     hs = xs!!0
  --     ls = xs!!1
  -- displayHeader "diagrams/LV.png"
  --               (diagSirGen "Hares and Lynxes"
  --                              (zip [0,1..] hs)
  --                              (zip [0,1..] ls))
  -- let xs = map toList $ toRows $ tr solPz
  --     ps = xs!!0
  --     zs = xs!!1
  -- displayHeader "diagrams/PZ.png"
  --               (diagSirGen "Phyto and Zoo Plankton"
  --                              (zip ps zs)
  --                              (zip ps zs))

  let xs = map toList $ toRows $ tr $ solPp h0 l0 -- 100.0 50.0
      ps = xs!!0
      zs = xs!!1
  let ys = map toList $ toRows $ tr $ solPp h0 l0 -- 100.0 25.0
      ps' = ys!!0
      zs' = ys!!1
  displayHeader "diagrams/PP.png"
                (diagSirGen "Hares and Lynxes"
                               (zip ps zs)
                               (zip ps' zs'))

  putStrLn "Hello"