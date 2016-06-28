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

import OdeExample
import Numeric.LinearAlgebra


denv :: DEnv Double
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

diagSirGen :: String ->
                 [(Double, Double)] ->
                 [(Double, Double)] ->
                 [(Double, Double)] ->
                 Diagram Cairo
diagSirGen t l xs es =
  fst $ runBackend denv (render (chartSirGen t l xs es) (600, 500))

chartSirGen :: String ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              [(Double, Double)] ->
              Renderable ()
chartSirGen title acts obs ests = toRenderable layout
  where

    actuals = plot_lines_values .~ [acts]
            $ plot_lines_style  . line_color .~ opaque red
            $ plot_lines_title .~ "Susceptible"
            $ plot_lines_style  . line_width .~ 2.0
            $ def

    measurements = plot_lines_values .~ [obs]
                 $ plot_lines_style  . line_color .~ opaque blue
            $ plot_lines_style  . line_width .~ 2.0
                 $ plot_lines_title .~ "Infected"
                 $ def

    estimas = plot_lines_values .~ [ests]
            $ plot_lines_style  . line_color .~ opaque black
            $ plot_lines_title .~ "Recovered"
            $ plot_lines_style  . line_width .~ 2.0
            $ def

    layout = layout_title .~ title
           $ layout_plots .~ [toPlot actuals, toPlot measurements, toPlot estimas]
           $ layout_y_axis . laxis_title .~ "Pupil Numbers"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

diagSirParticles :: String ->
                    String ->
                    [(Double, Double)] ->
                    Diagram Cairo
diagSirParticles t y l =
  fst $ runBackend denv (render (chartSirParticles t y l) (600, 500))

chartSirParticles :: String ->
                     String ->
                     [(Double, Double)] ->
                     Renderable ()
chartSirParticles title ytitle acts = toRenderable layout
  where

    actuals = plot_points_values .~ acts
            $ plot_points_style  . point_color .~ opaque red
            $ plot_points_title .~ "Particles"
            $ def

    layout = layout_title .~ title
           $ layout_plots .~ [toPlot actuals]
           $ layout_y_axis . laxis_title .~ ytitle
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
  let xs = map toList $ toRows $ tr sol
      ss = xs!!0
      is = xs!!1
      rs = xs!!2
  displayHeader "diagrams/Sir.png"
                (diagSirGen "Influenza Outbreak"
                               (zip [0,1..] ss)
                               (zip [0,1..] is)
                               (zip [0,1..] rs))
  let xs = testFilteringF
      ss = fst $ fst xs
      is = snd $ fst xs
      rs = snd xs
  displayHeader "diagrams/Sir1.png"
                (diagSirGen "Influenza Outbreak"
                               (zip [0,1..] ss)
                               (zip [0,1..] is)
                               (zip [0,1..] rs))
  let ys = testFilteringS
      is = fst $ snd ys
      js = concat $ zipWith (\i ps -> zip (repeat i) ps) [0,1..] is
  displayHeader "diagrams/SirParts.png"
                (diagSirParticles "Influenza Outbreak" "Infected" js)

  let ds = fst $ snd $ snd $ snd ys
      js = concat $ zipWith (\i ps -> zip (repeat i) ps) [0,1..] ds
  displayHeader "diagrams/SirDeltaParts.png"
                (diagSirParticles "Influenza Outbreak" "Infection Parameter" js)

  let gs = snd $ snd $ snd $ snd ys
      js = concat $ zipWith (\i ps -> zip (repeat i) ps) [0,1..] gs
  displayHeader "diagrams/SirGammaParts.png"
                (diagSirParticles "Influenza Outbreak" "Recovery Parameter" js)
  putStrLn "Hello"