% Modelling an Epidemic
% Dominic Steinitz
% 28th March 2016

---
bibliography: Stochastic.bib
---

Introduction
============

This is a bit different from my usual posts (well apart from my write
up of hacking at Odessa) in that it is a log of how I managed to get
[LibBi] (http://libbi.org)(Library for Bayesian Inference) to run on my
MacBook and then not totally satisfactorily (as you will see if you
read on).

The intention is to try a few more approaches to the same problem, for
example, [Stan](http://mc-stan.org),
[monad-bayes](https://github.com/adscib/monad-bayes) and hand-crafted.

@RSPSA1927:115 give a simple model of the spread of an infectious
disease. Individuals move from being susceptible ($S$) to infected
($I$) to recovered ($R$).

$$
\begin{eqnarray}
\frac{dS}{dt} & = & - \delta S(t) I(t) \label{eq2a} \\
\frac{dI}{dt} & = & \delta S(t) I(t) - \gamma I(t) \label{eq2b} \\
\frac{dR}{dt} & = & \gamma I(t) . \label{eq2c}
\end{eqnarray}
$$

In 1978, anonymous authors sent a note to the British Medical Journal
reporting an influenza outbreak in a boarding school in the north of
England (@bmj-influenza). The chart below shows the solution of the
SIR (Susceptible, Infected, Record) model with parameters which give
roughly the results observed in the school.

![](diagrams/Sir.png)

LibBi
=====

Step 1
------

~~~
~/LibBi-stable/SIR-master $ ./init.sh
error: 'ncread' undefined near line 6 column 7
~~~

The README says this is optional so we can skip over it. Still it
would be nice to fit the bridge weight function as described in
@Moral2015.

The README does say that
[GPML](http://www.gaussianprocess.org/gpml/code/matlab/doc/) is
required but since we don't (yet) need to do this step, let's move on.

~~~
~/LibBi-stable/SIR-master $ ./run.sh
./run.sh

Error: ./configure failed with return code 77. See
.SIR/build_openmp_cuda_single/configure.log and
.SIR/build_openmp_cuda_single/config.log for details
~~~

It seems the example is configured to run on CUDA and it is highly
likely that my installation of LibBI was not set up to allow this. We
can change `config.conf` from

~~~
--disable-assert
--enable-single
--enable-cuda
--nthreads 2
~~~

to

~~~
--nthreads 4
--enable-sse
--disable-assert
~~~

On to the next issue.

~~~
~/LibBi-stable/SIR-master $ ./run.sh
./run.sh
Error: ./configure failed with return code 1. required QRUpdate
library not found. See .SIR/build_sse/configure.log and
.SIR/build_sse/config.log for details
~~~

But QRUpdate is installed!

~~~
~/LibBi-stable/SIR-master $ brew info QRUpdate
brew info QRUpdate
homebrew/science/qrupdate: stable 1.1.2 (bottled)
http://sourceforge.net/projects/qrupdate/
/usr/local/Cellar/qrupdate/1.1.2 (3 files, 302.6K)
/usr/local/Cellar/qrupdate/1.1.2_2 (6 files, 336.3K)
  Poured from bottle
/usr/local/Cellar/qrupdate/1.1.2_3 (6 files, 337.3K) *
  Poured from bottle
From: https://github.com/Homebrew/homebrew-science/blob/master/qrupdate.rb
==> Dependencies
Required: veclibfort ✔
Optional: openblas ✔
==> Options
--with-openblas
	Build with openblas support
--without-check
	Skip build-time tests (not recommended)
~~~

Let's look in the log as advised. So it seems that a certain symbol
cannot be found.

~~~
checking for dch1dn_ in -lqrupdate
~~~

Let's try ourselves.

~~~
nm -g /usr/local/Cellar/qrupdate/1.1.2_3/lib/libqrupdate.a | grep dch1dn_
0000000000000000 T _dch1dn_
~~~

So the symbol is there! What gives? Let's try setting one of the
environment variables.

~~~
export LDFLAGS='-L/usr/local/lib'
~~~

Now we get further.

~~~
./run.sh
Error: ./configure failed with return code 1. required NetCDF header
not found. See .SIR/build_sse/configure.log and
.SIR/build_sse/config.log for details
~~~

So we just need to set another environment variable.

~~~
export CPPFLAGS='-I/usr/local/include/'
~~~

This is more mysterious.

~~~
./run.sh
Error: ./configure failed with return code 1. required Boost header
not found. See .SIR/build_sse/configure.log and
.SIR/build_sse/config.log for details ~/LibBi-stable/SIR-master
~~~

Let's see what we have.

~~~
brew list | grep -i boost
~~~

Nothing! I recall having some problems with `boost` when trying to use
a completely different package. So let's install `boost`.

~~~
brew install boost
~~~

Now we get a different error.

~~~
./run.sh
Error: make failed with return code 2, see .SIR/build_sse/make.log for details
~~~

Fortunately at some time in the past [sbfnk](https://github.com/sbfnk)
took pity on me and advised me
[here](https://github.com/libbi/LibBi/issues/5) to use `boost155`, a
step that should not be lightly undertaken.

~~~
/usr/local/Cellar/boost155/1.55.0_1: 10,036 files, 451.6M, built in 15 minutes 9 seconds
~~~

Even then I had to say

~~~
brew link --force boost155
~~~

Finally it runs.

~~~
./run.sh 2> out.txt
~~~

And produces a lot of output

~~~
wc -l out.txt
   49999 out.txt

ls -ltrh results/posterior.nc
   1.7G Apr 30 19:57 results/posterior.nc
~~~

Rather worringly, `out.txt` has all lines of the form

~~~
1: -51.9191 -23.2045 nan beats -inf -inf -inf accept=0.5
~~~

`nan` beating `-inf` does not sound goo.

Now we are in a position to analyse the results.

~~~
octave --path oct/ --eval "plot_and_print"
error: 'bi_plot_quantiles' undefined near line 23 column 5
~~~

I previously found an Octave package(?) called `OctBi` so let's create
an `.octaverc` file which adds this to the path. We'll also need to
load the `netcdf` package which we previously installed.

```
addpath ("../OctBi-stable/inst")
pkg load netcdf
```

~~~
~/LibBi-stable/SIR-master $ octave --path oct/ --eval "plot_and_print"
octave --path oct/ --eval "plot_and_print"
warning: division by zero
warning: called from
    mean at line 117 column 7
    read_hist_simulator at line 47 column 11
    bi_read_hist at line 85 column 12
    bi_hist at line 63 column 12
    plot_and_print at line 56 column 5
warning: division by zero
warning: division by zero
warning: division by zero
warning: division by zero
warning: division by zero
warning: print.m: fig2dev binary is not available.
Some output formats are not available.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
warning: opengl_renderer: x/y/zdata should have the same dimensions. Not rendering.
sh: pdfcrop: command not found
~~~

I actually get a chart from this so some kind of success.

![](diagrams/posterior.png)

This does *not* look like the chart in the @Moral2015, the fitted
number of infected patients looks a lot smoother and the "rates"
parameters also vary in a much smoother manner. For reasons I haven't
yet investigated, it looks like over-fitting.



Haskell Preamble
================

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}

> {-# LANGUAGE BangPatterns                 #-}
> {-# LANGUAGE DataKinds                    #-}
> {-# LANGUAGE ExplicitForAll               #-}

> module OdeExample where

> import Numeric.GSL.ODE
> import Numeric.LinearAlgebra hiding ( R, vector, matrix, sym )

> import           Data.Random hiding ( StdNormal, Normal, gamma )
> import           Data.Random.Source.PureMT ( pureMT )
> import           Control.Monad.State ( evalState )
> import           Control.Monad.Writer ( tell, WriterT, lift,
>                                         runWriterT
>                                       )
> import           Numeric.LinearAlgebra.Static
>                  ( R, vector, Sym,
>                    headTail, matrix, sym
>                  )
> import           GHC.TypeLits ( KnownNat )
> import           MultivariateNormal ( MultivariateNormal(..) )
> import qualified Data.Vector as V
> import qualified Data.Vector.Storable as VS
> import           Data.Bits ( shiftR )
> import           Control.Parallel.Strategies


Solving the Equations
=====================

It is straightforward to solve the equations. The parameters are set to
something like those which give a result something like that observed
in @bmj-influenza and shown in the chart above.

> sirOde :: Double -> Double -> Double -> [Double] -> [Double]
> sirOde delta gamma _t [s, i, _r] =
>   [
>     negate (delta * i * s)
>   , (delta * i * s) - (gamma * i)
>   , gamma * i
>   ]
> sirOde _b _g _t vars = error $ "sirOde called with: " ++ show (length vars) ++ " variable"

> delta, gamma :: Double
> delta = 0.0026
> gamma = 0.5

> initS, initI, initR :: Double
> initS = 762.0
> initI = 1.0
> initR = 0.01

> sol :: Matrix Double
> sol = odeSolve (sirOde delta gamma) [initS, initI, initR] (fromList [0.0,deltaT..14.0])


Bootstrap Particle Filter Naive Haskell Implementation
======================================================

> nParticles :: Int
> nParticles = 10000

The usual Bayesian update step.

> type Particles a = V.Vector a

> oneFilteringStep ::
>   MonadRandom m =>
>   (Particles a -> m (Particles a)) ->
>   (Particles a -> Particles b) ->
>   (b -> b -> Double) ->
>   Particles a ->
>   b ->
>   WriterT [Particles a] m (Particles a)
> oneFilteringStep stateUpdate obsUpdate weight statePrevs obs = do
>   tell [statePrevs]
>   stateNews <- lift $ stateUpdate statePrevs
>   let obsNews = obsUpdate stateNews
>   let weights       = V.map (weight obs) obsNews
>       cumSumWeights = V.tail $ V.scanl (+) 0 weights
>       totWeight     = V.last cumSumWeights
>   vs <- lift $ V.replicateM nParticles (sample $ uniform 0.0 totWeight)
>   let js = indices cumSumWeights vs
>       stateTildes = V.map (stateNews V.!) js
>   return stateTildes

> data SystemState a = SystemState { stateS     :: a
>                                  , stateI     :: a
>                                  , stateR     :: a
>                                  , stateDelta  :: a
>                                  , stateGamma :: a
>                                  }
>   deriving Show

> deltaT :: Double
> deltaT = 1.0

> stateUpdate :: Particles (SystemState Double) ->
>                Particles (SystemState Double)
> stateUpdate xPrevs =
>   V.zipWith5 SystemState ss is rs deltas gammas
>   where
>     sPrevs     = V.map stateS     xPrevs
>     iPrevs     = V.map stateI     xPrevs
>     rPrevs     = V.map stateR     xPrevs
>     deltaPrevs  = V.map stateDelta  xPrevs
>     gammaPrevs = V.map stateGamma xPrevs
>
>     f b g xs = odeSolve (sirOde b g) xs
>                (0.0 `VS.cons` (deltaT `VS.cons` VS.empty))
>
>     ms :: V.Vector (Matrix Double)
>     ms = V.zipWith3 f deltaPrevs gammaPrevs
>                       (V.zipWith3 (\s i r -> [s, i, r]) sPrevs iPrevs rPrevs)
>
>     ns :: V.Vector (VS.Vector Double)
>     ns = V.map (\m -> (toRows m)!!1) ms
>
>     ss     = V.map (VS.! 0) ns
>     is     = V.map (VS.! 1) ns
>     rs     = V.map (VS.! 2) ns
>     deltas  = deltaPrevs
>     gammas = gammaPrevs



> (.+) :: (Num a) => V.Vector a -> V.Vector a -> V.Vector a
> (.+) = V.zipWith (+)

> stateUpdateNoisy :: MonadRandom m =>
>                     Sym 5 ->
>                     Particles (SystemState Double) ->
>                     m (Particles (SystemState Double))
> stateUpdateNoisy bigQ xPrevs = do
>   let xs = stateUpdate xPrevs
>
>       sPrevs :: V.Vector Double
>       sPrevs     = V.map stateS     xs
>       iPrevs     = V.map stateI     xs
>       rPrevs     = V.map stateR     xs
>       deltaPrevs = V.map stateDelta xs
>       gammaPrevs = V.map stateGamma xs
>
>   let mus :: V.Vector (R 5)
>       mus = V.zipWith5 (\a b c d e -> vector (map log [a, b, c, d, e]))
>                        sPrevs iPrevs rPrevs deltaPrevs gammaPrevs
>
>   nus <- mapM (\mu -> fmap exp $ sample $ rvar (MultivariateNormal mu bigQ)) mus
>
>   let nu1s, nu2s, nu3s, nu4s, nu5s :: V.Vector Double
>       nu1s = V.map (fst . headTail) nus
>       nu2s = V.map (fst . headTail . snd . headTail) nus
>       nu3s = V.map (fst . headTail . snd . headTail . snd . headTail) nus
>       nu4s = V.map (fst . headTail . snd . headTail . snd . headTail .
>                     snd . headTail) nus
>       nu5s = V.map (fst . headTail . snd . headTail . snd . headTail .
>                     snd . headTail . snd . headTail) nus
>
>   return (V.zipWith5 SystemState nu1s nu2s nu3s nu4s nu5s)

> newtype SystemObs a = SystemObs { obsI  :: a }
>   deriving Show

> obsUpdate :: Particles (SystemState Double) ->
>              Particles (SystemObs Double)
> obsUpdate xs = V.map (SystemObs . stateI) xs

> priorMu :: R 5
> priorMu = vector [log initS, log initI, log initR, log {- 0.005 -} delta, log {- 0.4 -} gamma]

> bigP :: Sym 5
> bigP = sym $ matrix [
>     1e-6, 0.0, 0.0, 0.0, 0.0
>   , 0.0, 1e-6, 0.0, 0.0, 0.0
>   , 0.0, 0.0, 1e-6, 0.0, 0.0
>   , 0.0, 0.0, 0.0, 5e-3, 0.0
>   , 0.0, 0.0, 0.0, 0.0, 5e-3
>   ]

> initParticles :: MonadRandom m =>
>                  m (Particles (SystemState Double))
> initParticles = V.replicateM nParticles $ do
>   r <- sample $ rvar (MultivariateNormal priorMu bigP)
>   let x1 = exp $ fst $ headTail r
>       x2 = exp $ fst $ headTail $ snd $ headTail r
>       x3 = exp $ fst $ headTail $ snd $ headTail $ snd $ headTail r
>       x4 = exp $ fst $ headTail $ snd $ headTail $ snd $ headTail $
>            snd $ headTail r
>       x5 = exp $ fst $ headTail $ snd $ headTail $ snd $ headTail $
>            snd $ headTail $ snd $ headTail r
>   return $ SystemState { stateS = x1, stateI = x2, stateR = x3,
>                          stateDelta = x4, stateGamma = x5}

> gens :: [[Double]]
> gens = map toList $ toRows $ tr sol

> obs :: [Double]
> -- obs = [3, 8, 28, 75, 221, 291, 255, 235, 190, 125, 70, 28, 12, 5]
> obs  = gens!!1

> bigR :: Sym 1
> bigR  = sym $ matrix [2.0]

> bigQ :: Sym 5
> bigQ = sym $ matrix
>        [ 1e-4, 0.0, 0.0, 0.0, 0.0
>        , 0.0, 1e-4, 0.0, 0.0, 0.0
>        , 0.0, 0.0, 1e-4, 0.0, 0.0
>        , 0.0, 0.0, 0.0, 1e-3, 0.0
>        , 0.0, 0.0, 0.0, 0.0, 1e-2
>        ]

> weight :: forall a n . KnownNat n =>
>           (a -> R n) ->
>           Sym n ->
>           a -> a -> Double
> weight f bigR obs obsNew = pdf (MultivariateNormal (f obsNew) bigR) (f obs)

> runFilter :: [Particles (SystemState Double)]
> runFilter = snd $ evalState action (pureMT 19)
>   where
>     action = runWriterT $ do
>       xs <- lift $ initParticles
>       V.foldM
>         (oneFilteringStep (stateUpdateNoisy bigQ) obsUpdate (weight f bigR))
>         xs
>         (V.fromList $ map SystemObs obs)

> testFiltering :: [Double]
> testFiltering = map ((/ (fromIntegral nParticles)). sum . V.map stateGamma) runFilter

> testFilteringF :: (([Double], [Double]), [Double])
> testFilteringF = ((s, i), r)
>   where
>     ps = runFilter
>     s = map ((/ (fromIntegral nParticles)). sum . V.map stateS) ps
>     i = map ((/ (fromIntegral nParticles)). sum . V.map stateI) ps
>     r = map ((/ (fromIntegral nParticles)). sum . V.map stateR) ps

> type ParticleStream = [[Double]]

> testFilteringS ::
>   (ParticleStream, (ParticleStream, (ParticleStream, (ParticleStream, ParticleStream))))
> testFilteringS = (s, (i, (r, (d, g))))
>   where
>     ps = runFilter
>     s = map (take 200 . V.toList . V.map stateS) ps
>     i = map (take 200 . V.toList . V.map stateI) ps
>     r = map (take 200 . V.toList . V.map stateR) ps
>     d = map (take 200 . V.toList . V.map stateDelta) ps
>     g = map (take 200 . V.toList . V.map stateGamma) ps

Notes
=====

Helpers for Converting Types
----------------------------

> f :: SystemObs Double -> R 1
> f = vector . pure . obsI


Helpers for the Inverse CDF
---------------------------

That these are helpers for the inverse CDF is delayed to another blog
post.

> indices :: V.Vector Double -> V.Vector Double -> V.Vector Int
> indices bs xs = V.map (binarySearch bs) xs

> binarySearch :: (Ord a) =>
>                 V.Vector a -> a -> Int
> binarySearch vec x = loop 0 (V.length vec - 1)
>   where
>     loop !l !u
>       | u <= l    = l
>       | otherwise = let e = vec V.! k in if x <= e then loop l k else loop (k+1) u
>       where k = l + (u - l) `shiftR` 1

Bibliography
============