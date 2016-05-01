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

`nan` beating `-inf` does not sound good.

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

Bibliography
============