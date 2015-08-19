% Stochastic Integration
% Dominic Steinitz
% 25th July 2015

---
bibliography: Stochastic.bib
---

Introduction
============

Suppose we wish to model a process described by a differential
equation and initial condition

$$
\begin{aligned}
\dot{x}(t) &= a(x, t) \\
x(0) &= a_0
\end{aligned}
$$

But we wish to do this in the presence of noise. It's not clear how do
to this but maybe we can model the process discretely, add noise and
somehow take limits.

Let $\pi = \{0 = t_0 \leq t_1 \leq \ldots \leq t_n = t\}$ be a
partition of $[0, t]$ then we can discretise the above, allow the
state to be random and add in some noise which we model as samples of
Brownian motion at the selected times multiplied by $b$ so that we can
vary the amount noise depending on the state. We change
the notation from $x$ to $X(\omega)$ to indicate that the variable is
now random over some probability space.

$$
\begin{aligned}
{X}(t_{i+1}, \omega) - {X}(t_i, \omega)  &= a({X}(t_i, \omega))(t_{i+1} - t_i) +
                                            b({X}(t_i, \omega))(W(t_{i+1}, \omega) - W(t_i, \omega)) \\
X(t_0, \omega) &= A_{0}(\omega)
\end{aligned}
$$

We can suppress explicit mention of $\omega$ and use subscripts to
avoid clutter.

$$
\begin{aligned}
{X}_{t_{i+1}} - {X}_{t_i}  &= a({X}_{t_i})(t_{i+1} - t_i) +
                              b({X}_{t_i})(W_{t_{i+1}} - W_{t_i}) \\
X(t_0) &= A_{0}(\omega)
\end{aligned}
$$

We can make this depend continuously on time specifying that

$$
X_t = X_{t_i} \quad \mathrm{for} \, t \in (t_i, t_{i+1}]
$$

and then telescoping to obtain

$$
\begin{aligned}
{X}_{t} &= X_{t_0} + \sum_{i=0}^{k-1} a({X}_{t_i})(t_{i+1} - t_i) +
                     \sum_{i=0}^{k-1} b({X}_{t_i})(W_{t_{i+1}} - W_{t_i})
                     \quad \mathrm{for} \, t \in (t_k, t_{k+1}]
\end{aligned}
$$

In the limit, the second term on the right looks like an ordinary
integral with respect to time albeit the integrand is stochastic but
what are we to make of the the third term? We know that Brownian
motion is nowhere differentiable so it would seem the task is
impossible. However, let us see what progress we can make with
so-called simple proceses.

Simple Processes
================

Let

$$
X(t,\omega) = \sum_{i=0}^{k-1} B_i(\omega)\mathbb{I}_{(t_i, t_{i+1}]}(t)
$$

where $B_i$ is ${\cal{F}}(t_i)$-measurable. We call such a process
*simple*. We can then define

$$
\int_0^\infty X_s \mathrm{d}W_s \triangleq \sum_{i=0}^{k-1} B_i{(W_{t_{i+1}} - W_{t_{i+1}})}
$$

So if we can produce a sequence of simple processes, $X_n$ that
converge in some norm to $X$ then we can define

$$
\int_0^\infty X(s)\mathrm{d}W(s) \triangleq \lim_{n \to \infty}\int_0^\infty X_n(s)\mathrm{d}W(s)
$$

Of course we need to put some conditions of the particular class of
stochastic processes for which this is possible and check that the
limit exists and is unique.

We consider the ${\cal{L}}^2(\mu \times \mathbb{P})$, the space of
square integrable functions with respect to the product measure $\mu
\otimes \mathbb{P}$ where $\mu$ is Lesbegue measure on
${\mathbb{R}^+}$ and $\mathbb{P}$ is some given probability
measure. We further restrict ourselves to [progressively
measurable](http://mathoverflow.net/questions/176622/progressively-measurable-vs-adapted)
functions. More explicitly, we consider the latter class of stochastic
processes such that

$$
\mathbb{E}\int_0^\infty X^2_s\,\mathrm{d}s < \infty
$$

Less Simple Processes
=====================

Bounded, Almost Surely Continuous and Progressively Adapted
-----------------------------------------------------------

Let $X$ be a bounded, almost surely continuous and progressively
measurable process which is (almost surely) $0$ for $t > T$ for some
positive constant $T$. Define

$$
X_n(t, \omega) \triangleq X\bigg(T\frac{i}{n}, \omega\bigg) \quad \mathrm{for} \quad T\frac{i}{n} \leq t \lt T\frac{i + 1}{n}
$$

These processes are cleary progressively measurable and by bounded
convergence ($X$ is bounded by hypothesis and $\{X_n\}_{n=0,\ldots}$
is uniformly bounded by the same bound).

$$
\lim_{n \to \infty}\|X - X_n\|_2 = 0
$$

Bounded and Progressively Measurable
------------------------------------

Let $X$ be a bounded and progressively measurable process which is
(almost surely) $0$ for $t > T$ for some positive constant $T$. Define

$$
X_n(t, \omega) \triangleq \frac{1}{1/n}\int_{t-1/n}^t X(s, \omega) \,\mathrm{d}s
$$

Then $X^n(s, \omega)$ is bounded, continuous and progressively
measurable and it is well known that $X^n(t, \omega) \rightarrow X(t,
\omega)$ as $n \rightarrow 0$. Again by bounded convergence

$$
\lim_{n \to \infty}\|X - X_n\|_2 = 0
$$

Progressively Measurable
------------------------

Firstly, let $X$ be a progressively measurable process which is (almost
surely) $0$ for $t > T$ for some positive constant $T$. Define $X_n(t,
\omega) = X(t, \omega) \land n$. Then $X_n$ is bounded and by
dominated convergence

$$
\lim_{n \to \infty}\|X - X_n\|_2 = 0
$$

Finally let $X$ be a progressively measurable process. Define

$$
X_n(t, \omega) \triangleq
\begin{cases}
X(t, \omega) & \text{if } t \leq n \\
0            & \text{if } \mathrm{otherwise}
\end{cases}
$$

Clearly

$$
\lim_{n \to \infty}\|X - X_n\|_2 = 0
$$

The Itô Isometry
====================

Let $X$ be a simple process such that

$$
\mathbb{E}\int_0^\infty X^2_s\,\mathrm{d}s < \infty
$$

then

$$
\mathbb{E}\bigg(\int_0^\infty X_s\,\mathrm{d}W_s\bigg)^2 =
\mathbb{E}\bigg(\sum_{i=0}^{k-1} B_i{(W_{t_{i+1}} - W_{t_{i}})}\bigg)^2 =
\sum_{i=0}^{k-1} \mathbb{E}(B_i)^2({t_{i+1}} - {t_{i}}) =
\mathbb{E}\int_0^\infty X^2_s\,\mathrm{d}s
$$

Now suppose that $\{H_n\}_{n \in \mathbb{N}}$ is a Cauchy sequence of
progressively measurable simple functions in ${\cal{L}}^2(\mu \times
\mathbb{P})$ then since the difference of two simple processes is again a simple process we can apply the Itô Isometry to deduce that

$$
\lim_{m,n \to \infty}\mathbb{E}\bigg(\int_0^\infty (H_n(s) - H_m(s))\,\mathrm{d}W(s)\bigg)^2 =
\lim_{m,n \to \infty}\mathbb{E}\int_0^\infty (H_n(s) - H_m(s))^2\,\mathrm{d}s =
0
$$

In other words, $\int_0^\infty H_n(s)\,\mathrm{d}W(s)$ is also Cauchy
in ${\cal{L}}^2(\mathbb{P})$ and since this is complete, we can conclude that

$$
\int_0^\infty X(s)\mathrm{d}W(s) \triangleq \lim_{n \to \infty}\int_0^\infty X_n(s)\mathrm{d}W(s)
$$

exists (in ${\cal{L}}^2(\mathbb{P})$). Uniqueness follows using the
triangle inequality and the Itô isometry.

Notes
=====

 1. We defer proving the definition also makes sense almost surely to
another blog post.

 2. This approach seems fairly standard see for example @VonHandel2007
and @morters2010brownian.

 3. @rogers2000diffusions takes a more general approach.

 4. @protter2004stochastic takes a different approach by defining
stochastic processes which are good integrators, a more abstract
motivation than the one we give here.

5. The requirement of progressive measurability can be relaxed.

Bibliography
============