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
\dot{x}_t &= a(x_t) \\
x_0 &= a_0
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
now random over some probability space although we suppress explicit
mention of $\omega$ to avoid clutter.

$$
\begin{aligned}
X_{t_0} &= X_{0} \\
{X}_{t_{i+1}} - {X}_{t_i} &= a({X}_{t_i})(t_{i+1} - t_i) +
                             b({X}_{t_i})(W_{t_{i+1}} - W_{t_i})
\end{aligned}
$$

We can make this depend continuously on time specifying that

$$
X_t = X_{t_i} \quad \mathrm{for} \, t \in (t_i, t_{i+1}]
$$

and then telescoping to obtain

$$
\begin{aligned}
{X}_{t} &= X_0 + \sum_{i=0}^{k-1} a({X}_{t_i})(t_{i+1} - t_i) +
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

where $B_i$ is ${\cal{F}}(t_i)$-measurable. We call such a proces
*simple*. We can then define

$$
\int_0^\infty X(s) \mathrm{d}W(s) \triangleq \sum_{i=0}^{k-1} B_i{(W(t_{i+1}) - W(t_{i+1}))}
$$

So if we can produce a sequence of simple processes, $X_n$ that
converge in some norm to $X$ then we can define

$$
\int_0^\infty \triangleq \lim_{n \to \infty}\int_0^\infty X_n(s)\mathrm{d}W(s)
$$

Bibliography
============