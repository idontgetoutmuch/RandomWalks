% Stochastic Integration
% Dominic Steinitz
% 25th July 2015

---
bibliography: Stochastic.bib
---

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
partition of $[0, t]$ then we can discretise the above and add in some
noise which we model as samples of Brownian motion at the selected
times multiplied by $b$ so that we can vary the amount noise depending
on the time and the state.

$$
{x}_{t_{i+1}} - {x}_{t_i} = \sum_{i=0}^{n-1} a({x}_{t_i})(t_{i+1} - t_i) +
                            \sum_{i=0}^{n-1} b({x}_{t_i})(W_{t_{i+1}} - W_{t_i})
$$

where $\xi_n \sim {\cal{N}}(0, \delta t)$. We would like to create a continuous model.
One possibility is to write

$$
X_t = \int_0^t x(s)\mathrm{d}s = \int_0^t a(s)\mathrm{d}s + W_t
$$

where $W_t$ is Brownian motion (which we define below). Even though
Brownian motion is nowhere differentiable, we can adopt a convention
that we can write this as

$$
\mathrm{d}X_t = a_t + \mathrm{d}W_t
$$


Let $S$ an $T$ be stopping times and let the filtration on which they
are defined be right continuous. Then

1. $S \land T = \min(S, T)$,
2. $S \lor T = \max(S, T)$,
3. $S + T$ and
4. $\alpha T$

are stopping times where $\alpha > 1$.

For the first we have $\{S \land T \leq t\} = \{S \leq t\} \cup \{T
\leq t\}$ and both the latter are in ${\mathcal{F}_t}$ by the definition
of a [stopping time](http://en.wikipedia.org/wiki/Stopping_time).

Similarly for the second $\{S \lor T \leq t\} = \{S \leq t\} \cap \{T
\leq t\}$.

For the fourth we have $\{\alpha T \leq t\} = \{T \leq \frac{t}{\alpha}\}
\in {\mathcal{F}}_{t / \alpha} \subset {\mathcal{F}}_{t}$ since $\alpha > 1$.

The third is slightly trickier. For $\omega \in \Omega$, $S(\omega)
+T(\omega) < t$ if and only if for some rational $q$, we have
$S(\omega) + T(\omega) < q < t$. We can thus we can find $r \in
\mathbb{Q}$ such that $S(\omega) < r < q - T(\omega)$. Writing $s
\triangleq q - r \in \mathbb{Q}$ we also have $T(\omega) < q - r =
s$. Thus we have $S(\omega) + T(\omega) < t$ if and only if there
exist $r, s \in \mathbb{Q}$ and $r, s > 0$ such that $r + s < t$ and
$S(\omega) < r$ and $T(\omega) < s$. In other words

$$
\{S +T < t\} = \bigcup_{r, s \in \mathbb{Q}^+} (\{S < r\} \cap \{T < s\}
$$

By right continuity [@protter2004stochastic Theorem 1] of the
filtration, we know the terms on the right hand side are in
${\mathcal{F}}_r \subset {\mathcal{F}}_t$ and ${\mathcal{F}}_s \subset
{\mathcal{F}}_t$ so that the whole right hand side is in
${\mathcal{F}}_t$. We thus know that the left hand side is in
${\mathcal{F}}_t$ and using right continuity again that therefore $S +
T$ must be a stopping time.

Bibliography
============