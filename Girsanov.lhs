% Girsanov's Theorem
% Dominic Steinitz
% 18th August 2015

---
bibliography: Stochastic.bib
---

If we have a probability space $(\Omega, {\mathcal{F}}, \mathbb{P})$ and a
non-negative random variable $Z$ then we can define a new probability
measure $\mathbb{Q}$ on the same $\sigma$-algebra by
$$
\mathbb{Q} A \triangleq \int_A Z \,\mathrm{d} \mathbb{P}
$$
For any two probability measures when such a $Z$ exists, it is called
the Radon-Nikod\'ym derivative of $\Q$ with respect to $\P$ and
denoted $\frac{\dif \Q}{\dif \P}$

An Example
----------

Let $\xi$ be a normally distributed random variable with zero mean and
unit variance under a probability measure $\P$ and let $\xi' = a + \xi$. Then if we set
$$
\zeta = e^{-a\xi - \frac{a^2}{2}}
$$
and calculate
\begin{align*}
\Q \{\xi' \le b\} &= \int_{\{\xi' \le b\}} \zeta \dif\P \\
                  &= \int_{\{\xi \le b - a\}} e^{-\xi a -
                  \frac{a^2}{2}} \dif\P \\
                  &= \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{b - a} e^{-xa -
                  \frac{a^2}{2}} e^{\frac{-x^2}{2}} \dif x \\
                  &= \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{b - a}
                  e^{-(x+a)^2} \dif x \\
                  &= \frac{1}{\sqrt{2\pi}} \int_{-\infty}^b
                  e^{-y^2} \dif y
\end{align*}
We have that $\xi'$ is a normal random variable with zero mean and unit
variance under the probability measure $\Q$.

This suggests that we ought to be able to \lq\lq shift\rq\rq\, Brownian
Motion with a drift under a probability measure $\P$ to be pure
Brownian Motion under another probability measure $\Q$.

Girsanov Theorem Statement
--------------------------

Let $W_t$ be Brownian Motion on a probability space $(\Omega, \calF,
\P)$ and let $\{\calF_t\}_{t \in [0,T]}$ be a filtration for this
Brownian Motion and let $\gamma_t$ be an adapted process such that
$$
\EX{\exp (\frac{1}{2}\int_0^T \gamma_t^2\dif t)} < \infty 
$$
then there exists a probability measure $\Q$ such that
\begin{enumerate}
\item
$\Q$ is equivalent to $\P$;
\item
$\displaystyle {\frac{\dif \Q}{\dif \P} = \exp \Bigg(-\int_0^T \gamma_t \dif W_t -
\frac{1}{2} \int_0^T \gamma^2_t \dif t\Bigg)}$;
\item
$\tilde W_t = W_t + \int_0^t \gamma_s \dif s$ is Brownian Motion on the
probabiity space $(\Omega, \calF, \Q)$ also with the filtration $\{\calF_t\}_{t \in [0,T]}$.
\end{enumerate}
