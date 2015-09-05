% Girsanov's Theorem
% Dominic Steinitz
% 18th August 2015

---
bibliography: Stochastic.bib
---

We previously used importance sampling in the case where we did not
have a sampler available for the distribution from which we wished to
sample. There is an even more compelling case for using importance sampling.

Suppose we wish to estimate the probability of a rare event. For
example, suppose we wish to estimate $\mathbb{P}(X > 5)$ where $X \sim
{\mathcal{N}}(0,1)$. In this case, we can look up the answer
$\mathbb{P}(X > 5) \approx 2.86710^{-7}$. But suppose we couldn't look
up the answer. One strategy that might occur to us is to sample and
then estimate the probability by counting the number of times out of
the total that the sample was bigger than 5.

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

The Novikov Sufficiency Condition
=================================

Let $X \in {\cal{L}}^2_{\mathrm{LOC}}[0,T]$ and further let it
satisfy the [Novikov
condition](https://en.wikipedia.org/wiki/Novikov%27s_condition)

$$
\mathbb{E}\bigg[\exp{\bigg(\frac{1}{2}\int_0^T X^2(s, \omega) \mathrm{d}s\bigg)}\bigg] < \infty
$$

then the process defined by

$$
M_t(X) = \exp{\bigg(\int_0^t X(t, \omega) \mathrm{d}W_s  -
                      \frac{1}{2}\int_0^t X^2(t, \omega) \mathrm{d}s\bigg)}
$$

is a martingale.

**Proof**

Since $M_t$ is a local martingale (FIXME: we haven't defined this yet!), so is

$$
M_t(\sqrt{\alpha} X) = \exp{\bigg(\int_0^t \sqrt{\alpha} X(t, \omega) \mathrm{d}W_s  -
                      \frac{1}{2}\int_0^t \alpha X^2(t, \omega) \mathrm{d}s\bigg)}
$$

for any $0 < \alpha < 1$.

Lemma
-----

Let $M_t$ for $t \in [0,t]$ be a non-negative local martingale then
$M_t$ is a super-martingale and if further $\mathbb{E}M_T =
\mathbb{E}M_0$ then $M_t$ is a strict martingale.

**Proof**

Let $\{\tau_n\}_{n \in \mathbb{N}}$ be a localizing sequence for $M_t$
then for $0 < s < t < T$ and using [Fatou's
lemma](http://math.stackexchange.com/questions/242920/tricks-to-remember-fatous-lemma)
and the fact that the stopped process is a strict martingale

$$
\mathbb{E}(M_t \,|\, {\mathcal{F}_s}) =
\mathbb{E}(\liminf_{n \rightarrow \infty} M_{t \land \tau_m} \,|\, {\mathcal{F}_s}) \leq
\liminf_{n \rightarrow \infty} \mathbb{E}(M_{t \land \tau_m} \,|\, {\mathcal{F}_s}) =
\liminf_{n \rightarrow \infty} M_{s \land \tau_m} = M_s
$$

Thus $M_t$ is a super-martingale and therefore

$$
\mathbb{E}M_T \leq \mathbb{E}M_t \leq \mathbb{E}M_s \leq \mathbb{E}M_0
$$

By assumption we have $\mathbb{E}M_T \leq \mathbb{E}M_0$ thus $M_t$ is
a strict martingale.

Lemma
-----

Let $M_t$ be a non-negative martingale. If $\{\tau_n\}_{n \in
\mathbb{N}}$ is a localizing sequence such that $\sup_n \|M_{T \land
\tau_n}\|_p < \infty$ for some $p > 1$ then $M_t$ is a strict martingale.

**Proof**

$$
\mathbb{E}(|M_T - M_{T \land \tau_n}|) \leq
\mathbb{E}(|M_T - r \land M_T) +
\mathbb{E}(|r \land M_T - r \land M_{T \land \tau_n}|) +
\mathbb{E}(M_{T \land \tau_n} - r \land M_{T \land \tau_n})
$$

By the super-martingale property $\mathbb{E}(M_T) < \mathbb{E}(M_0) <
\infty$ and thus by bounded convergence we have that

$$
\lim_{r \rightarrow \infty} \mathbb{E}(r \land M_T) = \mathbb{E}(M_T) \quad \mathrm{and} \quad
\lim_{r \rightarrow \infty}\lim_{n \rightarrow \infty}\mathbb{E}(|r \land M_T - r \land M_{T \land \tau_n}|) = 0
$$

Notes
=====

We have already used [importance
sampling](https://idontgetoutmuch.wordpress.com/2014/08/23/importance-sampling/)
and also touched on [changes of
measure](https://idontgetoutmuch.wordpress.com/2015/07/13/conditional-expectation-under-change-of-measure/).

Bibliography
============
