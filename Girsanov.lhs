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

Let $M_t$ be a non-negative local martingale. If $\{\tau_n\}_{n \in
\mathbb{N}}$ is a localizing sequence such that $\sup_n \|M_{T \land
\tau_n}\|_p < \infty$ for some $p > 1$ then $M_t$ is a strict
martingale.

**Proof**

$$
\mathbb{E}(|M_T - M_{T \land \tau_n}|) \leq
\mathbb{E}(|M_T - r \land M_T) +
\mathbb{E}(|r \land M_T - r \land M_{T \land \tau_n}|) +
\mathbb{E}(M_{T \land \tau_n} - r \land M_{T \land \tau_n})
$$

By the super-martingale property $\mathbb{E}(M_T) \leq \mathbb{E}(M_0) <
\infty$ and thus by dominated convergence we have that

$$
\lim_{r \rightarrow \infty} \mathbb{E}(r \land M_T) = \mathbb{E}(M_T) \quad \mathrm{and} \quad
\lim_{r \rightarrow \infty}\lim_{n \rightarrow \infty}\mathbb{E}(|r \land M_T - r \land M_{T \land \tau_n}|) = 0
$$

We also have that

$$
\begin{aligned}
\mathbb{E}(M_{T \land \tau_n} - r \land M_{T \land \tau_n}) &=
\mathbb{E}((M_{T \land \tau_n} - r \land M_{T \land \tau_n}){I}_{(M_{T \land \tau_n} > r)}) +
\mathbb{E}((M_{T \land \tau_n} - r \land M_{T \land \tau_n}){I}_{(M_{T \land \tau_n} \leq r)}) \\
&= \mathbb{E}((M_{T \land \tau_n} - r \land M_{T \land \tau_n}){I}_{(M_{T \land \tau_n} > r)}) \\
&= \mathbb{E}(M_{T \land \tau_n}{I}_{(M_{T \land \tau_n} > r)}) - r\mathbb{P}({M_{T \land \tau_n} > r})
\end{aligned}
$$

By Chebyshev's inequality (see note (2) below), we have

$$
r\mathbb{P}({M_{T \land \tau_n} > r}) \leq \frac{r\mathbb{E}|X|^p}{r^p} \leq
\frac{\sup_n{\mathbb{E}(M_{T \land \tau_n})^p}}{r^{p-1}}
$$

Taking limits first over $n \rightarrow \infty$ and then over $r
\rightarrow \infty$ we see that

$$
\lim_{r \rightarrow \infty}\lim_{n \rightarrow \infty} r\mathbb{P}({M_{T \land \tau_n} > r}) \rightarrow 0
$$

For $0 \leq r \leq x$ and $p > 1$ we have $x \leq r^{1-p}x^p$. Thus

$$
\mathbb{E}(M_{T \land \tau_n}{I}_{(M_{T \land \tau_n} > r)}) \leq
r^{1-p}\mathbb{E}(M_{T \land \tau_n}^p{I}_{(M_{T \land \tau_n} > r)}) \leq
r^{1-p}\sup_n(M_{T \land \tau_n}^p)
$$

Again taking limits over $n \rightarrow \infty$ and then over $r
\rightarrow \infty$ we have

$$
\lim_{r \rightarrow \infty}\lim_{n \rightarrow \infty} \mathbb{E}(M_{T \land \tau_n}{I}_{(M_{T \land \tau_n} > r)}) \rightarrow 0
$$

These two conclusions imply

$$
\lim_{r \rightarrow \infty}\lim_{n \rightarrow \infty} \mathbb{E}(M_{T \land \tau_n} - r \land M_{T \land \tau_n}) \rightarrow 0
$$

We can therefore conclude (since $M_{T \land \tau_n}$ is a martingale)

$$
\mathbb{E}(M_T) = \lim_{n \rightarrow \infty}\mathbb{E}(M_{T \land \tau_n}) =
\mathbb{E}(M_0)
$$

Thus by the preceeding lemma $M_t$ is a strict as well as a local martingale.

The Novikov Sufficient Condition
--------------------------------

For any $\mu \in {\cal{L}}^2_{\mathrm{LOC}}[0,T]$ define the local
martingale

$$
M_t(\mu) = \exp{\bigg(\int^t_0 \mu(\omega, s)\,\mathrm{d}W_s -
           \frac{1}{2}\int^0_t \mu^2(\omega, s)\,\mathrm{d}s\bigg)}
$$

and suppose that $\mu$ satisfies the Novikov condition

$$
\mathbb{E}\bigg[\exp{\bigg(\frac{1}{2}\int^T_0 \mu^2(\omega, s)\,\mathrm{d}s\bigg)\bigg]} < \infty
$$

then $M_t(\mu)$ is a strict martingale.

**Proof**

Firstw e note that $M_t(\lambda\mu)$ is a local martingale for $0 <
\lambda < 1$. Let us show that it is a strict martingale. We can
do this if for any localizing sequence $\{\tau_n\}_{n \in \mathbb{N}}$
we can show

$$
\sup_n\mathbb{E}(M_{T \land \tau_n}(\lambda\mu))^p < \infty
$$

using the preceeding lemma where $p > 1$.

We note that

$$
\begin{aligned}
M_t(\lambda\mu) &=
\exp{\bigg(\int^t_0 \lambda\mu(\omega, s)\,\mathrm{d}W_s -
\frac{1}{2}\int^t_0 \lambda^2\mu^2(\omega, s)\,\mathrm{d}s\bigg)} \\
&= {(M_t(\mu))}^{\lambda^2}\exp{\bigg((\lambda - \lambda^2)\int^t_0 \mu(\omega, s)\,\mathrm{d}W_s\bigg)}
\end{aligned}
$$

Now apply Hölder's inequality with conjugates $({p\lambda^2})^{-1}$
and $({1 - p\lambda^2})^{-1}$ where $p$ is chosen to ensure that the
conjugates are both strictly greater than 1 (otherwise we cannot apply
the inequality).

$$
\begin{aligned}
\mathbb{E}((M_t(\lambda\mu))^p)
&=
\mathbb{E}\bigg[{(M_t(\mu))}^{p\lambda^2}\exp{\bigg(p(\lambda - \lambda^2)\int^t_0 \mu(\omega, s)\,\mathrm{d}W_s\bigg)}\bigg] \\
&\le
\bigg|\bigg|{M_t(\mu)}^{p\lambda^2}\bigg|\bigg|_{p\lambda^2}
\bigg|\bigg|\exp{\bigg(p(\lambda - \lambda^2)\int^t_0 \mu(\omega, s)\,\mathrm{d}W_s\bigg)}\bigg|\bigg|_{1 - p\lambda^2} \\
&=
\mathbb{E}{\bigg[M_t(\mu)}\bigg]^{p\lambda^2}
\mathbb{E}\bigg[\exp{\bigg(p\frac{\lambda - \lambda^2}{1 - p\lambda^2}\int^t_0 \mu(\omega, s)\,\mathrm{d}W_s\bigg)}\bigg]^{1 - p\lambda^2}
\end{aligned}
$$

Now let us choose

$$
p\frac{\lambda - \lambda^2}{1 - p\lambda^2} = \frac{1}{2}
$$

then

$$
\begin{aligned}
2p(\lambda - \lambda^2) &= 1 - p\lambda^2 \\
p & = \frac{1}{2(\lambda - \lambda^2) + \lambda^2} \\
p &= \frac{1}{(2 - \lambda)\lambda}
\end{aligned}
$$

In order to apply Hölder's inequality we need to check that
$(p\lambda^2)^{-1} > 1$ and that $(1 - p\lambda^2)^{-1} > 1$ but this
amounts to checking that $p\lambda^2 > 0$ and that $1 > \lambda$. We
also need to check that $p > 0$ but this amounts to checking that $(2
- \lambda)\lambda < 1$ for $0 < \lambda < 1$ and this is easily
checked to be true.


Notes
=====

1. We have already used [importance
sampling](https://idontgetoutmuch.wordpress.com/2014/08/23/importance-sampling/)
and also touched on [changes of
measure](https://idontgetoutmuch.wordpress.com/2015/07/13/conditional-expectation-under-change-of-measure/).

2. [Chebyshev's
inequality](http://mathoverflow.net/questions/28296/analog-of-chebyshevs-inequality-for-higher-moments)
is usually stated for the second moment but the proof is easily
adapted:

$$
\mathbb P( |X| > u ) = \int 1_{|X| > u} ~d\mathbb P = \frac 1 {u^p} \int u^p 1_{|X| > u} ~d\mathbb P < \frac 1 {u^p} \int |X|^p 1_{|X| > u} ~ d\mathbb P \le \frac 1 {u^p} \mathbb E|X|^p.
$$

Bibliography
============
