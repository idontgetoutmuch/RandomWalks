% Dynamical Systems and Bayesian Inference
% Dominic Steinitz
% 28th June 2016

---
bibliography: Stochastic.bib
---

Introduction
============

Dynamical systems and Bayesian inference are two disciplines which on
the surface appear to have very little in common. The purpose of this
blog post is to show these can be related. Of we have to be
restrictive in which Bayesian inference problems we consider and limit
ourselves to state space models while noting that many general
problems can be re-formulated as state space models.

Dynamical System Specification
==============================

Let $X$ be a set (usually taken to be a complete metric space), $T$ a
monoid and $S : T \times X \rightarrow X$ a family of operators, often
written as $S_t(x) \triangleq S(t,x)$.

A [dynamical
system](https://en.wikipedia.org/wiki/Dynamical_system_(definition)#General_definition)
is ${\mathcal{D}}$ is a triple of objects $(X, S, T)$ (@Chueshov1999,
@Temam1997 and @Robinson2003) such that

$$
\begin{aligned}
S_0(x) = x \\
S_t(S_s(x)) = S_{t+s}(x)
\end{aligned}
$$

$T$ is called a semi-group in dynamical systems literature and in
Markov process literature (we will shortly see that a state space
model is defined in terms of Markov processes). We use the designation
monoid to stress that there is a unit element as well as an
associative binary operation.

$S$ is called the evolutionary operator and sometimes the family
$\{S_t}_{t \in T}\}$ is called the flow of the dynamical system.

Examples
--------



Bibliography
============