\documentclass{beamer}
%include polycode.fmt
%options ghci -fglasgow-exts

%format rho = "\rho"
%format y1 = "y_1"
%format y2 = "y_2"
%format var = "(1 - \rho^2)"
%format oldTheta2 = "{\theta}_2"
%format newTheta1 = "\tilde{\theta}_1"
%format newTheta2 = "\tilde{\theta}_2"
%format oldTau    = "{\tau}"
%format newTau    = "\tilde{\tau}"
%format initTau   = "{\tau}_0"
%format newMu     = "\tilde{\mu}"
%format xBar      = "\bar{x}"
%format x2Sum     = "\sum x_i^2"

\usepackage{media9}

\usepackage[latin1]{inputenc}
\usepackage{listings}
\usepackage{lstbayes}
\usepackage{color}
\usepackage{ulem}

\newcommand {\framedgraphic}[2] {
    \begin{frame}{#1}
        \begin{center}
            \includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{#2}
        \end{center}
    \end{frame}
}

\definecolor{darkgreen}{rgb}{0,0.5,0}

\usepackage{hyperref}
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,
    urlcolor=cyan,
}

\urlstyle{same}

\usetheme{Warsaw}
\setbeamertemplate{footline}{\insertframenumber/\inserttotalframenumber}

\title[HMC / SMC]{Hamiltonian and Sequential Monte Carlo \\ An Ecosystem Example}
\author{Dominic Steinitz}
\date{27th September 2016}
\begin{document}

\begin{frame}
\titlepage
\end{frame}

\begin{frame}{Outline}
  \tableofcontents[subsectionsonly, pausesections]
\end{frame}

\section{Introduction to Markov Chain Methods}

\subsection{Typical Problems}

\begin{frame}[fragile]{Some Application Areas}

Suppose you have a model described by differential equations, e.g.

\begin{itemize}

\item

Epidemiology: susceptible, infected, recovered (SIR)

\item

Pharmokinetics / pharmodynamics (PK / PD)

\item

Ecology: predator / prey

\end{itemize}

\begin{itemize}

\item

How to fit parameters?

\item

How to capture second order effects?

\end{itemize}

\end{frame}

\framedgraphic{Hudson Bay Hares and Lynxes}{diagrams/HaresLynxes.png}

\begin{frame}[fragile]{Lotka and Volterra Original Model}

$$
\begin{aligned}
\frac{\mathrm{d}N_1}{\mathrm{d}t} & = & \rho_1 N_1  - c_1 N_1 N_2 \label{eq2a} \\
\frac{\mathrm{d}N_2}{\mathrm{d}t} & = & c_2 N_1 N_2 - \rho_2 N2 \label{eq2b}
\end{aligned}
$$

\end{frame}

\framedgraphic{A Typical Path}{diagrams/LV.png}

\begin{frame}[fragile]{Structurally Unstable}

\begin{itemize}

\item
So it seems that this might be a good model.

\item
But can hare population really grow without a limit?

\item
Let us add carrying capacities.

\end{itemize}

$$
\begin{aligned}
\frac{\mathrm{d}N_1}{\mathrm{d}t} & = & \rho_1 N_1 \bigg(1 - \frac{N_1}{K_1}\bigg) - c_1 N_1 N_2 \\
\frac{\mathrm{d}N_2}{\mathrm{d}t} & = & -\rho_2 N_2 \bigg(1 + \frac{N_2}{K_2}\bigg) + c_2 N_1 N_2
\end{aligned}
$$

\end{frame}

\framedgraphic{Typical Paths}{diagrams/PPviaLibBi.png}

\begin{frame}[fragile]{Capturing our Lack of Knowledge}

Hare growth parameter

\begin{itemize}

\item
Uncertainty grows over time

\item
Positive

\end{itemize}

$$
\begin{aligned}
\frac{\mathrm{d}N_1}{\mathrm{d}t} & = & \rho_1 N_1 \bigg(1 - \frac{N_1}{K_1}\bigg) - c_1 N_1 N_2 \\
\frac{\mathrm{d}N_2}{\mathrm{d}t} & = & -\rho_2 N_2 \bigg(1 + \frac{N_2}{K_2}\bigg) + c_2 N_1 N_2 \\
\mathrm{d} \rho_1 & = & \rho_1 \sigma_{\rho_1} \mathrm{d}W_t
\end{aligned}
$$

$W_t$ being a Wiener process.

\end{frame}

\begin{frame}[fragile]{Capturing our Lack of Knowledge II}

By It\^{o} we have

$$
\mathrm{d} (\log{\rho_1}) = - \frac{\sigma_{\rho_1}^2}{2} \mathrm{d} t + \sigma_{\rho_1} \mathrm{d}W_t
$$

We can use this to generate paths for $\rho_1$.

$$
\rho_1(t + \Delta t) = \rho_1(t)\exp{\bigg(- \frac{\sigma_{\rho_1}^2}{2} \Delta t + \sigma_{\rho_1} \sqrt{\Delta t} Z\bigg)}
$$

where $Z \sim {\mathcal{N}}(0,1)$.

\end{frame}

\framedgraphic{Log Brownian Paths}{diagrams/LogBrownianPaths.png}

\framedgraphic{A Typical Noisy System Path}{diagrams/LVdata.png}

\subsection{The Stochastic Model}

\framedgraphic{Posterior for $\mu$}{diagrams/PP.png}



\begin{frame}{LibBI}

\href{http://libbi.org}{LibBi} is used for state-space modelling and
Bayesian inference on high-performance computer hardware, including
multi-core CPUs, many-core GPUs (graphics processing units) and
distributed-memory clusters.

Default method: PMMH / PMCMC C. Andrieu, A.D. \& R. Holenstein,
Particle Markov chain Monte Carlo for Efficient Numerical Simulation

\end{frame}



\begin{frame}[fragile]{LibBI Constants}

\lstset{language=Stan,
  basicstyle=\ttfamily\scriptsize,
  keywordstyle=\color{blue}\ttfamily,
  stringstyle=\color{red}\ttfamily,
  commentstyle=\color{darkgreen}\ttfamily,
  breaklines=true
  }

\begin{lstlisting}[language=Stan]
model PP {
  const h         = 0.1;    // time step
  const delta_abs = 1.0e-3; // absolute error tolerance
  const delta_rel = 1.0e-6; // relative error tolerance

  const a  = 5.0e-1 // Hare growth rate - superfluous
                    // for inference but a reminder
                    // of what we should expect
  const k1 = 2.0e2  // Hare carrying capacity
  const b  = 2.0e-2 // Hare death rate per lynx
  const d  = 4.0e-1 // Lynx death rate
  const k2 = 2.0e1  // Lynx carrying capacity
  const c  = 4.0e-3 // Lynx birth rate per hare
\end{lstlisting}
\end{frame}

\begin{frame}[fragile]{LibBI Parameters}

\lstset{language=Stan,
  basicstyle=\ttfamily\scriptsize,
  keywordstyle=\color{blue}\ttfamily,
  stringstyle=\color{red}\ttfamily,
  commentstyle=\color{darkgreen}\ttfamily,
  breaklines=true
  }

\begin{lstlisting}
  state P, Z       // Hares and lynxes
  state ln_alpha   // Hare growth rate (log of)
  obs P_obs        // Observations of hares
  param mu, sigma  // Mean and standard deviation of
                   // hare growth rate
  noise w          // Noise
\end{lstlisting}
\end{frame}

\begin{frame}[fragile]{LibBI Initialisation}

\lstset{language=Stan,
  basicstyle=\ttfamily\scriptsize,
  keywordstyle=\color{blue}\ttfamily,
  stringstyle=\color{red}\ttfamily,
  commentstyle=\color{darkgreen}\ttfamily,
  breaklines=true
  }

\begin{lstlisting}
  sub parameter {
    mu ~ uniform(0.0, 1.0)
    sigma ~ uniform(0.0, 0.5)
  }

  sub proposal_parameter {
     mu ~ truncated_gaussian(mu, 0.02, 0.0, 1.0);
     sigma ~ truncated_gaussian(sigma, 0.01, 0.0, 0.5);
   }

  sub initial {
    P ~ log_normal(log(100.0), 0.2)
    Z ~ log_normal(log(50.0), 0.1)
    ln_alpha ~ gaussian(log(mu), sigma)
  }
\end{lstlisting}
\end{frame}

\begin{frame}[fragile]{LibBI Model}

\lstset{language=Stan,
  basicstyle=\ttfamily\scriptsize,
  keywordstyle=\color{blue}\ttfamily,
  stringstyle=\color{red}\ttfamily,
  commentstyle=\color{darkgreen}\ttfamily,
  breaklines=true
  }

\begin{lstlisting}

  sub transition(delta = h) {
    w ~ normal(0.0, sqrt(h));
    ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg = 'RK4(3)') {
      dP/dt =  exp(ln_alpha) * P * (1 - P / k1) - b * P * Z
      dZ/dt = -d * Z * (1 + Z / k2) + c * P * Z
      dln_alpha/dt = -sigma * sigma / 2 - sigma * w / h
    }
  }

  sub observation {
    P_obs ~ log_normal(log(P), 0.1)
  }
}
\end{lstlisting}
\end{frame}

\framedgraphic{Posteriors}{diagrams/LvPosterior.png}

\begin{frame}{Stan}

\href{http://mc-stan.org}{Stan} implements gradient-based Markov chain
Monte Carlo (MCMC) algorithms for Bayesian inference, stochastic,
gradient-based variational Bayesian methods for approximate Bayesian
inference, and gradient-based optimization for penalized maximum
likelihood estimation.

Stan implements reverse-mode automatic differentiation to calculate
gradients of the model, which is required by HMC, NUTS, L-BFGS, BFGS,
and variational inference. The automatic differentiation within Stan
can be used outside of the probabilistic programming language.

\end{frame}

\begin{frame}[fragile]{Stan Constants}

\lstset{language=Stan,
  basicstyle=\ttfamily\scriptsize,
  keywordstyle=\color{blue}\ttfamily,
  stringstyle=\color{red}\ttfamily,
  commentstyle=\color{darkgreen}\ttfamily,
  breaklines=true
  }

\begin{lstlisting}
data {
  int<lower=1> T;   // Number of observations
  real y[T];        // Observed hares
  real k1;          // Hare carrying capacity
  real b;           // Hare death rate per lynx
  real d;           // Lynx death rate
  real k2;          // Lynx carrying capacity
  real c;           // Lynx birth rate per hare
  real deltaT;      // Time step
}
\end{lstlisting}

\lstset{language=R,
  basicstyle=\ttfamily\scriptsize,
  keywordstyle=\color{blue}\ttfamily,
  stringstyle=\color{red}\ttfamily,
  commentstyle=\color{darkgreen}\ttfamily,
  breaklines=true
  }

\begin{lstlisting}
data=list(T  = length(rdata_PP$P_obs$value),
          y  = rdata_PP$P_obs$value,
          k1 = 2.0e2,
          b  = 2.0e-2,
          d  = 4.0e-1,
          k2 = 2.0e1,
          c  = 4.0e-3,
          deltaT = rdata_PP$P_obs$time[2] -
                   rdata_PP$P_obs$time[1]
          ),

\end{lstlisting}
\end{frame}

\section{Appendices}

\begin{frame}[fragile]{Acknowledgements}
\begin{itemize}
\item \textcolor{blue}{http://education.mrsec.wisc.edu/463.htm}
\item \textcolor{blue}{http://www.inference.phy.cam.ac.uk/itila/book.html}
\item \textcolor{blue}{http://idontgetoutmuch.wordpress.com}
\end{itemize}
\end{frame}

\end{document}
