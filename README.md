# TRIUMF-2019-Summer
This is from my work on 'Cosmic Phase Transitions as a Source of Gravitational Waves'.

My code(s) can take in the parameters of the Theory of a given field undergoing a first-order phase transition, and return the parameters conventionally used to find the gravitational wave power spectra.

To do this, the Bubble (or Bounce) approximation was used i.e. the Stationary Phase approximation was applied to the Path Integral of the field tunneling from False to True Vacuum, and the **Instanton Solution** found this way, was used.

A preliminary analysis of a typical 2-3-4 Potential was performed, using the WKB Approximation; this was primarily to build intuition and familiarize oneself with different aspects of the problem at hand.


### Code descriptions:
\begin{enumerate}
\renewcommand{\labelenumii}{\arabic{enumii}}
    \item{.} WKB Approximation for quartic potential
    \item{.} Bounce Approximation for static quartic potential
    \item{.} Bounce Approximation for general static potential
    \item{.} 2-3-4 Potential
    \begin{itemize}
        \item{.} Var\textunderscore Quar\textunderscore Pot\textunderscore Act
        \item{.} S vs. alpha - comparisons        
        \item{.} Act\textunderscore Temp
        \item{.} S vs. T - brute force method
        \item{.} S vs. alpha - template
        \item{.} Act\textunderscore Temp\textunderscore quick
        \item{.} S vs. T - quick method
        \item{.} Tunneling Rate vs. Hubble Rate
    \end{itemize}
    \item{.} 2-4-6 Potential
\end{enumerate}