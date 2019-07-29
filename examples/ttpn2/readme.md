# Model

We analyze the time to the first grade 2+ peripheral neuropathy (PN) event in patients treated with an antibody-drug conjugate (ADC) delivering monomethyl auristatin E (MMAE). We will simulate and analyze data using a simplified version of the model reported in <sup id="0d9ff650cd3741dc904e6c83788ba3a2"><a href="#lu_time--event_2017" title="Lu, Gillespie, Girish, Agarwal, Li, Hirata, Chu, Kagedal, Leon, Maiya \&amp; Jin, Time-to-{Event} {Analysis} of {Polatuzumab} {Vedotin}-{Induced} {Peripheral} {Neuropathy} to {Assist} in the {Comparison} of {Clinical} {Dosing} {Regimens}, {CPT: pharmacometrics \&amp; systems pharmacology}, v(6), 401--408 (2017).">lu_time--event_2017</a></sup>.

-   Fauxlatuzumab vedotin 1.2 mg/kg IV boluses q3w \(\times\) 6 does.
-   19 patients with 6 right-censored (simulated data).
-   To keep things simpler, we use the simulated individual CL and V values, and only model PD part of the problem.

<./lu2017Model.pdf>

-   PN hazard is substantially delayed relative to PK exposure.
-   Hazard increases over time to an extent not completely described by PK.

Likelihood for time to first PN \(\ge\) 2 event in the \(i^{th}\) patient:

\begin{align*}
\lefteqn{L\left(\theta | t_{\text{PN},i}, \text{censor}_i, X_i\right)} \\
  &= \left\{ \begin{array}{ll}
     h_i\left(t_{\text{PN},i} | \theta, X_i\right) e^{-\int_0^{t_{\text{PN},i}} h_i\left(u | \theta, X_i\right) du}, &
    \text{censor}_i = 0 \\
     e^{-\int_0^{t_{\text{PN},i}} h_i\left(u | \theta, X_i\right) du}, &
     \text{censor}_i = 1
\end{array} \right.
\end{align*}

where

 \begin{align*}
   t_{\text{PN}} &\equiv \text{time to first PN $\ge$ 2 or right
     censoring event} \\
 \theta &\equiv \text{model parameters} \\
 X &\equiv \text{independent variables / covariates} \\
 \text{censor} &\equiv \left\{ \begin{array}{ll}
     1, & \text{PN $\ge$ 2 event is right censored} \\
     0, & \text{PN $\ge$ 2 event is observed} 
 \end{array} \right.
\end{align*}

One can see the expression

\begin{equation*}
  e^{-\int_0^{t_{\text{PN},i}} h_i\left(u | \theta, X_i\right) du}
\end{equation*}

as the survival function at time \(t\).

-   Hazard of PN grade 2+ based on the Weibull distribution, with drug effect proportional to effect site concentration of MMAE:

\begin{align*}
  h_j(t) &= \beta E_{\text{drug}j}(t)^\beta t^{(\beta - 1)} \\
  E_{\text{drug}j}(t) &= \alpha c_{ej}(t) \\
  c^\prime_{ej}(t) &= k_{e0} \left(c_j(t) - c_{ej}(t)\right).
\end{align*}

Overall ODE system including integration of the hazard function:

\begin{align*}
  x_1^\prime &= -\frac{CL}{V} x_1 \\
  x_2^\prime &= k_{e0} \left(\frac{x_1}{V} - x_2\right) \\
  x_3^\prime &= h(t)
  \end{align*}

where \(x_2(t) = c_e(t)\) and \(x_3(t) = \int_0^t h(u) du\) aka cumulative hazard.


# Build


## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
CXXFLAGS += -isystem /usr/local/include                 # path to MPI library's headers
```


## Build in `cmdstan`

```sh
make ../example-models/ttpn2/ttpn2_group
```


# Run

```sh
mpiexec -n 4 -l ttpn2_group sample num_warmup=500 num_samples=500 data file=ttpn2.data2.R init=ttpn2.init.R
```

\bibliography{ttpn2}
