# Model

We consider the kinetics of an autocatalytic reaction <sup id="13f5fad7ac685adc4d9e0f6c9690956c"><a href="#robertson_numerical_1966" title="Robertson, Numerical analysis, an introduction, chapitre {The} solution of a set of reaction rate equations, Academic Press (1966).">robertson_numerical_1966</a></sup>. The structure of the reactions is

\begin{align*}
A &\xrightarrow{k_1} B\\
B+B &\xrightarrow{k_2} C + B\\
B+C&\xrightarrow{k_3} C + A,
\end{align*}

where \(k_1\), \(k_2\), \(k_3\) are the rate constants and \(A\), \(B\) and \(C\) are the chemical species involved. The corresponding ODEs are

\begin{align*}
x_1' &= -k_1x_1 + k_3x_2x_3\\
x_2' &=  k_1x_1 - k_2y_2^2 - k_3x_2x_3\\
x_3' &=  k_2y_2^2
\end{align*}

Given \(k_1=0.04, k_2=1.0e4, k_3=3.0e7\), we seek the initial condition for \(x_1(t=0)\), combined with initial condition \(x_2(t=0)=0.0\) and \(x_3(t=0)=0.0\).


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

\bibliography{autocatalysis}
