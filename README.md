# MarkovLumping

This README provides the methods to reproduce the results of the paper.

# Dependancies

All the scripts are compatible with Python in version 3.x.
Some packages must be installed to execute the scripts. You can use pip in order to install them as follow:

```
pip install -U numpy tqdm more_itertools networkx 
pip install git+git://github.com/riccardoscalco/Pykov@master
```

# Scripts

The Scripts_for_paper/ directory that contains all the scripts nedeed to obtain 
the results. The tables 1 and 2 are based on the matrix stored in Matrix\bench\diag_max\uniform directory. 

## Matrix Generation
The main_matrix.py Python script can be used to randomly generate into a file a normalized ergotic Markov chain with a size n (number of states) and an distribution drift term (uniform, rayleigh, binomial, weibfull and beta). It generate a stochastic matrix filled with random numbers, given some conditions:

    * The rows have to sum up to 1.
    * The values on the diagonal should be significantly higher than the other values.

For example, if you want to generate uniform matrix with n=4 and high=0.1:

```
python main_matrix.py 4 4x4_0.1.dat uniform 0.1
```

See the header of the main_matrix.py to have more details concerning the options. 

The Matrix/ directory contains all of the matrix used in the experiments.

## Best partition for k=n-1 

![best_partition.py python execution trace](https://user-images.githubusercontent.com/233341/175769664-dd6dfdf3-b56b-4bd8-ab78-f884e65b10ba.gif)

The best partition for k=n-1 results presented in the FOP(LP) column of the table 1 are obtained by executing the folowing command from the Script_for_paper/ directory:

```
python best_partition.py
```

### Proof by contradiction of proposition P(n)

The goal is to prove that given an initial n-ordered Markov chain, the optimum partition has $n-1$ elements. This proposition has been called $P(n)$.

The proof proceeds as follows:

Basic step: we first prove that $P(n)$ is true for the first value of $n$, namely, $n=2$.
The proposition $P(2)$ is obvious.

Inductive step: we assume that $P(k)$ is true and prove that, as a consequence of this, $P(k+1)$ is true.

We have to demonstrate that if $(N,P,\pi)$ 
is a given stationary Markov chain n-ordered, there is a partition function $\phi: N \rightarrow M$ 
where $M$
is $n-1$ 
ordered such as for each other partition m-ordered with $m \textless n-1$, 
$\phi'$,
$R^{(\phi)}(P || \widehat{Q}) < R^{(\phi')}(P || \widehat{Q})$.

We have to demonstrate that if $(N,P,\pi)$ 
is a given stationary Markov chain n+1-ordered, there is a partition function $\phi_2: N \rightarrow M$
where $M$ is $n$ 
ordered such as for each other partition m-ordered with $m \textless n$, 
$\phi_2'$, 
$R^{(\phi_2)}(P || \widehat{Q}) < R^{(\phi_2')}(P || \widehat{Q})$.

Let us suppose that it is not true.
So we have $(N,P,\pi)$ 
being a given stationary n+1-ordered Markov chain and for each $n$ ordered partition ($\phi2: N \rightarrow M$
where $M$ is n ordered), there is a m-ordered partition with $m \textless n$, 
$\phi2'$ such as 
$R^{(\phi2)}(P || \widehat{Q}) \geq R^{(\phi2')}(P || \widehat{Q})$.

But in this case it is true when the n-partition $\phi_2$ has the following sets $C_1, C_2, \ldots, C_n$ 
where $C_n$ is a singleton.

Then we have: 

$$\sum_{i,j=1}^{n+1}{\pi_{i}P_{i,j} log (\frac{P_{i,j}}{\widehat{Q}_{i,j}})} 
\geq
\sum_{i,j=1}^{n+1} {\pi_{i}P_{i,j} log (\frac{P_{i,j}}{\widehat{Q}_{i,j}'})},
$$

then,
 
$$
\sum_{i,j=1}^{n} {\pi_{i}P_{i,j} log (\frac{P_{i,j}}{\widehat{Q}_{i,j}})} + log (\frac{P_{n+1,n+1}}{\widehat{Q}_{n+1,n+1}})
\geq 
\sum_{i,j=1}^{n+1} {\pi_{i}P_{i,j} log (\frac{P_{i,j}}{\widehat{Q}_{i,j}'})}. 
$$

But concerning $\phi_2'$ 
we have $m$ 
set for the corresponding partition $(C'_1, C'_2, \ldots, C'_m)$.

We have to considerate two cases depending on the nature of $C'_m$
which is a singleton containing the state 
$P_{n+1,n+1}$.

In the first case it is impossible to have the inequality since the proposition $P(n)$ is true so if $(N,P,\pi)$ 
be a given stationary n-ordered Markov chain, there is a partition function $\phi: N \rightarrow M$ 
where M is $n-1$ ordered such as for each other partition m-ordered with $m \textless n-1$, 
$\phi'$, 
$R^{(\phi)}(P || \widehat{Q}) < R^{(\phi')}(P || \widehat{Q})$ and since

$$
\sum_{i,j=1}^{n} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}}\right)} + log\left(\frac{P_{n+1,n+1}}{\widehat{Q}_{n+1,n+1}}\right)
\geq 
\sum_{i,j=1}^{n+1} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}'}\right)}, 
$$

we can find a partition that contradicts the inequality of Proposition $P(n)$ since,  

$$
R^{(\phi_2')}(P || \widehat{Q}) = \sum_{i,j=1}^{n} {\pi_{i}P_{i,j}log(\frac{P_{i,j}}{\widehat{Q}_{i,j}'})}  + log(\frac{P_{n+1,n+1}}{\widehat{Q}_{n+1,n+1}'}).
$$

In the second case, $C'_m$ is not a singleton so it contains $n+1-m$ states (at least $m=n-2$). 
So the partition have $m$ sets ($C'_1, C'_2, \ldots, C'_m$) and the inequality $R^{\left(\phi_2\right)}\left(P || \widehat{Q}\right) \geq R^{\left(\phi_2'\right)}\left(P || \widehat{Q}\right)$ 
where $\phi_2$ 
is a n-partition having the following sets $C_1, C_2, \ldots, C_n$ 
and $\phi_2'$
is a partition having $(m \textless n)$ sets 
$C_1, C_2, \ldots, C_m$.

So we can write:

$$
\begin{split}
R^{\left(\phi_2\right)}\left(P || \widehat{Q}\right) = \sum_{i,j=1}^{n} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}}\right)} + \\ 
\sum_{i=n+1,j=1}^{n} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}}\right)} + \\ 
\sum_{i=1,j=n+1}^{n} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}}\right)} + \\ 
{\pi_{n+1}P_{n+1,n+1}log\left(\frac{P_{n+1,n+1}}{\widehat{Q}_{n+1,n+1}}\right)}.
\end{split} 
$$

We should have: $R^{\left(\phi_2\right)}\left(P || \widehat{Q}\right) \geq R^{\left(\phi_2'\right)}\left(P || \widehat{Q}\right)$ 
whatever the partition $\phi_2$.

But we have for m, 
$m \textless n-1$:

$$
\begin{split}
R^{\left(\phi_2'\right)}\left(P || \widehat{Q}\right) = \sum_{i,j=1}^{m} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q`}_{i,j}}\right)}  + \\
\sum_{k=1}^{n} \left(\sum_{i=m+k,j=1}^{n+1} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q'}_{i,j}}\right)} +  \\
\sum_{i=1,j=m+k}^{n+1} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q'}_{i,j}}\right)} \right) + \\
{\pi_{n+1}P_{n+1,n+1}log\left(\frac{P_{n+1,n+1}}{\widehat{Q`}_{n+1,n+1}}\right)}.
\end{split} 
$$

So the inequality $R^{\left(\phi_2\right)}\left(P || \widehat{Q}\right) \geq R^{\left(\phi_2'\right)}\left(P || \widehat{Q}\right)$ 
whatever the partition $\phi_2$ can be simplified and becomes Equation 1 :

$$
\begin{split}
R^{\left(\phi_2\right)}\left(P || \widehat{Q}\right) = \sum_{i,j=1}^{m} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}}\right)}+\\
\sum_{k=1}^{n} \left(\sum_{i=m+k,j=1}^{n+1} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}}\right)} + \sum_{i=1,j=m+k}^{n+1} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}}\right)} \right) + \\ 
{\pi_{n+1}P_{n+1,n+1}log\left(\frac{P_{n+1,n+1}}{\widehat{Q}_{n+1,n+1}}\right)},
\end{split} 
$$

and Equation 2 is 

$$
\begin{split}
R^{\left(\phi_2'\right)}\left(P || \widehat{Q}\right) = \sum_{i,j=1}^{m} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q'}_{i,j}}\right)}+\\
\sum_{k=1}^{n} \left(\sum_{i=m+k,j=1}^{n+1} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q'}_{i,j}}\right)} + \sum_{i=1,j=m+k}^{n+1} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q'}_{i,j}}\right)} \right) + \\ 
{\pi_{n+1}P_{n+1,n+1}log\left(\frac{P_{n+1,n+1}}{\widehat{Q'}_{n+1,n+1}}\right)}.
\end{split} 
$$

So for the partition $\phi_2'$ 
($C'_1, C'_2, \ldots, C'_m$) 
corresponding to the classes involved in the partition, we should have
$R^{\left(\phi_2\right)}\left(P || \widehat{Q}\right) \geq R^{\left(\phi_2'\right)}\left(P || \widehat{Q}\right)$ 
whatever the partition $\phi_2$. 

But once again since proposition $P(n)$ is true (for example when $n=m$):

$$
\begin{split}
\sum_{i,j=1}^{m} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q}_{i,j}}\right)} 
\leq
\sum_{i,j=1}^{m} {\pi_{i}P_{i,j}log\left(\frac{P_{i,j}}{\widehat{Q'}_{i,j}}\right)}.
\end{split} 
$$

Furthermore we can notice that:

$$
\begin{split}
{\widehat{Q'}_{i,j}} = \frac{\pi_j}{\sum\limits_{l \in \psi(j)}{\pi_l}}Q_{\phi(i)\phi(j)},
\end{split} 
$$

when $i$ varies from $m+k$ to $n+1$ and $i$ from 1 to $n+1$ as well as when $i$ varies from 1 to $n+1$ and $j$ to $m+k$ to $n+1$ (with $m \textless n-1$). 
 
That implies that all the terms of the sum involved in the Equation 2 are inferior to the terms of the sum involved the Equation 1. So there is a counter example for the proposed proposition (for the partition $\phi_2'=(C'_1, C'_2, \ldots, C'_m)$ 
corresponding to the classes involved in the 
$R^{\left(\phi_2\right)}\left(P || \widehat{Q}\right) \geq R^{\left(\phi_2'\right)}\left(P || \widehat{Q}\right)$ 
whatever the partition $\phi_2$). 
The proof of proposition $P(n+1)$ is done.

## Validation of deterministic and heuristic improvements (Table 1)

![table1 python execution trace](https://user-images.githubusercontent.com/233341/175769986-60410c41-e392-4be0-8b61-c8dccfb3b77e.gif)

Table 1 is obtained by executing the folowing command from the Script_for_paper/ directory:

```
python table1.py
```
## Validation of BESTA algorithm (Table 2)

![table2 python execution trace](https://user-images.githubusercontent.com/233341/175770132-636fdfdc-33eb-47c9-8b64-cd1b9502e4f7.gif)

Table 2 is obtained by executing the folowing command from the Script_for_paper/ directory:

```
python table2.py
```
table2.py scipt depends on the main_loop.py script that present in the header some boolean constants as:
    - PLOT for plotting reduced graph at each step and the trace of the KL value
    - WRITE_FILE for write reduced matrices in files
    - STAT for displaying some state (using pandas)