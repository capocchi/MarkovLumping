# MarkovLumping

This README provides the methods to reproduce the results of the paper.

# Dependancies

All the scripts are compatible with Python in version 3.x.
Some packages must be installed to execute the scripts. You can use pip in order to install them as follow:

```
pip install -U numpy tqdm more_itertools networkx pandas matplotlib
```

# Scripts

The Scripts_for_paper/ directory that contains all the scripts nedeed to obtain 
the results. The tables 1 and 2 are based on the matrix stored in Matrix\bench\diag_max\uniform directory. 

## Matrix Generation
The *main_matrix.py* Python script can be used to randomly generate into a file a normalized ergotic Markov chain with a size n (number of states) and an distribution drift term (uniform, rayleigh, binomial, weibfull and beta). It generate a stochastic matrix filled with random numbers, given some conditions:

* The rows have to sum up to 1.
* The values on the diagonal should be significantly higher than the other values.

For example, if you want to generate uniform matrix with n=4 and high=0.1:

```
python main_matrix.py 4 4x4_0.1.dat uniform 0.1
```

See the header of the *main_matrix.py* to have more details concerning the options. 

The Matrix/ directory contains all of the matrix used in the experiments.

## Best partition for k=n-1 

![best_partition.py python execution trace](https://user-images.githubusercontent.com/233341/175769664-dd6dfdf3-b56b-4bd8-ab78-f884e65b10ba.gif)

The best partition for k=n-1 results presented in the **FOP(LP)** column of the table 1 are obtained by executing the folowing command from the Script_for_paper/ directory:

```
python best_partition.py
```

Results can be resumed inside the following table where the last row is reported in the Table 1 as a column.

| **k/n**                        | **3**  | **4**  | **5**  | **6**  | **7**  | **8**  | **9**  | **10** | **11** | **12** | **13** | **14** |
|--------------------------------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| 1                              | 0.465  | 0.8911 | 0.8902 | 0.9003 | 1.0478 | 1.1119 | 1.1016 | 1.0571 | 1.1444 | 1.0176 | 1.0868 | 1.0739 |
| 2                              | **0.0927** | 0.4497 | 0.5242 | 0.5591 | 0.7028 | 0.7588 | 0.7887 | 0.7797 | 0.8619 | 0.785  | 0.8406 | 0.8464 |
| 3                              |        | **0.1679** | 0.3008 | 0.3618 | 0.4699 | 0.523  | 0.6014 | 0.5954 | 0.6677 | 0.6348 | 0.6798 | 0.6923 |
| 4                              |        |        | **0.1403** | 0.197  | 0.319  | 0.3642 | 0.4409 | 0.4696 | 0.5273 | 0.5153 | 0.5526 | 0.5743 |
| 5                              |        |        |        | **0.0829** | 0.1792 | 0.239  | 0.3384 | 0.3596 | 0.4242 | 0.4214 | 0.4556 | 0.4345 |
| 6                              |        |        |        |        | **0.0732** | 0.1494 | 0.2392 | 0.2649 | 0.3318 | 0.3408 | 0.375  | 0.3567 |
| 7                              |        |        |        |        |        | **0.0694** | 0.1483 | 0.1862 | 0.2486 | 0.5153 | 0.3046 | 0.3011 |
| 8                              |        |        |        |        |        |        | **0.0622** | 0.1111 | 0.1812 | 0.2135 | 0.2399 | 0.2761 |
| 9                              |        |        |        |        |        |        |        | **0.0524** | 0.1168 | 0.1555 | 0.1856 | 0.229  |
| 10                             |        |        |        |        |        |        |        |        | **0.0583** | 0.0965 | 0.1365 | 0.1731 |
| 11                             |        |        |        |        |        |        |        |        |        | **0.0466** | 0.0886 | 0.1274 |
| 12                             |        |        |        |        |        |        |        |        |        |        | **0.0435** | 0.0818 |
| 13                             |       |        |        |        |        |        |        |        |        |        |        | **0.0396** |
| **FOP(LP)** Time[s] | **1.0987** | **1.1106** | **1.0815** | **1.1545** | **1.465**  | **3.295**  | **14.974** | **93.988** | **707.23** | **4746.9** | **35588**  | **>10e3**  |

The above table confirms that the deterministic improvement is efficient by pointing out that the optimum solutions (best KL in bold) are belonging to partitions having k=n-1 elements if the initial Markov chain was n-order using a eleven benchmark n-ordered Markov chains uniformly distributed with 0.1 as high parameter.

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

**Table 1** is obtained by executing the folowing command from the Script_for_paper/ directory:

```
python table1.py
```

The execution of the *table1.py* Python script allows to build the following table which gives elements to validate the efficiency of both the deterministic and heuristic improvements on large normalized ergotic Markov chains uniformly distributed with high=0.1.


| **Uniform Matrix (n x n)** | **Nb. of partitions** | **Nb. of partitions for k=n-1** | **Time for k=n-1 [s]** | **Best KL** | **Optimal partition max length** | **FOP(LP) Time[s]** | **{REDUCED(LP) Time [s]** | **FOP(REDUCED(LP)) Time [s]** | **Reduction Rate (compared to k=n-1) [\%]** |
|----------------------------|-----------------------|---------------------------------|------------------------|-------------|----------------------------------|---------------------|---------------------------|-------------------------------|---------------------------------------------|
| 3x3                        | 5                     | 3                               | 1.091                  | 0.0927      | 2                                | 1.0987              | 6.9e^{-7}                 | 0.0028                        | 33.33                                       |
| 4x4                        | 14                    | 5                               | 1.1018                 | 0.1679      | 3                                | 1.1106              | 6.9e^{-7}                 | 0.0043                        | 40.00                                       |
| 5x5                        | 49                    | 8                               | 1.083                  | 0.1403      | 7                                | 1.0815              | 8.0e^{-7}                 | 0.0083                        | 12.50                                       |
| 6x6                        | 198                   | 13                              | 1.1003                 | 0.0829      | 2                                | 1.1545              | 7.0e^{-7}                 | 0.014                         | 7.69                                        |
| 7x7                        | 869                   | 17                              | 1.0891                 | 0.0732      | 14                               | 1.465               | 6.0e^{-7}                 | 0.020                         | 17.64                                       |
| 8x8                        | 4130                  | 24                              | 1.1319                 | 0.0694      | 18                               | 3.295               | 8.0e^{-7}                 | 0.046                         | 25.00                                       |
| 9x9                        | 2110                  | 32                              | 1.1628                 | 0.0622      | 24                               | 14.974              | 7.0e^{-7}                 | 0.053                         | 25.00                                       |
| 10x10                      | 11600                 | 40                              | 1.2344                 | 0.0524      | 30                               | 93.988              | 8.0e^{-7}                 | 0.092                         | 25.00                                       |
| 15x15                      | 1.38e^{+09}           | 99                              | 1.642                  | 0.037       | 68                               | 707.23              | 8.0e^{-7}                 | 0.390                         | 31.31                                       |
| 20x20                      | 5.17e^{+13}           | 181                             | 2.905                  | 0.0279      | 138                              | 4746.9              | 8.9e^{-7}                 | 1.485                         | 23.75                                       |
| 30x30                      | 8.47e^{+23}           | 422                             | 9.9782                 | 0.0148      | 259                              | 35588               | 8.0e^{-7}                 | 6.777                         | 38.62                                       |
| 40x40                      | 1.57e^{+35}           | 762                             | 31.6959                | 0.0115      | 496                              | >10e{^3}            | 1.5$e^{-6}                | 19.526                        | 34.90                                       |
| 50x50                      | 1.86e^{+47}           | 1200                            | 77.3586                | 0.0077      | 799                              | -                   | 8.0e^{-7}                 | 45.942                        | 33.41                                       |
| 60x60                      | 9.77e^{+59}           | 1740                            | 165.386                | 0.0046      | 1192                             | -                   | 9.0e^{-7}                 | 108.486                       | 31.49                                       |
| 70x70                      | 1.81e^{+73}           | 2380                            | 312.5046               | 0.004       | 1508                             | -                   | 8.9e^{-7}                 | 220.963                       | 36.63                                       |
| 80x80                      | 9.91e^{+86}           | 3120                            | 551.0799               | 0.0036      | 2049                             | -                   | 8.9e^{-7}                 | 361.596                       | 34.33                                       |
| 90x90                      | 1.42e^{+101}          | 3960                            | 883.8726               | 0.003       | 2628                             | -                   | 1.0e^{-6}                 | 558.582                       | 33.63                                       |
| 100x100                    | 4.76e^{+115}          | 4900                            | 1391                   | 0.004       | 3214                             | -                   | 9.9e^{-7}                 | 962.731                       | 34.40                                       |

## Validation of BESTA algorithm (Table 2)

![table2 python execution trace](https://user-images.githubusercontent.com/233341/175770132-636fdfdc-33eb-47c9-8b64-cd1b9502e4f7.gif)

**Table 2** is obtained by executing the folowing command from the Script_for_paper/ directory:

```
python table2.py
```
*table2.py* Python scipt depends on the *main_loop.py* Python script that presents in the header some boolean constants as:
* PLOT for plotting reduced graph at each step and the trace of the KL value
* WRITE_FILE for writing reduced matrices in files
* STAT for displaying some statistics (using pandas)