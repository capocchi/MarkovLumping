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
python main_matrix.py 4 uniform 4x4_0.1.dat uniform 0.1
```

See the header of the main_matrix.py to have more details concerning the options. 

The Matrix/ directory contains all of the matrix used in the experiments.

## Best partition k=n-1 (Table 1)

![Alt Text](https://media.giphy.com/media/vFKqnCdLPNOKc/giphy.gif)

Table 1 is obtained by executing the folowing command from the Script_for_paper/ directory:

```
python table1.py
```

## Validation of deterministic and heuristic improvements (Table 2)

Table 1 is obtained by executing the folowing command from the Script_for_paper/ directory:

```
python table2.py
```
