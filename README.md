# MarkovLumping

This README provides the methods to reproduce the results of the paper.

The Scripts_for_paper is the directory that contains all the scripts nedeed to obtain 
the result. The Matrix directory store all of the Matrix used to make the results resumed in table 1 and 2. 

The main_matrix.py Python script can be used to randomly generate into a file a normalized ergotic Markov chain with a size n (number of states) and an distribution drift term (uniform, rayleigh, binomial, weibfull and beta). It generate a stochastic matrix filled with random numbers, given some conditions:

    * The rows have to sum up to 1.
    * The values on the diagonal should be significantly higher than the other values.

For example, if you want to generate uniform matrix with n=4 and high=0.1:

'''
python main_matrix.py 4 uniform 4x4_0.1.dat uniform 0.1
'''

<!--- 
main_matrix.py is the main file that conducts the experiments from the beginning.


The rules and policies are stored in .txt file. Please refer to it when you have
questions or want to have your own rules or policies added.

Reproducing the results:

To generate table 4, simply run
python RemediotMain.py

The results might not be consistent for each run, because we randomize it for
performance issue. But it should not vary significantly.

To generate figure 6 and 7, you need to uncomment evalNumberOfRemedialActions()
in RemediotMain.py and run:
python RemediotMain.py conflict_rules.txt

The raw measurement data are stored in:
https://github.com/nesl/buildsys-19-code/tree/master/results

-->
