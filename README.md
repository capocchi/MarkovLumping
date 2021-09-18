# MarkovLumping

This README provides the methods to reproduce the results of the paper.

RemediotMain.py is the main file that conducts the experiments from the beginning.
However, it has multiple dependencies to be used, which can be tracked through
the code implementations.

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
