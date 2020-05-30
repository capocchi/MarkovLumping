#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing
import subprocess
import statistics

def run_command(cmd):
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError:
        output = "ERROR in {}".format(cmd)
    return int(output)

def run_commands(commands, n_parallel=4):
    worker = multiprocessing.Pool(n_parallel)
    print(min(worker.map(run_command, commands)), max(worker.map(run_command, commands)))

if __name__ == "__main__":
    n = 5
    run_commands([
        "python OptimalPartitionsList.py "+os.path.join("..","Matrix","bench","diag_max","uniform","20x20_0.1.dat")]*n, n_parallel=4)