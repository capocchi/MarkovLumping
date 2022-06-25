#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess

def run_command(cmd):
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError:
        output = f"ERROR in {cmd}"
    return output

if __name__ == "__main__":
    for i in range(3,14):
        print(f"\n---------------------------- Matrix {i}x{i}_0.1")
        pool = "python main_parallel.py "+os.path.join("..","Matrix","bench","diag_max","uniform",f"{i}x{i}_0.1.dat")
        run_command(pool)