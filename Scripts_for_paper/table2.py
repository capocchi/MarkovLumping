#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess

def run_command(cmd):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout_iterator = iter(process.stdout.readline, b"")
    for line in stdout_iterator:
        print(line)

if __name__ == "__main__":
    
    print("Column 2 'Best KL'\n")
    pool = "python main_loop.py "+os.path.join("..","Matrix","16x16.dat ")
    run_command(pool)
    
    print("Last Column 'KL Difference against 16x16 initial matrix'\n")
    
    for k,p in zip((3,4,5,6,7,8,9,10,11,12,13,14,15),("1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,3","1,1,1,2,2,2,2,2,3,3,4,4,4,4,4,4","1,1,1,2,2,2,2,3,4,4,4,4,4,5,5,5","1,1,1,2,2,2,2,3,4,4,5,5,5,6,6,6","1,1,1,2,2,3,3,4,5,5,6,6,6,7,7,7",
                                                     "1,1,1,2,2,3,3,4,5,5,6,6,7,8,8,8", "1,1,1,2,2,3,3,4,5,5,6,6,7,8,9,9","1,1,1,2,2,3,3,4,5,5,6,6,7,8,9,10",
                                                     "1,1,1,2,2,3,3,4,5,6,7,7,8,9,10,11","1,1,1,2,2,3,3,4,5,6,7,8,9,10,11,12", "1,1,1,2,2,3,4,5,6,7,8,9,10,11,12,13",
                                                      "1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14","1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15")):
        print(f"\n---------------------------- Matrix 16x16 for k={k}")
        pool = "python main_parallel.py "+os.path.join("..","Matrix","16x16.dat ")+ f"{k} "+f"{p}"
        run_command(pool)