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
        
    for i in (3,4,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50,60,70,80,90,100):
        print(f"\n---------------------------- Uniform Matrix {i}x{i}_0.1")
        pool = "python main_mftp.py "+os.path.join("..","Matrix","bench","diag_max","uniform",f"{i}x{i}_0.1.dat")
        run_command(pool)

        # print(f"\nTIME FOR ALL K")
        # os.system("python main_parallel.py "+os.path.join("..","Matrix","bench","diag_max","uniform",f"{i}x{i}_0.1.dat"))

        # print(f"\nTIME FOR ALL K=N-1")
        # os.system("python main_parallel.py "+os.path.join("..","Matrix","bench","diag_max","uniform",f"{i}x{i}_0.1.dat")+" "+str(i-1))
