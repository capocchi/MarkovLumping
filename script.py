#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess, sys

template = "python Partition.py {} {}"

n = int(sys.argv[1])

args = [[n,i] for i in range(1,n+1)]

processes = []
for arg in args:
    command = template.format(*[str(a) for a in arg])
    process = subprocess.Popen(command)
    processes.append(process)
    
output = [p.wait() for p in processes]
