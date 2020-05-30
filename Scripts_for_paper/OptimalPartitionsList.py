#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from main_mftp import get_ordered_partitions
import pykov

fn = sys.argv[1]

### Markov matrix
try:
    P = pykov.readmat(fn)
except AttributeError as info:
    with open(fn,'r') as f:
        for s in f.readlines():
            s1,s2,p = s.split(' ')
            P[(s1,s2)]=float(p.strip())

### list of labeled states (dont use S = list(P.states() due to the unordered resulting list)
S = []
for a in [i.strip('\n').split() for i in open(fn)]:
    if a[0] not in S:
        S.append(a[0])

count = 0
for c in get_ordered_partitions(S,P):
    count+=1

print(count)
#print(len(list(get_ordered_partitions(S,P))))