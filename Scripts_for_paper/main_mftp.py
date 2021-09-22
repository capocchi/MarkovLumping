#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is used to compute the best partition of a (pykov) ergotic Markov chain using the KL rate. See the __main__ for the usage.
# Copyright (C) 2021  Laurent Capocchi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Email: capocchi@univ-corse.fr

#######################################################################
### Name: main_mftp.py
### Author: L. Capocchi
### Version: 2.0
### Description: script to compute the correlation between Mean First Passage Time and partioning of a Markov chain. See the __main__ for the usage. 
### Dependencies: pykov, networkx, numpy, tdqm
### Python version: 3.9
### Date: 09/22/2021
#######################################################################

import pykov, sys, time
import networkx as nx
import statistics
from more_itertools import locate
from sys import platform
from multiprocessing import freeze_support
import multiprocessing as mp

from Partition import Partition
from Lifting import Lump, KL, Lifting
from nbPartitions import calculSbyGordon

def get_edges_with_weights(P:pykov.Chain)->tuple:
    """
    """  
    for s1 in P.states():
        for s2,v in P.mfpt_to(s1).items():
            yield (s1,s2,v)

__getMFTP_cache = {}
def getMFTP(P:pykov.Chain,s):
    if s not in __getMFTP_cache:
        __getMFTP_cache[s] = P.mfpt_to(s).items()
    return __getMFTP_cache[s]

def getMFTPs(P:pykov.Chain,p)->float:
    L = [v for s in p for s2,v in getMFTP(P,s) if s2 not in p]
    return min([abs(a-b) for a,b in zip(L[::2], L[1::2])])

def getMFTPs2(P:pykov.Chain,p)->float:
    for s in p:
        a = None
        b = None
        r = 0.0
        min = 100000000
        for s2,v in getMFTP(P,s):
            if s2 not in p:
                if a:
                    b=v
                else:
                    a=v
                if a and b:
                    r = abs(a-b)
                    if r < min:
                        min = r
                    a = None
                    b = None
    return min

    
def get_ordered_partitions_from_mftp(S:[str],P:pykov.Chain)->tuple:
    """ Get orderded list of partitions.
    """ 
    n = len(S)
    partitionObject = Partition(n)
    partitionObject.AddStateLabels(S)
    
    dd = {}
    
    for c in partitionObject.GetLabeled(k=n-1):
        listOfLengths = map(len,c)
        ### for all couple of states p / for k=n-1 the lenght of the list is 1
        ### consider only the partition with 2 states
        ### the aggregation of 2 state is the best choice... (only for these kind of matrix ?)
        for i in locate(listOfLengths, lambda a: a == 2):
            p = c[i]
            if p in dd:
               mfpt = getMFTPs(P,p)
               if mfpt < dd[d][-1]:
                   dd[p]=(c,mfpt)
            else:
                dd[p]=(c,getMFTPs(P,p))
                
    ### heristic is based on the mean of mftp values
    mean = statistics.mean([v[1] for v in dd.values()])

    for k,v in dd.items():
        if v[1] < mean:
            yield v[0]

if __name__ == '__main__':

    # for Windows support of tqdm
    if platform == "win32":
        freeze_support()

    ### filename of the P Matrix
    try:
        fn = sys.argv[1]
    except:
        sys.exit()
    else:
        # starting time
        start = time.time()
        
        ### Markov matrix
        try:
            P = pykov.readmat(fn)
        except AttributeError as info:
            P = pykov.Matrix()
            with open(fn,'r') as f:
                for s in f.readlines():
                    s1,s2,p = s.split(' ')
                    P[(s1,s2)]=float(p.strip())
        

        S = tuple(set([a[0] for a in [i.strip('\n').split() for i in open(fn)]]))

        ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
        G = nx.DiGraph(list(P.keys()), directed=True)
        assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergotic!"

        ###  Mean First Passage Times Analysis ###################################
        
        Pi = P.steady()
        count = 0
        ### result varaible - kl = 1.0 is the max value; we find the min.
        result = {'kl':1.0,'partition':None}
        ### loop on partitions to find the best from the mftp analysis
        for p in get_ordered_partitions_from_mftp(S,P):
            partition = {''.join(['NS',str(i)]):a for i,a in enumerate(p)}

            ### compute the kl divergence rate
            Q = Lump(partition, Pi, P)
            p = [(v, k1) for k1,v1 in partition.items() for v in v1]
            Q_mu = Lifting(Q, Pi, S, p)
            kl = KL(S, P, Pi, Q_mu)
        
            ### store the best kl and partition
            if kl < result['kl']:
                result['kl']=kl
                result['partition']=p
        
            count+=1

        # number of states        
        n = len(P.states())

        with mp.Pool(mp.cpu_count()) as pool:
            r = pool.starmap(calculSbyGordon, [(n,k) for k in range(n)])

        print(f"Nb of partition:{sum(r)}")
        print(f"Nb of partition for k=n-1=[{n-1}:{calculSbyGordon(n,n-1)}")
        print(f"Best KL :{result['kl']}")
        print(f"Optimal Partition Max Lenght:{count}") 
        print(f"Best partition :{result['partition']}")

        # end time
        end = time.time()

        # total time taken
        print(f"`\nRuntime of the program is {end - start}")
