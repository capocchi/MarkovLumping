#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
### Name: main_mftp.py
### Author: L. Capocchi
### Version: 2.0
### Description: script to compute the correlation between Mean First Passage Time and partioning of a Markov chain. See the __main__ for the usage. 
### Dependencies: pykov, networkx, numpy, tdqm
### Python version: 3.7
### Date: 04/23/2020
#######################################################################

import pykov, sys, time, random
import networkx as nx
import numpy as np
import statistics

from sys import platform
from multiprocessing import freeze_support
from pprint import pprint

from Partition3 import Partition

def get_mfpt(P:pykov.Chain)->[tuple]:
    """ Get mfpt of Markov chain P as list of tuple.
    """
    t = {(s1,s2):v for s1 in P.states() for s2,v in P.mfpt_to(s1).items()}
    r = {}
    for k,v in t.items():
        a = tuple(sorted(k))
        if a not in r.keys():
            r[a]=v
        else:
            r[a]-=v
        
    return [(k[0],k[1],abs(v)) for k,v in r.items()]

def get_edges_with_weights(P):
    """
    """  
    for s1 in P.states():
        for s2,v in P.mfpt_to(s1).items():
            yield (s1,s2,v)

def get_ordered_partitions(S:[str],P:pykov.Chain)->tuple:
    """ Get orderded list of partitions.
    """
    #edges = sorted(edges_with_weights, key=lambda tup: tup[-1])
    
    n = len(S)
    partitionObject = Partition(n)
    partitionObject.AddStateLabels(S)
    edges = set(get_edges_with_weights(P))

    dd = {}
    ### k/2 is the best choice ?
    for c in partitionObject.GetLabeled(k=n/2):
        for p in c:
            ### si == 2 uniquement que les partitions à deux états
            if (len(p)>=2) and (p not in dd):
                ### TODO ajouter c (pas p)
                dd[p]=[]
                for s in p:
                    for s1,s2,dist in edges:
                        if s == s1 and s2 not in p:
                            dd[p].append(dist)

    ddd = {}
    for k,v in dd.items():
        ddd[k] = min([abs(a-b) for a,b in zip(v[::2], v[1::2])]) 

    mean = statistics.mean(ddd.values())
    for k,v in ddd.items():
        if v < mean:
            #### les couples sont les meilleurs partitions !
            if len(k) == 2:
                yield k

if __name__ == '__main__':

    # for Windows support of tqdm
    if platform == "win32":
        freeze_support()

    ### $ python main_parallel_3.py 4x4.dat 2 1212 
    ### get the KL of the partiton 1212 from 4x4 matrix
    ### $ python main_parallel_3.py 4x4.dat 2
    ### get the best partition with bi-classes (k=2)

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
            with open(fn,'r') as f:
                for s in f.readlines():
                    s1,s2,p = s.split(' ')
                    P[(s1,s2)]=float(p.strip())
        
        ### list of labeled states (dont use S = list(P.states() due to the unordered resulting list)
        S = []
        for a in [i.strip('\n').split() for i in open(fn)]:
            if a[0] not in S:
                S.append(a[0])
        
        ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
        G = nx.DiGraph(list(P.keys()), directed=True)
        assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergotic!"

        ### Best K based on Generalized Degree
        #MFTP = get_mfpt(P)
        #pprint(sorted(MFTP, key=lambda tup: tup[-1]))

        #G = nx.Graph()
        #G.add_weighted_edges_from(MFTP)
        
        #l = list(set([v2 for v1 in nx.generalized_degree(G).values() for v2 in v1.values()]))
        #assert(len(l)==1)
        #k = l[0]
        #print("\nBest k based on Generalized Degree:",k)

        ###  Mean First Passage Times Analysis ###################################
        print(f"\nOrdered list of pairs (best to worst):")
        for p in get_ordered_partitions(S,P):
            print(p)
        
        # end time
        end = time.time()

        # total time taken
        print(f"`\nRuntime of the program is {end - start}")
