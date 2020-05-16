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

def get_mfpt(P:pykov.Chain)->[tuple]:
    """ Get mfpt of Markov chain P as list of tuple.
    """
    t = { (s1,s2):v for s1 in P.states() for s2,v in P.mfpt_to(s1).items()}
    r = {}
    for k,v in t.items():
        a = tuple(sorted(k))
        if a not in r.keys():
            r[a]=v
        else:
            r[a]-=v
        
    return [(k[0],k[1],abs(v)) for k,v in r.items()]

def reduce(l:[float])->[float]:
    while(len(l)>=2):
        l = [abs(a-b) for a,b in zip(l[::2], l[1::2])]
    return l

def partiton_gen(S:[str],n:int=2)->tuple:
    """ Get partitions of lenght n as list of tuple.
    """
    for p1 in range(len(S)):
        for p2 in range(p1+(n-1),len(S)):
            yield (S[p1],S[p2])

def get_ordered_partitions(S:[str],P:pykov.Chain, n:int=2)->[tuple]:
    """ Get orderded list of partitions.
    """
    edges_with_weights = [(s1,s2,v) for s1 in P.states() for s2,v in P.mfpt_to(s1).items()]
    edges = sorted(edges_with_weights, key=lambda tup: tup[-1])
    
    d = {}
    for c in edges_with_weights:
        if c[0] not in d:
            d[c[0]] = {'V':[c[-1]], 'S':[c[1]]}
        else:
            d[c[0]]['V'].append(c[-1])
            d[c[0]]['S'].append(c[1]) 
    
    ### compute all couple with distance regarding all other states
    dd = {}
    for i in partiton_gen(S,n):
        dd[i]=[]
        for s in i:
            for s1,s2,dist in edges:
                if s == s1:
                    dd[i].append(dist)

    #ddd = {}
    r = 10000000
    for k,v in dd.items():
        l = [abs(a-b) for a,b in zip(v[::2], v[1::2])]
        var = statistics.variance(l)
        if var < r:
            r = var
            K = k
        #ddd[k]=statistics.variance(l)

    return K
    
    #a1_sorted_keys = sorted(ddd, key=ddd.get, reverse=False)

    #return a1_sorted_keys

    #d = {}
    #for c in edges_with_weights:
    #    num_cls = int(round(c[-1]/interval))
    #    if num_cls == 0:
    #        num_cls = 1
    #    if c[0] not in d:
    #        d[c[0]] = set([num_cls-1])
    #    else:
    #        d[c[0]].add(num_cls-1)

    #print("Possible partitions:")
    #print(d)
    #pprint([a for a in itertools.product(*[list(d[s]) for s in S]) if k-1 in a])

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
        MFTP = get_mfpt(P)
        #pprint(sorted(MFTP, key=lambda tup: tup[-1]))

        G = nx.Graph()
        G.add_weighted_edges_from(MFTP)
        
        l = list(set([v2 for v1 in nx.generalized_degree(G).values() for v2 in v1.values()]))
        assert(len(l)==1)
        k = l[0]
        print("\nBest k based on Generalized Degree:",k)

        ###  Mean First Passage Times Analysis ###################################
        print(f"\nOrdered list of pairs (best to worst): {get_ordered_partitions(S,P,2)}")
        
        # end time
        end = time.time()

        # total time taken
        print(f"`\nRuntime of the program is {end - start}")
