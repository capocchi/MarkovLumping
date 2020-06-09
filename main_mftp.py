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

import pykov, sys, time
import networkx as nx
import numpy as np
import statistics

from sys import platform
from multiprocessing import freeze_support
from pprint import pprint

from Partition3 import Partition

# def get_mfpt(P:pykov.Chain)->[tuple]:
#     """ Get mfpt of Markov chain P as list of tuple.
#     """
#     t = {(s1,s2):v for s1 in P.states() for s2,v in P.mfpt_to(s1).items()}
#     r = {}
#     for k,v in t.items():
#         a = tuple(sorted(k))
#         if a not in r.keys():
#             r[a]=v
#         else:
#             r[a]-=v
        
#     return [(k[0],k[1],abs(v)) for k,v in r.items()]

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
    #L = []
    for s in p:
        a = None
        b = None
        r = 0.0
        min = 100000000
        #for s1,s2,dist in get_edges_with_weights(P):
        for s2,v in getMFTP(P,s):
            #if s == s1 and s2 not in p:
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
                #dd[p].append(dist)
                #L.append(v)
    return min
    #return min([abs(a-b) for a,b in zip(L[::2], L[1::2])])
    
def get_ordered_partitions(S:[str],P:pykov.Chain)->tuple:
    """ Get orderded list of partitions.
    """
    #edges = sorted(edges_with_weights, key=lambda tup: tup[-1])
    
    n = len(S)
    partitionObject = Partition(n)
    partitionObject.AddStateLabels(S)
    #edges = get_edges_with_weights(P)

    # dd = {}
    # ### k=2 is the best choice ?
    # for c in partitionObject.GetLabeled(k=2):
        
    #     ### consider only the partition with 2 states
    #     for p in c:
    #         ### si == 2 uniquement que les partitions à deux états
    #         ### si >=2 toutes les partitions à plusieurs états
    #         if (len(p)>=2) and (p not in dd):
    #             ### TODO ajouter c (pas p)
    #             dd[p]=getMFTPs(P,p)
        
    #             # dd[p]=[]
    #             # for s in p:
    #             #     #for s1,s2,dist in get_edges_with_weights(P):
    #             #     for s2,v in getMFTP(P,s):
    #             #         #if s == s1 and s2 not in p:
    #             #         if s2 not in p:
    #             #             dd[p].append(dist)

    dd = {}
    ### k=2 is the best choice ?
    for c in partitionObject.GetLabeled(k=n-1):
        length_list = list(map(len,c))
        ### consider only the partition with 2 states
        if 2 in length_list:
            p = c[length_list.index(2)]
            dd[p]=getMFTPs(P,p)

    #print(len(dd))

    mean = statistics.mean(dd.values())
    #print(mean)
    for k,v in dd.items():
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
        
        count = 0
        for p in get_ordered_partitions(S,P):
            print(p)
            count+=1

        print(f"Number of possible best partitions:{count}") 
        
        # end time
        end = time.time()

        # total time taken
        print(f"`\nRuntime of the program is {end - start}")
