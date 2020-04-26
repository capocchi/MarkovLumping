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

from functools import partial
from multiprocessing import Pool, freeze_support, Lock
from tqdm import tqdm
from sys import platform

from pprint import pprint
from Partition3 import Partition
from Lifting import *

def get_mfpt(P):
    t = { (s1,s2):v for s1 in P.states() for s2,v in P.mfpt_to(s1).items()}
    r = {}
    for k,v in t.items():
        a = tuple(sorted(k))
        if a not in r.keys():
            r[a]=v
        else:
            r[a]-=v

    return [(k[0],k[1],abs(v)) for k,v in r.items()]

def get_best_partition(k,partitionObject,S,P,Pi,coordinate_choice):

    result = {k:[sys.maxsize, None, None]}

    #print(next(partitionObject.GetLabeled(k, 'NS', coordinate_choice)))
    #print(next(partitionObject.GetLabeled(k, 'NS', coordinate_choice)))

    ### compute all partitions for a given k (optional coordinate_choice)
    pbar = tqdm(partitionObject.GetLabeled(k, 'NS', coordinate_choice=coordinate_choice), unit='part', position=k, total=None if coordinate_choice else partitionObject.Count(k))
    for p in pbar:
    
        ### one partition for k classes and formated to be computed by Lump and more functions...
        partition = {}
        for a, b in p:
            partition.setdefault(b, []).append(a)

        ### to test
        #partition = {'NS0':"('S')", 'NS1':"('R','C')"}

        par = {str(partition):[(v, k1) for k1,v1 in partition.items() for v in v1]}
        
        ### compute the kl rate
        Q = Lump(partition, Pi, P)
        p = par[str(partition)]
        Q_mu = Lifting(Q, Pi, S, p)
        kl = KL(S, P, Pi, Q_mu)

        #M = 0.0005
        #new_states = ['NS%i'%i for i in range(n-1)]
        #for M in np.linspace(0, 0.1, num=1000):
        #    tetha = np.ones(n)
        #    neta = Neta(p, tetha, M, new_states, S)
        #    print(np.gradient(neta)*kl)
    
        ### store the best kl, the best partition and the lumped matrix Q
        if kl<result[k][0]:
            result[k] = (kl, p, Q)
            pbar.set_description(f"Processing for k={k}, kl={kl :.4f}, par={partitionObject.GetNumberedPartition(p)}")

    return result

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
        ### define the number of classes
        try:
            k = int(sys.argv[2])
        except:
            k = None
        finally:
            ### define cordinate choice
            try:
                z = sys.argv[3]
            except:
                z = None
            finally:

                # starting time
                start = time.time()
                
                ### Markov matrix
                P = pykov.readmat(fn)

                ### Steady
                Pi = P.steady()

                ### list of labeled states (dont use S = list(P.states() due to the unordered resulting list)
                S = []
                for a in [i.strip('\n').split() for i in open(fn)]:
                    if a[0] not in S:
                        S.append(a[0])
                
                ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
                G = nx.DiGraph(list(P.keys()), directed=True)
                assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergotic!"

                ### number of state
                n = len(S)


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
                print("\nMean First Passage Times of P:")

                ### try to find the k from generalized degree
                edges_with_weights = [ (s1,s2,v) for s1 in P.states() for s2,v in P.mfpt_to(s1).items()]
                pprint(sorted(edges_with_weights, key=lambda tup: tup[-1]))
                
                max_mftp = sorted(edges_with_weights, key=lambda tup: tup[-1])[-1][-1]
                interval = max_mftp/k
                
                d = {}
                for c in edges_with_weights:
                    num_cls = int(round(c[-1]/interval))
                    if num_cls == 0:
                        num_cls = 1
                    if c[0] not in d:
                        d[c[0]] = set([num_cls-1])
                    else:
                        d[c[0]].add(num_cls-1)

                print(S)
                print(d)
                
                print("Possible partitions:")
                pprint([a for a in itertools.product(*[list(d[s]) for s in S]) if k-1 in a])

                ############################################################################

                # end time
                end = time.time()

                # total time taken
                print(f"`\nRuntime of the program is {end - start}")
