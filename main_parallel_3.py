#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
### Name: main_parallel_3.py
### Author: L. Capocchi
### Version: 2.0
### Description: script to compute the best partition of a (pykov) ergotic Markov chain using the KL rate. See the __main__ for the usage. 
### Dependencies: pykov, networkx, numpy, tdqm
### Python version: 3.7
### Date: 04/22/2020
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
                #G = nx.DiGraph(list(P.keys()), directed=True)
                #assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergotic!"

                ### number of state
                n = len(S)
                    
                if z is not None:
                    zz = list(map(int, str(z)))
                    assert len(zz)==n, f'The number of coordinates len({z})={len(zz)} must be equal to n ({n})'

                ### parallel processing
                K = range(k,k+1) if k else range(1,n)

                ### Partition object
                partitionObject = Partition(n)
                partitionObject.AddStateLabels(S)

                ### for Windows support
                if platform == "win32":
                    pool = Pool(initializer=tqdm.set_lock, initargs=(Lock(),))
                else: 
                    pool = Pool()

                mapped_values = tqdm(pool.imap_unordered(partial(get_best_partition, partitionObject=partitionObject, S=S, P=P, Pi=Pi, coordinate_choice=z), K), desc='Optimal partitioning process', total=len(K))

                #pool = Pool()
                #mapped_values = pool.imap_unordered(partial(get_best_partition,  partitionObject=partitionObject, S=S, P=P, Pi=Pi, splitting_choice=c,coordinate_choice=z), K)

                ### find best partition by minimizing
                kl = sys.maxsize
                result = None
                for d in mapped_values:
                    for k,v in d.items():
                        if v[0]<kl:
                            kl=v[0]
                            result=v
            
                ### Print section #########################################################
                print("\nStates:",S)
                print(f"Best KL rate: {result[0]: .4f}")                   
                print(f"Best partition: {result[1]} ({partitionObject.GetNumberedPartition(result[1])})")
                #print("Best lumped matrix: ")
                #pprint(result[-1])

                # ### test
                #print("\nSteady of the original: ",P.steady())
                #print("Steady of the lumped: ",result[-1].steady())
                ##########################################################################
                #print("\nMean First Passage Times of P:")

                ###  Mean First Passage Times Analysis ###################################

                ### try to find the k from generalized degree
                #edges_with_weights = [ (s1,s2,v) for s1 in P.states() for s2,v in P.mfpt_to(s1).items()]
                #pprint(sorted(edges_with_weights, key=lambda tup: tup[-1]))
                
                #max_mftp = sorted(edges_with_weights, key=lambda tup: tup[-1])[-1]

                #l = round(max_mftp/k)

                
                #print("\nAggreate Mean First Passage Time of P:")
                #MFTP = get_mfpt(P)
                #pprint(sorted(MFTP, key=lambda tup: tup[-1]))

                #G = nx.Graph()
                #G.add_weighted_edges_from(MFTP)
                
                #l = list(set([v2 for v1 in nx.generalized_degree(G).values() for v2 in v1.values()]))
                #assert(len(l)==1)
                #print("\nk from Generalized_degree:",l[0]) 
                #### two state ca be aggregated
                #if l[0] == n-1:
                #    pass
                
                #print("Mean First Passage Times of Q_mu:",list(get_mfpt(result[-1])))

                ############################################################################

                # end time
                end = time.time()

                # total time taken
                print(f"`\nRuntime of the program is {end - start}")
