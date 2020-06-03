#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is used to compute the best partition of a (pykov) ergotic Markov chain using the KL rate. See the __main__ for the usage.
# Copyright (C) 2020  Laurent Capocchi
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

"""Script documentation.
.. module:: Computes the best partition of a (pykov) ergotic Markov chain using the KL rate.
   :platform: Unix, Windows, Mac
.. moduleauthor::
   Laurent Capocchi <capocchi@univ-corse.fr>
"""
__date__ = 'Mai 2020'

__version__ = 1.0

__license__ = 'GNU General Public License Version 3'

__authors__ = 'Laurent Capocchi'

__many_thanks_to__ = 'Jean-François Santucci'

import pykov, sys, time, random
import networkx as nx
import numpy as np

from functools import partial
from multiprocessing import Pool, freeze_support, Lock
from tqdm import tqdm
from sys import platform

from pprint import pprint
from Partition3 import Partition, partition_k
from Lifting import *
from main_mftp import *

def _get_best_partition(k,partitionObject,S,P,Pi,coordinate_choice):
    """
    """

    result = {k:[sys.maxsize, None, None]}

    #print(next(partitionObject.GetLabeled(k, 'NS', coordinate_choice)))
    #print(next(partitionObject.GetLabeled(k, 'NS', coordinate_choice)))

#    partitions = []
    
    ### compute all partitions for a given k (optional coordinate_choice)
    pbar = tqdm(partitionObject.GetLabeled(k, 'NS', coordinate_choice=coordinate_choice), unit='part', position=k, total=None if coordinate_choice else partitionObject.Count(k))
    for p in pbar:

 #       partitions.append(p)
 
        ### one partition for k classes and formated to be computed by Lump and more functions...
        partition = {}
        for a, b in p:
            partition.setdefault(b, []).append(a)

        ### to test
        #partition = {'NS0':"('S')", 'NS1':"('R','C')"}

        par = {str(partition):[(v, k1) for k1,v1 in partition.items() for v in v1]}

        ### compute the kl rate
        Q = Lump(partition, Pi, P)
        
        ### Q must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
        #G = nx.DiGraph(list(Q.keys()), directed=True)
        #assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergotic!"

        p = par[str(partition)]
        Q_mu = Lifting(Q, Pi, S, p)
        kl = KL(S, P, Pi, Q_mu)
  
        ### store the best kl, the best partition and the lumped matrix Q
        if kl<result[k][0]:
            result[k] = (kl, p, Q)
            pbar.set_description(f"Processing for k={k}, kl={kl :.4f}, par={partitionObject.GetNumberedPartition(p)}")

#    init_state = S[0]
#    M = 0.1
#    N = 100
#    new_states = ['NS%i'%i for i in range(len(S)-1)]
#    tetha = np.ones(len(S))
#    T = Tetha(tetha, init_state, partitions, M, new_states, P, Pi, S, N)
#    print(T[-1])
#    Analyse(dict(zip(S,T[-1])))

    return result

if __name__ == '__main__':
    """
        >>> python main_parallel_3.py 4x4.dat 2 1212
        Processing for k=2, kl=0.4429, par=1212: : 1part [00:00, 705.99part/s]                                               | 0/1 [00:00<?, ?it/s]
        Optimal partitioning process: 100%|██████████████████████████████████████████████████████████████████████████| 1/1 [00:01<00:00,  1.13s/it]
        Processing for k=2, kl=0.4429, par=1212: : 0part [00:00, ?part/s]
        States: ['S', 'R', 'C', 'D']
        Best KL rate:  0.4429
        Best partition: [('S', 'NS1'), ('C', 'NS1'), ('R', 'NS2'), ('D', 'NS2')] (1212)
    """

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
                try:
                    P = pykov.readmat(fn)
                except AttributeError as info:
                    with open(fn,'r') as f:
                        for s in f.readlines():
                            s1,s2,p = s.split(' ')
                            P[(s1,s2)]=float(p.strip())
                    
                Pi = P.steady()

                assert(Pi.sum()!=0)

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

                mapped_values = tqdm(pool.imap_unordered(partial(_get_best_partition, partitionObject=partitionObject, S=S, P=P, Pi=Pi, coordinate_choice=z), K), desc='Optimal partitioning process', total=len(K))

                #pool = Pool()
                #mapped_values = pool.imap_unordered(partial(_get_best_partition,  partitionObject=partitionObject, S=S, P=P, Pi=Pi, splitting_choice=c,coordinate_choice=z), K)

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

                # end time
                end = time.time()

                # total time taken
                print(f"`\nRuntime of the program is {end - start}")

                # ### test
                #print("\nSteady of the original: ",P.steady())
                #print("Steady of the lumped: ",result[-1].steady())

#                ### Best K based on Generalized Degree
#                MFTP = get_mfpt(P)
#                #pprint(sorted(MFTP, key=lambda tup: tup[-1]))

#                G = nx.Graph()
#                G.add_weighted_edges_from(MFTP)
                
#                l = list(set([v2 for v1 in nx.generalized_degree(G).values() for v2 in v1.values()]))
#                assert(len(l)==1)
#                k = l[0]
#                print("\nBest k based on Generalized Degree:",k)

                ###  Mean First Passage Times Analysis ###################################
                #print("\nMean First Passage Times of P:")

                ### try to find the k from generalized degree
                #edges_with_weights = [ (s1,s2,v) for s1 in P.states() for s2,v in P.mfpt_to(s1).items()]
                #edges = sorted(edges_with_weights, key=lambda tup: tup[-1])
                
                #pprint(edges)
                
                # starting time
                start = time.time()

                ### extract the best partiton with the good format ('S1','S2')
                import statistics
                from statistics import mode 
                d = dict(result[1])
                NS = mode(d.values())

                best_partition = tuple(k for k,v in d.items() if v == NS)

                ###  Mean First Passage Times Analysis ###################################
                ### find to know if the best partition is in the list of best pair finded using mftp!
                
                count = 0
                L = get_ordered_partitions(S,P)
                flag = 'not finded...'
                for c in L:
                    s1,s2 = best_partition
                    if s1 in c and s2 in c:
                        flag = 'finded'
                        break

                print(f"\nNew state: {best_partition} ({NS}) {flag} amoung {len(list(get_ordered_partitions(S,P)))} pairs")

                #print(f"\nOrdered list of pairs (best to worst): {list(get_ordered_partitions(S,P))}")

                # end time
                end = time.time()

                # total time taken
                print(f"`\nRuntime of the program is {end - start}")
