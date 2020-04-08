#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pykov, sys, time, random
from lxml import etree
import networkx as nx
from functools import partial
from multiprocessing import Pool, freeze_support, Lock
from tqdm import trange, tnrange, tqdm

from pprint import pprint
from Partition import Partition
from nbPartitions import calculSbyGordon
from Lifting import *

def do_work(k,P):

    ### Steady
    Pi = P.steady()

    ### list of labeled states
    S = list(P.states())

    ### number of state
    n = len(S)

    ### find a partition randomly with k classes
    comb = Partition(n)
    comb.AddStateLabels(S)

    ### number of total partition per 1000 in order to have a certitute of exploring max of partitions for k classes
    N = calculSbyGordon(n,k)*1000
    #print('\nProcessing with k=%d for %d loop...'%(k,N))
    
    result = {k:[sys.maxsize, None, None]}
    
    text = "processing for k={}".format(k)
    for i in trange(N, desc=text, position=k-1):

        ### one partition for k classes generated randomly
        for p in comb.GetLabeled(k, 'NS'):
        
            ### one partition for k classes and formated to be computed by Lump and more functions...
            partition = {}
            for a, b in p:
                partition.setdefault(b, []).append(a)

            ### to test
            #partition = {'NS0':"('S')", 'NS1':"('R','C')"}

            par = {str(partition):[(v, k1) for k1,v1 in partition.items() for v in v1]}
            
            #print("Partition",par)
            
            Q = Lump(partition, Pi, P)
            p = par[str(partition)]
            Q_mu = Lifting(Q, Pi, S, p)
            
            #M = 0.0005
            #new_states = ['NS%i'%i for i in range(n-1)]
            #for M in np.linspace(0, 0.1, num=1000):
            #    tetha = np.ones(n)
            #    neta = Neta(p, tetha, M, new_states, S)
            #    print(np.gradient(neta)*kl)    
            
            init_state = random.choice(S)
            kl = KL(init_state, Q_mu, P, Pi, S)

            if kl<result[k][0]:
                result[k] = [kl, list(p), Q]
    
    return result

if __name__ == '__main__':

    # for Windows support
    freeze_support()

    ### $ python PartitionsSchemes 3x3.txt 2
    ### return the lits of partitioning schemes for 3x3 Markov chain and a partitioning in k classes
    
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

            # starting time
            start = time.time()
            
            ### Markov matrix
            P = pykov.readmat(fn)

            #Current=None

            #while(P!=Current):
            ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
            G = nx.DiGraph(list(P.keys()))
            assert(nx.is_strongly_connected(G) and nx.is_aperiodic(G))

            n = len(P.states())

            ### parallel processing
            K = range(k,k+1) if k else range(1,n)
            #pool = Pool()

            # for Windows support
            pool = Pool(initializer=tqdm.set_lock, initargs=(Lock(),))
            mapped_values = tqdm(pool.imap_unordered(partial(do_work, P=P), K), total=len(K))

            kl=10000
            result = None
            for d in mapped_values:
                for k,v in d.items():
                    if v[0]<kl:
                        result=v
                
            #    Current = P
            #    P = result[-1]

            print("KL divergence: %f"%result[0])
            # #print(neta)
            #print("\nBest M parameter: %f"%result[1])
            print("Best partition: ",result[1])
            print("Best lumped matrix: ")
            pprint(result[-1])

            # ### test
            print("Steady of the original: ",P.steady())
            print("Steady of the lumped: ",result[-1].steady())
            
            # end time
            end = time.time()

            # total time taken
            print(f"`\nRuntime of the program is {end - start}")
