#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pykov, sys, time, random
from lxml import etree
import networkx as nx

from pprint import pprint
from Partition import Partition
from nbPartitions import calculSbyGordon
from Lifting import *

if __name__ == '__main__':

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

            ### parameter of the lifting
            #M = .0005
            
            ### Markov matrix
            P = pykov.readmat(fn)

            #tree = etree.parse(os.path.join("C:\\","Users","Laurent","Dropbox","devsimpy","py2x","Domain","Markov","xml","LumpedAggregationComputer.xml"))
            #startStates = [ StartState.text for StartState in tree.xpath("/ProbDEVS/TransitionInfo/StartState")]
            #endStates = [ EndState.text for EndState in tree.xpath("/ProbDEVS/TransitionInfo/EndState")]
            #probValue = [ float(ProbValue.text) for ProbValue in tree.xpath("/ProbDEVS/TransitionInfo/ProbValue")]
        
            ### list of states
            #S = list(set(startStates+endStates))

            ### assoicate new state like 'S0', 'S1', etc.
            #D={}
            #i=0
            #for s in S:
            #    D[s] = 'S'+str(i)
            #    i += 1

            #assert(len(S)==len(D))

            #S = D.values()
            #startStates = [D[ss] for ss in startStates]
            #endStates = [D[ss] for ss in endStates]

            #P = pykov.Chain(dict(zip(zip(startStates,endStates),probValue)))
            
            G = nx.DiGraph(list(P.keys()))

            ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.
            assert(nx.is_strongly_connected(G) and nx.is_aperiodic(G))

            ### states of P as dict
            states = P.states()

            ### list of labeled states
            S = list(states)

            ### number of state
            n = len(S)

            ### Steady
            Pi = P.steady()
            print('Steady: ',Pi)

            ### find a partition randomly with k classes
            comb = Partition(n)
            comb.AddStateLabels(S)

            ### to store the best kl, the best partition and the lumped matrix
            result = [sys.maxsize, None, None]

            ### list of partitions for k classes
            #print(comb.GetList(k))

            K = range(k,k+1) if k else range(1,n)

            pool = Pool()
            mapped_values = list(tqdm.tqdm(pool.imap_unordered(partial(do_work, n=n, S=S,P=P), K), total=len(K)))

            kl=10000
            result = None
            for d in mapped_values:
                for k,v in d.items():
                    if v[0]<kl:
                        result=v

            for k in K:
                ### number of total partition per 1000 in order to have a certitute of exploring max of partitions for k classes
                N = calculSbyGordon(n,k)*100
                print('\nProcessing with k=%d for %d loop...'%(k,N))
               
                for i in range(N):

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

                        new_states = ['NS%i'%i for i in range(n-1)]
                        init_state = random.choice(S)
                        
                        Q = Lump(partition, Pi, P)
                        p = par[str(partition)]
                        Q_mu = Lifting(Q, Pi, S, p)
                        
                        M = 0.0005
                        #for M in np.linspace(0, 0.1, num=1000):
                        #    tetha = np.ones(n)
                        #    neta = Neta(p, tetha, M, new_states, S)
                        #    print(np.gradient(neta)*kl)    
                        
                        kl = KL(init_state, Q_mu, P, Pi, S)

                        if kl<result[0]:
                            result = [kl, M, list(p), Q]
            
                    print("\r{:2.12%}".format(i/N), end="", flush=True)

            #print("KL divergence: %f"%result[0])
            #print(neta)
            print("\nBest M parameter: %f"%result[1])
            print("Best partition: ",result[2])
            print("Best lumped matrix: ")
            pprint(result[-1])

            # ### test
            print("Steady of the lumped: ",result[-1].steady())
            
            # end time
            end = time.time()

            # total time taken
            print(f"`\nRuntime of the program is {end - start}")
