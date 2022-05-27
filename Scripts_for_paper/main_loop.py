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
### Name: main_loop.py
### Author: L. Capocchi
### Version: 1.0
### Description: script to compute the loop lumping from the matric given in argv[1]. See the __main__ for the usage. 
### Dependencies: pykov, networkx, numpy, matplotlib
### Python version: 3.9
### Date: 11/13/2021
#######################################################################

import os
import pykov, sys, time
import networkx as nx
import statistics

from more_itertools import locate
from sys import platform
from multiprocessing import freeze_support
import multiprocessing as mp

from Partition import Partition
from Lifting import Lump, KL, Lifting
#from nbPartitions import calculSbyGordon

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np 

PLOT = False
WRITE_FILE = False
STAT = True

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
    
    
#    for s in p:
#        for s2,v in sorted(getMFTP(P,s)):
#            print(s2,p,v)
#            sys.exit()           
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
        ### the aggregation of 2 states is the best choice... (only for these kind of matrix ?)
        
        for i in locate(listOfLengths, lambda a: a == 2):
            p = c[i]
            mfpt = getMFTPs(P,p)
            if p in dd:
               if mfpt < dd[d][-1]:
                   dd[p]=(c,mfpt)
            else:
                dd[p]=(c,mfpt)
              
    ### heristic is based on the mean of mftp values
    mean = statistics.fmean([v[1] for v in dd.values()])

    for k,v in dd.items():
        if v[1] <= mean:
            yield v[0]

def getMFTPAnalysis(S:[str],P:pykov.Chain)->tuple:
    Pi = P.steady()
    ### result varaible - kl = 1.0 is the max value; we find the min.
    result = {'kl':1000.0,'partition':None, 'Q':None}

    ### loop on partitions to find the best from the mftp analysis
    for p in get_ordered_partitions_from_mftp(S,P):
        #partition = {"/".join(a):a for a in p}
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
            result['Q']=Q
    
    return result

def getMFTPAnalysis2(S:[str],P:pykov.Chain)->tuple:
    
    Pi = P.steady()
    
    #Pi = pykov.Vector()
    #for a,b in P.steady().items():
    #    Pi[a] = round(b,11)
    
    ### result varaible - kl = 1.0 is the max value; we find the min.
    result = []

   
    ### loop on partitions to find the best from the mftp analysis
    for p in get_ordered_partitions_from_mftp(S,P):
        #partition = {"/".join(a):a for a in p}
        partition = {''.join(['NS',str(i)]):a for i,a in enumerate(p)}
         
        ### compute the kl divergence rate
        Q = Lump(partition, Pi, P)
        pp = [(v, k1) for k1,v1 in partition.items() for v in v1]
        Q_mu = Lifting(Q, Pi, S, pp)
        kl = KL(S, P, Pi, Q_mu)
         
        result.append((kl,pp,Q))
    
    for n in sorted(result,key=lambda x: x[0]):
        yield n
        
def displayGraph(P):
    G = nx.Graph()
   
    for k,v in P.items():
        G.add_edge(k[0], k[1], weight=float(v))
        
    elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] > 0.5]
    esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= 0.5]

    pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility

    # nodes
    nx.draw_networkx_nodes(G, pos, node_size=700)

    # edges
    nx.draw_networkx_edges(G, pos, edgelist=elarge, width=3)
    nx.draw_networkx_edges(
        G, pos, edgelist=esmall, width=3, alpha=0.5, edge_color="b", style="dashed"
    )

    # labels
    nx.draw_networkx_labels(G, pos, font_size=14, font_family="sans-serif")

    ax = plt.gca()
    ax.margins(0.08)
    plt.axis("off")
    plt.tight_layout()
    plt.show()

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
        #try:
        #    P = pykov.readmat(fn)
        #except AttributeError as info:
        P = pykov.Chain()
        with open(fn,'r') as f:
            for s in f.readlines():
                try:
                    s1,s2,p = s.split(' ')
                    P[(s1,s2)]=float(p.strip())
                except Exception as e:
                    pass
        ### extract the matrix dim from the filename passed as input of the script
        #n = int(fn.split('x')[-1].split('_')[0].split('.dat')[0])
        
        S = tuple(sorted(list(P.states())))
        
        #SS = S
        #PP = P
         
        n = len(S)
        
        if PLOT:
            X = []
            K_L = []
        
        if STAT:    
            STEADY = {n:P.steady()}
            #print(STEADY)
            if not PLOT:
                K_L=[]
        
        
        #e2p, p2e = P._el2pos_()
        #p = np.arange(len(e2p)*len(p2e),dtype=float).reshape(len(e2p),len(p2e))
        
        #for k,v in P._dok_(e2p).items():
        #    p[k[0],k[1]]=v
        #    print(k,v)
        
        a = getMFTPAnalysis2(S,P)
        kl,p,Q = next(a)
        
        d_old=10000000
        
        while (n>=2) :
            
            ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
            G = nx.DiGraph(list(P.keys()), directed=True)
            nx.strongly_connected_components(G)
            assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergotic!"
            
            ### Mean First Passage Times Analysis -----------------------------------
            #result = getMFTPAnalysis(S,P)    
            ### ----------------------------------------------------------------------
            
            print(f"----------------------------------------------------------{n}x{n}")
            print(f"Best KL:{kl}")
            
            ### transform states in order to avoid losing the fusion of states
            #p = result['partition']
            lst = list(map(lambda a: a[-1],p))
            new_state_for_aggregation = max(lst,key=lst.count)
            d = dict(p)
            state_to_aggregate = [k for k,v in d.items() if v == new_state_for_aggregation]
            d = {value: key for key, value in d.items()}
            d[new_state_for_aggregation]= "/".join(state_to_aggregate)
            
            
            print(f"Aggregate states: {'/'.join(state_to_aggregate)}")
            
            if n>=3:
                P = pykov.Chain()
                for k,v in Q.items():
                    P[(d[k[0]],d[k[1]])] = v
                
                S = tuple(sorted(P.states()))
                    
                a = getMFTPAnalysis2(S,P)     
                new_kl,p,Q = next(a)
                
                d = abs(new_kl-kl)
                #if d > 10*d_old:
                    ### exit !
                #    n = 2
                #else:
                    # Number of states
                n = len(S)
                kl = new_kl
                d_old = d
            else:
                ### exit !
                n=1
            
            if PLOT:
                X.append(n)
                K_L.append(kl)
                #   displayGraph(dict(P))
            if STAT:
                STEADY[n]=P.steady()
                if not PLOT:
                    K_L.append(kl)
                
            if WRITE_FILE:
                fn = os.path.join(os.pardir,'Matrix',f"{n}x{n}.dat")
                if os.path.exists(fn):
                    os.remove(fn)
                f = open(fn,'w')
                for k,v in dict(P).items():
                    f.write(f"{k[0]} {k[1]} {v} \n")
                f.close()
         
        # end time
        end = time.time()

        # total time taken
        print(f"\nRuntime of the program is {end - start}s")
     
        if PLOT:
            plt.plot(X, K_L)

            # show a legend on the plot
            #plt.legend()
            plt.axis([max(X),min(X),min(K_L),max(K_L)])
            plt.grid()
            
            plt.ylabel("KL")
            plt.xlabel("dim")

            # function to show the plot
            plt.show()
            
            #import statistics
            #st_dev = statistics.stdev(Y)
            #var  = statistics.variance(Y)
            #mean_h = statistic.median_high(Y)
            #print("Standard deviation of the KL: " + str(st_dev))
            #print("Variance of the KL: " + str(var))
        
        if STAT:    
            import pandas as pd
            import pprint
            s = pd.Series(K_L)
            print(s.describe())
            pprint.pprint(STEADY)
