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

import numpy as np
import os,sys
import pykov
import itertools

def Lifting(Q, Pi, S, partition):
    d = dict(partition)
    return {(i,j): Pi[j]*Q[d[i], d[j]]/np.sum([Pi[c[0]] for c in partition if d[j] in c]) for i in S for j in S}

def Lump2(partition, Pi, P):
    """ exact lumpability
    """
    ### new pi vector
    new_pi = pykov.Vector()
    for c in partition:
        state, new_state = c
        if new_state not in new_pi:
            new_pi[new_state] = Pi[state]
        else:
            new_pi[new_state] += Pi[state]

    ### new lumped markov chain
    new_chain = pykov.Chain()
    for pair in itertools.product(new_pi.keys(), repeat=2):
        # print pair
        new_p_ij = 0
        p0,p1 = pair
        L1 = [ c[0] for c in partition if c[1] == p1]
        L2 = [ c[0] for c in partition if c[1] == p0]
        for state in L2:
            sum2 = np.sum([P[(state,s)] for s in L1])
            new_p_ij += Pi[state]*sum2

        new_chain[pair] = new_p_ij/new_pi[p0]

    return new_chain

def Lump(partition, Pi, P):
    """ exact lumpability
    """
    ### new pi vector
    new_pi = {new_state:np.sum([Pi[s] for s in sub_partition]) for new_state, sub_partition in partition.items()}
    return { pair:np.sum([Pi[state]*np.sum([P[(state,s)] for s in partition[pair[1]]]) for state in partition[pair[0]]])/new_pi[pair[0]] for pair in itertools.product(partition.keys(), repeat=2)}

### neta
def Neta(partition, tetha, M, new_states, S):
    L =[]
    d = dict(partition)
    nb_state = len(tetha)
    nb_new_state = len(new_states)

    for t in range(nb_state):
        c = M*tetha[t]
        e = 1+np.exp(c)
        a = [0]*nb_new_state
        index = new_states.index(d[S[t]])
        a[index] = 1
        #print a[0]/e, map(lambda g: (((index)*np.exp(c)/(nb_new_state-1))/e)*g, a[1:]), index, d, t
        r = a[0]/e + np.sum([(np.exp(c)/e)*g for g in a[1:]])

        L.append(r)

    return L

def Tetha(tetha, init_state, partitions, M, new_states, P, Pi, S, N=500):
    l = len(tetha)
    T = []
    for t in range(N):
        state = P.move(init_state)
        init_state = state

        r = np.zeros(l)
        D = {}
        for p in partitions:
        
            Q = Lump2(p, Pi, P)
            Q_mu = Lifting(Q, Pi, S, p)
            neta = Neta(p, tetha, M, new_states, S)
            assert(len(neta)>=0)
            kl = KL2(state, Q_mu, P, Pi, S)
            D[kl] = neta

        Sum = np.zeros(len(S))
        for i in D.values():
            Sum = np.sum([Sum,i], axis=0)
            
        for kl,v in D.items():
            neta = np.divide(v, Sum)
            tmp = np.gradient(neta)*kl
            r = np.sum([r,tmp], axis=0)

        for k in range(l):
            tetha[k] -= 1*r[k]/(t+1)

        T.append(tetha)

    return T

def KL(S, P, Pi, Q_mu):
    ### To ensure R(P k Q) is finite, we require P to be absolutely
    ### continuous w.r.t. Q, i.e. Qij = 0 ⇒ Pij = 0.
    r = 0.0
    for i in S:
        r += Pi[i]*np.sum([P[i,j]*np.log(P[i,j]/Q_mu[i,j]) for j in S if Q_mu[i,j] != 0.0 and P[i,j] != 0.0])
    
    return r

def KL2(i, Q_mu, P, Pi, S):
    ### To ensure R(P k Q) is finite, we require P to be absolutely
    ### continuous w.r.t. Q, i.e. Qij = 0 ⇒ Pij = 0.
    L = filter(lambda a: Q_mu[i,a] != 0.0 and P[i,a] != 0.0, S)
    return np.sum(list(map(lambda j: Pi[j]*P[i,j]*np.log(P[i,j]/Q_mu[i,j]), L)))

def process(partition, M, P, new_states, S, par):
    Q = Lump(partition, Pi, P)
    p = par[str(partition)]
    Q_mu = Lifting(Q, Pi, S, p)
    tetha = np.ones(len(S))
    neta = Neta(p, tetha, M, new_states, S)
    kl = KL(P, Q_mu)
    return np.gradient(neta)*kl

def Analyse(T):
    d = {}
    for c,i in [('S'+str(i),10**i) for i in range(10)]:
        for k,t in T.items():
            a,b = str(i*t).split(".")
            if b[0] != '0' and k not in d.keys():
                d.update({k:c})
    print(d)

if __name__ == '__main__':

    from Partition import Partition

    #######################################################

    S = ['S', 'R', 'C']
    n = len(S)
    k = 2  ### classes
    P = pykov.readmat(os.path.join('Matrix','3x3.dat'))
    #P = pykov.Chain({(S[0],S[0]):0.97,
    #                 (S[0],S[1]):0.01,
    #                 (S[0],S[2]):0.02,
    #                 (S[1],S[0]):0.02,
    #                 (S[1],S[1]):0.48,
    #                 (S[1],S[2]):0.50,
    #                 (S[2],S[0]):0.01,
    #                 (S[2],S[1]):0.75,
    #                 (S[2],S[2]):0.24
    #                 })

    M = .0005

    ########################################################

    Pi = P.steady()
    ### should be = [0.3471, 0.3883, 0.2646].
    print("steady: ", Pi)

    ### states of P
    states = P.states()

    ### find a partition randomly with k classes
    comb = Partition(n)
    comb.AddStateLabels(S)

    result = [sys.maxsize,None, None]

    for i in range(100):

        ### one partition for k classes generated randomly
        p = next(comb.GetLabeled(k,'NS'))
        
        ### list of partitions for k classes
        #print(comb.GetList(k))

        ### one partition for k classes and formated to be computed by Lump and more functions...
        partition = {}
        for a, b in p:
            partition.setdefault(b, []).append(a)

        ### to test
        #partition = {'NS0':"('S')", 'NS1':"('R','C')"}

        par = {}
        par[str(partition)] = [(v, k1) for k1,v1 in partition.items() for v in v1]
        
        #print("Partition",par)

        new_states = ['NS%i'%i for i in range(n-1)]
        init_state = 'S'
        
        Q = Lump(partition, Pi, P)
        p = par[str(partition)]
        Q_mu = Lifting(Q, Pi, S, p)
        tetha = np.ones(n)
        neta = Neta(p, tetha, M, new_states, S)
        kl = KL(S, P,Pi,Q_mu)

        if kl<result[0]:
            result = [kl,list(p), Q_mu]
        #print(np.gradient(neta)*kl)
    
    from pprint import pprint

    print(result[1])
    pprint(result[-1])
    # ##### Test lifting

    # from PartitionsSchemes import *
    # from nbPartitions import *
    # from Partition import getPartitions
    # from Markov import Markov
    # import sys

    # k = 2
    # n = 3

    # ### list of states
    # S = np.array(range(n))

    # ### Markov matrix
    # P = Markov(np.array([ [0.97,0.01,0.02],
    #                 [0.02, 0.48,0.50],
    #                 [0.01,0.75,0.24]]))
    
    # P.SetState(S,['S', 'R', 'C'])

    # Pi = P.LabeledSteady()

    # print("steady: ", Pi)

    # #assert(np.array_equal(Pi, np.array([0.34707904,0.38831615,0.26460481])))
    
    # M = P.GetPartition(k)

    # #print(M)

    # ### choose the partition 0 -> {'NS0':"['S']", 'NS1':"['R','C']"}
    # #partition = P.GetLabeledPartition(0)
    
    # partition = {'NS0':"['S']", 'NS1':"['R','C']"}
    
    # par = {}
    # par[str(partition)] = [(v, k1) for k1,v1 in partition.items() for v in eval(v1)]

    # Q = Lump(partition, Pi, P)
    # p = par[str(partition)]
    # S = P.GetStateAsLabel()
    # Q_mu = Lifting(Q, Pi, S, p)

    # print(Q_mu)