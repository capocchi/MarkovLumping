#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import islice
import numpy as np
import random

from PartitionsSchemes import *
from nbPartitions import *

def randomStringDigits(stringLength=4):
    """Generate a random string of letters and digits """
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))

def Populate(S,M,C):
    
    Input = iter(S)
    Output = tuple(sorted(tuple(sorted(islice(Input, elem))) for elem in C))

    ### we add a partition only of there are not in the list of no empty cell of M
    indexes_of_empty_cell = np.where(M=='')
    if len(indexes_of_empty_cell[0]) > 0:
        ### the list of not empty cell is M slicing with the empty cell
        list_of_not_empty_cell = M[:indexes_of_empty_cell[0][0]]
        partition = str(Output)
        if partition not in list_of_not_empty_cell:
            M[len(list_of_not_empty_cell)] = partition

def getPartitionsList(M,S,C,nbPartition,verbose=False):
    ### test
    while('' in M):
        ### shuffle the list on place
        np.random.shuffle(S)
        for c in C:
            Populate(S,M,c)

        if verbose: print("\r{:2.12%}".format((len(np.where(M!='')[0])/nbPartition)), end="", flush=True)

def getPartitions(S,C):
    """ S is the states set
        C is the partitions schemes list
    """
    ### shuffle the list on place
    np.random.shuffle(S)

    for c in C:
        Input = iter(S)
        #Output = tuple(sorted(tuple(sorted(islice(Input, elem))) for elem in c))
        Output = tuple(tuple(islice(Input, elem)) for elem in c)
        yield str(Output)

def getPartitionSchemes(n, k):
    """ Return the list of partition schemes
    """
    for c in findCombinations(n): 
        if len(c) == k:
            yield c

class Partition():

    def __init__(self, n):

        self._n = n
        ### list of states
        self._s = np.array(range(self._n))

        self.state_to_label = None
        self.label_to_state = None

    def Count(self,k):
        return calculSbyGordon(self._n,k)

    def AddStateLabels(self,S):
        assert(len(S)==self._n)
        self.state_to_label=dict(zip(self._s,S))
        self.label_to_state=dict(zip(S,self._s))

    def GetStateLabel(self,i):
        assert(isinstance(i,int))
        return self.state_to_label[i]

    def GetStateIndex(self,s):
        assert(isinstance(s,str))
        return self.label_to_state[s]

    def GetClasses(self,k):
        return getPartitionSchemes(self._n,k)

    def GetLabeled(self,k, new_state=None):
        """ Return the random partition for k classes
            The partition is translated with the state labels
        """

        ### AddStateLabels must be called first
        if self.state_to_label:
            partitions_as_string = [eval(a) for a in getPartitions(self._s,self.GetClasses(k))]
            np.random.shuffle(partitions_as_string)

            for partition_as_tuple in partitions_as_string:
                r = []
                for i,c in enumerate(partition_as_tuple):
                    if new_state:
                        for a in c:
                            r.append((*tuple([self.state_to_label[a]]),new_state+str(i)))
                    else:
                        r.append(tuple(self.state_to_label[a] for a in c))

                yield r
        
    def Get(self,k):
        """ Return the random partition for k classes.
        """ 
        ### get partitions depending on the partition schemes C that depends on k!
        return getPartitions(self._s,self.GetClasses(k))

    def GetList(self, k, verbose=False):
        
        nbPartition = self.Count(k)

        ### list of partitions
        _m = np.array(['']*nbPartition)
        #M = M.astype('U256')
        _m = _m.astype('object')

        getPartitionsList(_m, self._s, self.GetClasses(k), nbPartition, verbose)
        
        return _m

if __name__ == '__main__':

    import time, sys, string

    ### $ python PartitionsSchemes 8 2
    ### return the lits of partitioning schemes for k=2 classes in n=8 elements
    if len(sys.argv[1:]) == 2:

        # starting time
        start = time.time()

        ### number of states and desired number of classes
        n,k = map(int,sys.argv[1:])
        
        print(f"Partitioning of {n} states in {k} classes")
        
        ### list of states
        S = np.array(range(n))
        

        ###################################################### List based
        ### the total number possible partition (by M.Gordon)
        #nbPartition = calculSbyGordon(n,k)

        ### list of partitions
        #M = np.array(['']*nbPartition)
        ##M = M.astype('U256')
        #M = M.astype('object')

        #C = list(getPartitionSchemes(n, k))
        #print(f"There is {len(C)} partition scheme(s):",C)

        #print("\nProcessing of possible combination...")
        #getPartitionsList(M,S,C,nbPartition,True)
        #print('\n')
        #print(f"Number of partitions: {nbPartition}->",M)

        ####################################################################

        ### another use with the class Partition!
        P = Partition(n)
        P.AddStateLabels([randomStringDigits() for i in range(n)])
        
        print(f"Ramdom Labeled partitions from Generator (yield):")
        ###
        for a in P.GetLabeled(k):
            print(a)

        #print(f"Ramdom Labeled partitions from Generator (yield) with associated classes:")
        #for a in P.GetLabeled(k, 'NS'):
        #    print(a)

#        print(P.GetList(k))

#        print("Calling the getPartitions random generator you can have:")
#        ### get partitions depending on the partition schemes C that depends on k!
#        for partition in getPartitions(S,C):
#            print(partition)

#        ### another us of getPartitions(S,C)
#        while(True):
#            for partition in getPartitions(S,C):
#                print(partition)

        # end time
        end = time.time()

        # total time taken
        print(f"`\nRuntime of the program is {end - start}")

        #    D = {}
        #    for c in C:
        #        D[str(c)] = 0
        #        for elem in M:
        #            if sorted(map(len,eval(elem))) == sorted(c):
        #                D[str(c)]+=1
        #    print(nbPartition,D)
        #    assert(nbPartition==D)
