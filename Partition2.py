#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import random
import itertools
from tqdm import tqdm

from PartitionsSchemes import *
from nbPartitions import *

### https://codereview.stackexchange.com/questions/1526/finding-all-k-subset-partitions
############################################################
def subsets_k(collection, k): yield from partition_k(collection, k, k)
def partition_k(collection, min, k):
  if len(collection) == 1:
    yield [ collection ]
    return

  first = collection[0]
  for smaller in partition_k(collection[1:], min - 1, k):
    if len(smaller) > k: continue
    # insert `first` in each of the subpartition's subsets
    if len(smaller) >= min:
      for n, subset in enumerate(smaller):
        yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
    # put `first` in its own subset 
    if len(smaller) < k:
        yield [ [ first ] ] + smaller
###############################################################

def randomStringDigits(stringLength=4):
    """Generate a random string of letters and digits """
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))

class Partition():

    def __init__(self, n):

        self._n = n

        ### list of states as number
        self._s = np.array(range(self._n))

        self._state_to_label = None

    def Count(self,k):
        return calculS(self._n,k) #calculSbyGordon(self._n,k)

    def AddStateLabels(self,S):
        assert(len(S)==self._n)
        self._state_to_label = S
        #self.state_to_label=dict(zip(self._s,S))

    def GetStateLabel(self,i):
        #assert(isinstance(i,int))
        return self._state_to_label[i]

    def GetStateLabels(self):
        return self._state_to_label

    def GetLabeled(self,k, new_state=None, splitting_choice=None, coordinate_choice=None):
        """ Return the random partition for k classes
            The partition is translated with the state labels
            new_state: associate the labeled partition with the new partition (usually prefixed by the string NS)
            splitting_choice: specify a splitting choice for the class k (for example if n=4, k=2, splitting_choice can be np.array([1,3]) or np.array([2,2]))
        """

        ### AddStateLabels must be called first
        if self._state_to_label:  
            partitions = self.Get(k)
            
            if coordinate_choice:
                l = list(range(self._n))

            for partition in partitions:
                splitting_cond = int(''.join(str(a) for a in map(len,partition)))==splitting_choice
                coordinate_cond = coordinate_choice and list(itertools.chain.from_iterable(partition))==l
                
                if not splitting_choice or \
                    (splitting_choice and not coordinate_choice and splitting_cond) or \
                    (splitting_choice and coordinate_choice and splitting_cond and coordinate_cond):
                    r = []
                    
                    if new_state and coordinate_choice:
                        iterator = iter(coordinate_choice)
                            
                    for i,c in enumerate(partition):
                        if new_state:
                            for a in c:
                                NS = new_state+str(next(iterator)) if coordinate_choice else new_state+str(i)
                                r.append((*tuple([self.GetStateLabel(a)]),NS))
                        else:
                            r.append(tuple(self.GetStateLabel(a) for a in c))
                        
                    yield r
        
    def Get(self,k):
        """ Return the random partition for k classes.
        """ 
        ### get partitions depending on the partition schemes C that depends on k!
        return subsets_k(list(range(self._n)),k)

if __name__ == '__main__':

    import time, sys, string

    ### $ python Partition2.py 8 2
    ### return the lits of partitioning for k=2 classes in n=8 elements
    if len(sys.argv[1:]) >= 2:

        # starting time
        start = time.time()

        ### number of states and desired number of classes
        n,k = map(int,(sys.argv[1], sys.argv[2]))
        
        print(f"Partitioning of {n} states in {k} classes")
        
        ### list of states
        S = np.array(range(n))
        
        ### another use with the class Partition!
        P = Partition(n)
        P.AddStateLabels([randomStringDigits() for i in range(n)])

        print(f"Ramdom Labeled partitions from Generator (yield)...")

        ### if splitting_choice is expected
        try:
            splitting_choice = int(sys.argv[3])
        except:
            splitting_choice = None
        
        ### if coordinate_choice is expected
        try:
            coordinate_choice = sys.argv[4]
        except:
            coordinate_choice = None

        ###
        cpt=0
        for a in tqdm(P.GetLabeled(k,'NS',splitting_choice=splitting_choice,coordinate_choice=coordinate_choice), total= None if splitting_choice else P.Count(k)):
            cpt+=1
            #print(a)

        ### for n,k the number of partitions must be equal to the gordon formula implemented in the Count method!
        if not splitting_choice and not coordinate_choice:
            assert(P.Count(k)==cpt)

        print("\nOrederd state:",P.state_to_label)

        # end time
        end = time.time()

        # total time taken
        print(f"`\nRuntime of the program is {end - start}")
