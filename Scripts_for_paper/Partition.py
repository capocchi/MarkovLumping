#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import random
from tqdm import tqdm

from more_itertools import locate
 
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

def randomStringDigits(stringLength:int=4)->str:
    """Generate a random string of letters and digits """
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))

class Partition():

    def __init__(self, n:int)->None:

        self._n = n

        ### list of states as number
        self._s = np.array(range(self._n))

        self._state_to_label = None

    def Count(self,k:int)->int:
        return calculS(self._n,k) #calculSbyGordon(self._n,k)

    def AddStateLabels(self,S:list)->None:
        assert(len(S)==self._n)
        self._state_to_label = S
        #self.state_to_label=dict(zip(self._s,S))

    def GetStateLabel(self,i:int)->str:
        #assert(isinstance(i,int))
        return self._state_to_label[i]

    def GetStateLabels(self)->list:
        return self._state_to_label

    def GetNumberedPartition(self,p:list)->str:
        ### p is [('A', 'NS0'), ('B', NS1'),...]
        if set(map(len,p))==set([2]):
            ### d = {'A':'NS0', 'B':'NS1';...}
            d = dict(p)
            ### NS
            new_state = ''.join([n for n in list(ordered(d.values()))[1] if not n.isdigit()])
            return ''.join([d[s][len(new_state):] for s in self.GetStateLabels()])
        
    def GetLabeled(self,k:int, new_state:str=None, coordinate_choice:int=None)->list:
        """ Return the partitions for a classe k
            The partition is translated with the state labels
            new_state: associate the labeled partition with the new partition (usually prefixed by the string NS)
            coordinate_choice: specify a coordinate choice for the class k (for example if n=4, k=2, coordinate_choice can be 1212)
        """

        
        ### AddStateLabels must be called first
        if self._state_to_label:  
            if coordinate_choice:
                
                if new_state:
                    iterator = iter(coordinate_choice) if ',' not in coordinate_choice else iter(coordinate_choice.split(','))
                    r = [(*tuple([a]), new_state+str(next(iterator))) for a in self.GetStateLabels()]
                else:
                    r = [None]*k
                    l = list(map(int,coordinate_choice)) if ',' not in coordinate_choice else list(map(int,coordinate_choice.split(',')))
                    for i in range(k):
                        r[i] = tuple(self.GetStateLabel(j) for j in list(locate(l, lambda a: a == i+1)))
                yield r
            else:
                partitions = self.Get(k)
                
                for partition in partitions:
                    r = []
                    for i,c in enumerate(partition):
                        if new_state:
                            r.extend([(*tuple([self.GetStateLabel(a)]), ''.join([new_state,str(i)])) for a in c])
                        else:
                            r.append(tuple(self.GetStateLabel(a) for a in c))
                    
                    yield r
        
    def Get(self,k:int):
        """ Return the random partition for k classes.
        """ 
        ### get partitions depending on the partition schemes C that depends on k!
        return subsets_k(list(range(self._n)),k)

if __name__ == '__main__':

    import time, sys, string

    ### $ python Partition3.py 8 2
    ### return the list of partitioning for k=2 classes in n=8 elements
    ### $ python Partition3.py 4 2 1122
    ### return the partition for k=2 classes in n=4 elements for the coordinate scheme 1122
    if len(sys.argv[1:]) >= 2:

        # starting time
        start = time.time()

        ### number of states and desired number of classes
        n,k = map(int,(sys.argv[1], sys.argv[2]))
        
        print(f"Partitioning of n={n} states in k={k} classes")
        
        ### list of states
        S = np.array(range(n))
        
        ### another use with the class Partition!
        P = Partition(n)
        P.AddStateLabels([randomStringDigits() for i in range(n)])
 
        ### if coordinate_choice is expected
        try:
            coordinate_choice = sys.argv[3]
        except:
            coordinate_choice = None

        ###
        pbar = tqdm(P.GetLabeled(k,'NS',coordinate_choice=coordinate_choice), desc=f"Processing for k={k}", total=None if coordinate_choice else P.Count(k))
        partitions = []
        for p in pbar:            
            partitions.append(''.join([c[-1][-1] for c in p]))

        ### for n,k the number of partitions must be equal to the gordon formula implemented in the Count method!
        if not coordinate_choice:
            assert(P.Count(k)==len(partitions))

        print("\nOrederd state:",P.GetStateLabels())
        print("Partitions:",partitions)

        print(list(P.GetLabeled(k,coordinate_choice=coordinate_choice)))
        print(list(P.GetLabeled(k,'NS',coordinate_choice=coordinate_choice)))

        # end time
        end = time.time()

        # total time taken
        print(f"`\nRuntime of the program is {end - start}")
