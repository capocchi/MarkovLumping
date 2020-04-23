#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import chain, combinations, product, permutations
import time
from multiprocessing import Pool
from numba import jit, cuda 

import sys, os

# starting time
start = time.time()

def subsets(n,r):
    """ Note this only returns non empty subsets of arr"""
    L = range(len(n))[1:]
    r = chain(*[combinations(n, i) for i in L])
    #print(list(r))
    #print(len(list(r)))
    return r

def k_subset(n, r):
    s_arr = sorted(n)
    out = []
    n = subsets(n,r)
    for c in combinations(n, r):
        ### each values of c are in s_arr (not repetition of value)
        if sorted(chain(*c)) == s_arr:
        #item = list(chain(*c))
        #len(c) == r and
        #if len(set(item)) == len(item) and set(item) == set(range(r+1)):
            out.append(c)
    return out

import random

# Python3 program to find out all 
# combinations of positive 
# numbers that add upto given number 

# arr - array to store the combination 
# index - next location in array 
# num - given number 
# reducedNum - reduced number 
def findCombinationsUtil(arr, index, num, reducedNum, r): 

    # Base condition 
    if (reducedNum < 0): 
        return; 

    # If combination is 
    # found, print it 
    if (reducedNum == 0): 
        a = [arr[i] for i in range(index)]
        r.append(sorted(a))
        return

    # Find the previous number stored in arr[]. 
    # It helps in maintaining increasing order 
    prev = 1 if(index == 0) else arr[index - 1]; 

    # note loop starts from previous 
    # number i.e. at array location 
    # index - 1 
    for k in range(prev, num + 1): 
        
        # next element of array is k 
        arr[index] = k; 

        # call recursively with 
        # reduced number 
        findCombinationsUtil(arr, index + 1, num, 
                                reducedNum - k, r); 

# Function to find out all 
# combinations of positive numbers 
# that add upto given number. 
# It uses findCombinationsUtil() 
def findCombinations(n): 

    # array to store the combinations 
    # It can contain max n elements 
    out = []

    # find all combinations 
    findCombinationsUtil([0]*n, 0, n, n, out)
    
    return out

def calculS(n,k):
    #Tester les cas "convention".
    if n==0 and k==0:
        #On traite ici le cas n=k=0
        return 1
    elif n==0 or k==0:
        #On traite ici les cas 0=n<k et 0=k<n
        return 0
    elif k>n:
        #On traite ici le cas qui n’a pas de sens k>n
        return 0
        #cas recursif
    else:
        res = calculS(n-1,k-1) + k*calculS(n-1,k)
    
    return res

def func(S,M,s):
    """  Remark: length_to_split is known because fn is in a loop depending on length_to_split...
    """
    
    Input = iter(S)
    #Output = [sorted(S[:elem]) for elem in length_to_split]
    #Output = [sorted(islice(S, elem)) for elem in length_to_split]
    Output = tuple(sorted(tuple(sorted(islice(Input, elem))) for elem in length_to_split))

    ### we add a partition only of there are not in the list of no empty cell of M
    indexes_of_empty_cell = np.where(M==s)
    if len(indexes_of_empty_cell[0]) > 0:
        ### the list of not empty cell is M slicing with the empty cell
        list_of_not_empty_cell = M[:indexes_of_empty_cell[0][0]]
        partition = str(Output)
        if partition not in list_of_not_empty_cell:
            M[len(list_of_not_empty_cell)] = partition.encode("ascii", "ignore")
    
        #M = np.append(M, str(Output))
        #M = np.insert(M[1:], M.size-1, str(Output))


def slow_function(all_arrays, S, C):
    
    #for M in all_arrays:
    M=all_arrays[0]
    print(M[0])
    ### test
    while('' in M):
        ### shuffle the list on place
        np.random.shuffle(S)
        for length_to_split in C:
            Input = iter(S)
            Output = tuple(sorted(tuple(sorted(islice(Input, elem))) for elem in length_to_split))

            ### we add a partition only of there are not in the list of no empty cell of M
            indexes_of_empty_cell = np.where(M=='')
            if len(indexes_of_empty_cell[0]) > 0:
                ### the list of not empty cell is M slicing with the empty cell
                list_of_not_empty_cell = M[:indexes_of_empty_cell[0][0]]
                partition = str(Output)
                if partition not in list_of_not_empty_cell:
                    M[len(list_of_not_empty_cell)] = partition

        print("\r{:2.2%}".format((len(np.where(M!='')[0])/nbPartition)), end="", flush=True)
    
    return M

    ### save partitions in file
    #path = os.path.join('C:\\','Users', 'Laurent', 'Downloads', 'P', 'partition_%s_%s'%(n,r))
    #print('\nPartitions saved in %s.npy'%path)
    #np.save(path,M)

#import gc
#gc.disable()

if __name__ == '__main__':
    
    import numpy as np
    from itertools import islice
    import math 

    ### number of states and requiered partitionning 
    n,r = map(int,sys.argv[1:])

    ### the total number possible partition (by M.Gordon)
    #nbPartition = calculS(n,r)
    
    nbPartition = int(sum([pow(-1,i)*pow(r-i,n)/(math.factorial(r-i)*math.factorial(i)) for i in range(r)])) if 1<= r < n else 1

    print(f"{n} states \n{r} partionning \nNumber of partitions: {nbPartition}")

    chunks = False
    if nbPartition > sys.maxsize:
        print("Number of states is too large and exceed the max size of int...")
        chunks = True
    else:
        try:
            ### list of partitions
            M = np.array([""]*nbPartition)
            M = M.astype('U256')
        except MemoryError:
            print("Number of states too large and exceed the possible memory allocated to an numpy array...")
            chunks = True

    if chunks:
        
        nbPart=500000
        nbArrays = round(nbPartition/nbPart)+1

        import tempfile
        tempdir = tempfile.gettempdir()
        #M = np.memmap(os.path.join(tempdir,'memmapped.dat'),dtype='U256', mode='w+', shape=(1,nbPartition))
        import h5py
        fn = os.path.join(tempdir,'partition_%s_%s.hdf5'%(n,r))
        os.remove(fn)
        #dt = h5py.special_dtype(vlen=str)
        hdf5_store = h5py.File(fn, "a")
        asciiList = [n.encode("ascii", "ignore") for n in ['']*nbPart]
        hdf5 = hdf5_store.create_dataset("M", (nbPart,), compression="gzip",dtype="S10",data=asciiList)
        M = hdf5
        s = b''
    else:
        s = ''
    
    ### list of states
    S = np.array(range(n))

    ### list of combination
    C = [c for c in findCombinations(n) if len(c) == r]
   
 #   if chunks:
 #       nbPart=50000
 #       nbArrays = round(nbPartition/nbPart)+1
 #       print(f"Trying to create and stack {nbArrays} sub arrays with the size {nbPart}...")
 #       a = np.stack([np.array([""]*nbPart).flat for a in range(nbArrays)])
 #       print("Done!")
 #       array_size_unit = 5000
 #       print(f"Split the array in {array_size_unit} arrays...")
 #       chunks = np.array_split(a, array_size_unit, axis=1)
 #       print("Done!")
 #       print(chunks[0])

        #np.apply_along_axis(slow_function, 0, chunks[0], S, C)
        
 #       sys.exit()

#    print(f"There is {len(C)} possible combination of split")
#    print(C)

#####################################################################

    #arg = (range(n),r)
    #out = k_subset(*arg)
    #r = k_subset(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],n)
    #print(out)

#######################################################################

#    print("List of possible combination...")
    
    ### test
    while(s in M):
    #while(np.size(M) < nbPartition):
    #while(len(M) < len(out)):
        ### shuffle the list on place
        np.random.shuffle(S)
        for length_to_split in C:
            func(S,M,s)

        print("\r{:2.12%}".format((len(np.where(M!=s)[0])/nbPartition)), end="", flush=True)

    if chunks:
        hdf5_store.close()
    else:
        ### save partitions in file
        path = os.path.join('C:\\','Users', 'Laurent', 'Downloads', 'P', 'partition_%s_%s'%(n,r))
        print('\nPartitions saved in %s.npy'%path)
        np.save(path,M)

#    print(f"There is {len(M)} possible combination")
#    print(M)

#    D = {}
#    for c in C:
#        D[str(c)] = 0
#        for elem in M:
#            if sorted(map(len,eval(elem))) == sorted(c):
#                D[str(c)]+=1
#    print(nbPartition,D)
#    assert(nbPartition==D)
############################################################################################################
       
    ### lists are the same ?
    #print(len(M))
    #print(sorted(map(sorted,out)) == sorted(map(sorted,M)))
            
    #from math import factorial

    #combi2 = lambda n,r: int(factorial(r+len(n)-1)/(factorial(r)*factorial(len(n)-1)))
    
    ### C(n,r)
    #combi = lambda n,r: int(factorial(len(n))/(factorial(r)*factorial(len(n)-r)))
    ### number of combinations
    #subset = [combi(n,i) for i in n]
    #print(sum(subset))

    #combi2 = lambda n,r: int(factorial(len(n))/(factorial(r)*factorial(len(n)-r)))
    #k_subset = [range(i) for i in range(r)]
    #print(sum(k_subset))

    
# end time
end = time.time()

# total time taken
print(f"Runtime of the program is {end - start}")

#(py37) C:\Users\Laurent\Dropbox\devsimpy\py3x\Domain\Markov>python Partition.py 3
#Runtime of the program is 1.8270726203918457