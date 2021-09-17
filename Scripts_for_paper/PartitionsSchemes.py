#!/usr/bin/env python
# -*- coding: utf-8 -*-

### from https://www.geeksforgeeks.org/find-all-combinations-that-adds-upto-given-number-2/

import numpy as np

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
        return 

    # If combination is 
    # found, print it 
    if (reducedNum == 0): 
        r.append(np.sort(arr[:index]))
        return

    # Find the previous number stored in arr[]. 
    # It helps in maintaining increasing order 
    prev = 1 if(index == 0) else arr[index - 1] 

    # note loop starts from previous 
    # number i.e. at array location 
    # index - 1 
    for k in range(prev, num + 1): 
        # next element of array is k 
        arr[index] = k

        # call recursively with 
        # reduced number 
        findCombinationsUtil(arr, index + 1, num, reducedNum - k, r) 

# Function to find out all 
# combinations of positive numbers 
# that add upto given number. 
# It uses findCombinationsUtil() 
def findCombinations(n): 

    # array to store the combinations 
    # It can contain max n elements 
    out = []

    # find all combinations 
    findCombinationsUtil(np.zeros((n,), dtype=int), 0, n, n, out)
    
    return out

if __name__ == '__main__':

    import sys

    ### $ python PartitionsSchemes.py 8 2
    ### return the lits of partitioning schemes for k=2 classes in n=8 elements
    if len(sys.argv[1:]) == 2:
        ### number of states and requiered partitionning 
        n,k = map(int,sys.argv[1:])

        C = [c for c in findCombinations(n) if len(c) == k]

    ### $ python PartitionsSchemes.py 8 
    ### return the list of all partitioning schemes for n=8 elements   
    elif len(sys.argv[1:]) == 1:

        n = int(sys.argv[1])

        ### list of combination
        C = findCombinations(n)
    
    print(C)