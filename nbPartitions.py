#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, math

def calculS(n,k):
    """ from https://www.prepamag.fr/concours/pdf/corriges.pdf.extraits/2017/PC_MATHS_CENTRALE_1_2017.extrait.pdf
    """
    if n==k:
        return 1
    elif (n==0 or k==0) or k>n:
        return 0
    else:
        return calculS(n-1,k-1) + k*calculS(n-1,k)
    
def calculSbyGordon(n,k):
    return int(sum([pow(-1,i)*pow(k-i,n)/(math.factorial(k-i)*math.factorial(i)) for i in range(k)])) if 1<= k < n else 1

if __name__ == '__main__':

    ### $ python nbPartitions 8 2
    ### return the number of k=2 classes in n=8 elements
    if len(sys.argv[1:]) == 2:
        ### number of states and requiered partitionning 
        n,k = map(int,sys.argv[1:])

        ### the total number possible partition (by M.Gordon)
        #nbPartition = calculS(n,k)
        #print(nbPartition)
        
        ### from P. Gordon
        print('{:.2e}'.format(calculSbyGordon(n,k)))

    ### $ python nbPartitions 8 
    ### return the list of numbers of k=1,2,3,...7 classes in n=8 elements   
    elif len(sys.argv[1:]) == 1:

        import multiprocessing as mp
        
        n = int(sys.argv[1])

        with mp.Pool(mp.cpu_count()) as pool:
            r = pool.starmap(calculSbyGordon, [(n,k) for k in range(n)])

        #print(['{:.2e}'.format(x) for x in r])

        import matplotlib.pyplot as plt
        r = list(map(float,r))
        plt.plot(r)
        plt.ylabel('Number of partitions')
        plt.xlabel('Number of classes k')
        plt.axis([0, n, 0, max(r)])
        plt.grid(True)
        plt.show()