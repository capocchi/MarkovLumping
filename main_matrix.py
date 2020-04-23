
# -*- coding: utf-8 -*-

import numpy as np
import random
import string
import sys

### if True, we test the ergoticity
ERGOTIC_TEST=True
### write the matrix in filename given as parmeter
WRITE_FILE = True
### if True, max values is on the diag, else we shuffle the matrix
MAX_DIAG = False

def randomStringDigits(stringLength=4):
    """Generate a random string of letters and digits """
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))

if __name__ == '__main__':

    ### Here's a method with numpy.identity, starting with a k x k identity matrix, adding a drift term to it, and then normalizing.
    ### Specifying a larger high parameter will get you a smaller ratio of "diagonals to the rest"

    ### $ python main_matrix.py 3 3x3_0.6.dat uniform
    ### store in the file 3x3.dat the markov ergotic matrix using uniform distribution with high=0.6

    ### $ python main_matrix.py 4 4x4_0.1.dat uniform 0.1
    ### store in the file 4x4.dat the markov ergotic matrix using uniform distribution with high=0.1

    ### $ python main_matrix.py 3 3x3.dat beta 1 3
    ### store in the file 3x3_1_3.dat the markov ergotic matrix using beta distribution with a=1 and b=3

    ### $ python main_matrix.py 3 3x3_1_0.5.dat binomial 1 0.5
    ### store in the file 3x3_1_0.5.dat the markov ergotic matrix using binomial distribution with m=1 and p=0.5

    ### $ python main_matrix.py 3 3x3_5.dat weibull 5
    ### store in the file 3x3_5.dat the markov ergotic matrix using weibull distribution with a=5

    ### $ python main_matrix.py 3 3x3_5.dat rayleigh 5
    ### store in the file 3x3_5.dat the markov ergotic matrix using raleigh distribution with a=5

    ### dim of matrix
    try:
        n = int(sys.argv[1])
    except:
        sys.exit()
    else:
        ### define filename to store the matrix
        try:
            fn = sys.argv[2]
        except:
            fn = None
        finally:
            ### define the distribution
            try:
                dist=sys.argv[3]
            except:
                dist='uniform'
            finally:

                if not fn: fn='out.dat'
                
                ### list of labels state
                state_labels = [randomStringDigits() for i in range(n)]

                ### matrix randomlly generated
                ### distribution is uniform but can be changed (https://docs.scipy.org/doc/numpy-1.15.0/reference/routines.random.html)
                if dist=='uniform':
                    ### Specifying a larger high parameter will get you a smaller ratio of "diagonals to the rest"
                    try:
                        high = float(sys.argv[4])
                    except:
                        high = 0.60
                    finally:
                        result = np.identity(n) + np.random.uniform(low=0., high=high, size=(n, n))
                ### https://fr.wikipedia.org/wiki/Loi_b%C3%AAta
                elif dist=='beta':
                    ### Specifying a larger high parameter will get you a smaller ratio of "diagonals to the rest"
                    try:
                        a = float(sys.argv[4])
                        b = float(sys.argv[5])
                    except:
                        a = 0.1
                        b = 0.6
                    finally:
                        result = np.identity(n) + np.random.beta(a=0.1, b=0.3, size=(n, n))
                ### 
                elif dist=='binomial':
                    ### Specifying a larger high parameter will get you a smaller ratio of "diagonals to the rest"
                    try:
                        m = float(sys.argv[4])
                        p = float(sys.argv[5])
                    except:
                        m = 1
                        p = 0.5
                    finally:
                        result = np.identity(n) + np.random.binomial(m, p, size=(n, n))
                elif dist=='weibull':
                    ### Specifying a larger high parameter will get you a smaller ratio of "diagonals to the rest"
                    try:
                        a = float(sys.argv[4])
                    except:
                        a = 5.
                    finally:
                        result = np.identity(n) + np.random.weibull(a, size=(n, n))
                elif dist=='rayleigh':
                    ### Specifying a larger high parameter will get you a smaller ratio of "diagonals to the rest"
                    try:
                        a = float(sys.argv[4])
                    except:
                        a = 5.
                    finally:
                        result = np.identity(n) + np.random.rayleigh(size=(n, n))
                else:
                    print('Distribution is unknown')
                    
                    sys.exit()
                result /= result.sum(axis=1, keepdims=1)

                if not MAX_DIAG:
                    perm = np.arange(result.shape[0])
                    np.random.shuffle(perm)
                    result = result[perm]
                
                print(result.round(2))

                if WRITE_FILE:
                    ### write in fn file
                    with open(fn,'w') as f:
                        for i,s1 in enumerate(state_labels):
                            for j,s2 in enumerate(state_labels):
                                f.write("%s %s %f\n"%(s1,s2,result[i,j]))

                if WRITE_FILE and ERGOTIC_TEST:
                    import pykov
                    import networkx as nx

                    P = pykov.readmat(fn)
                    G = nx.DiGraph(list(P.keys()), directed=True)
                    assert(nx.is_strongly_connected(G) and nx.is_aperiodic(G))
