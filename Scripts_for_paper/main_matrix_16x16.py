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
### Name: main_matrix.py
### Author: L. Capocchi
### Version: 2.0
### Description: script to compute matrix used for the benchmarking.
### Dependencies: numpy, pykov, networkx
### Python version: 3.9
### Date: 09/22/2021
#######################################################################

from re import A
import sys
import pykov
import networkx as nx
 
if __name__ == '__main__':

    ### $ python main_matrix_16x16.py 0.001 0.001 0.001 0.001 16x16_a.dat
    ### store in the file 16x16_a.dat the markov ergotic matrix with alpha=0.001, betha=0.001, gamma=0.001 tetha=0.001

    ### dim of matrix
    
    ### define filename to store the matrix
    try:
        fn = sys.argv[5]
    except:
        fn = None
    finally:
        ### define the distribution
        try:
            a=float(sys.argv[1])
        except:
            a=0.001
        finally:
            ### define the distribution
            try:
                b=float(sys.argv[2])
            except:
                b=0.001
            finally:
                ### define the distribution
                try:
                    c=float(sys.argv[3])
                except:
                    c=0.001
                finally:
                    ### define the distribution
                    try:
                        d=float(sys.argv[4])
                    except:
                        d=0.001
                    finally:
                        if not fn: fn='out.dat'
                        
                        P = pykov.Chain()
                        
                        S = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L','M', 'N', 'O', 'P')
                        
                        delta_a = a/15.0
                        ab=(a+b)/14.0
                        abc=(a+b+c)/13.0
                        abcd=(a+b+c+d)/12.0
                        
                        
                        
                        ### line 1 
                        P[('A','A')] = 0.6-a
                        P[('A','B')] = 0.4-b
                        for s in S[2:]:
                            P[('A',s)] = ab
                        
                        ### line 2
                        P[('B','A')] = abc
                        P[('B','B')] = 0.3-a
                        P[('B','C')] = 0.3-b
                        P[('B','D')] = abc
                        P[('B','E')] = 0.4-c
                        for s in S[5:]:
                            P[('B',s)] = abc
                            
                        ### line 3
                        P[('C','A')] = 0.2-a
                        P[('C','B')] = abcd
                        P[('C','C')] = 0.4-b
                        P[('C','D')] = abcd
                        P[('C','E')] = abcd
                        P[('C','F')] = abcd
                        P[('C','G')] = abcd
                        P[('C','H')] = 0.3-c
                        P[('C','I')] = 0.1-d
                        for s in S[9:]:
                            P[('C',s)] = abcd
                            
                        ### line 4
                        P[('D','A')] = abc
                        P[('D','B')] = abc
                        P[('D','C')] = abc
                        P[('D','D')] = abc
                        P[('D','E')] = abc
                        P[('D','F')] = 0.3-a
                        P[('D','G')] = 0.5-b
                        P[('D','H')] = 0.2-c
                        for s in S[8:]:
                            P[('D',s)] = abc
                            
                        ### line 5
                        P[('E','A')] = ab
                        P[('E','B')] = ab
                        P[('E','C')] = ab
                        P[('E','D')] = ab
                        P[('E','E')] = ab
                        P[('E','F')] = 0.3-a
                        P[('E','G')] = 0.7-b
                        for s in S[7:]:
                            P[('E',s)] = ab
                            
                        ### line 6
                        P[('F','A')] = abc
                        P[('F','B')] = abc
                        P[('F','C')] = abc
                        P[('F','D')] = 0.2-a
                        P[('F','E')] = 0.5-b
                        P[('F','F')] = abc
                        P[('F','G')] = abc
                        P[('F','H')] = 0.3-c
                        for s in S[8:]:
                            P[('F',s)] = abc
                            
                        ### line 7
                        P[('G','A')] = ab
                        P[('G','B')] = ab
                        P[('G','C')] = ab
                        P[('G','D')] = 0.2-a
                        P[('G','E')] = 0.8-b
                        for s in S[5:]:
                            P[('G',s)] = ab
                            
                        ### line 8
                        P[('H','A')] = ab
                        P[('H','B')] = ab
                        P[('H','C')] = ab
                        P[('H','D')] = ab
                        P[('H','E')] = ab
                        P[('H','F')] = ab
                        P[('H','G')] = ab
                        P[('H','H')] = ab
                        P[('H','I')] = 0.5-a
                        P[('H','J')] = ab
                        P[('H','K')] = ab
                        P[('H','L')] = ab
                        P[('H','M')] = 0.5-b
                        for s in S[13:]:
                            P[('H',s)] = ab
                            
                        ### line 9
                        P[('I','A')] = ab
                        P[('I','B')] = ab
                        P[('I','C')] = ab
                        P[('I','D')] = ab
                        P[('I','E')] = ab
                        P[('I','F')] = ab
                        P[('I','G')] = ab
                        P[('I','H')] = ab
                        P[('I','I')] = ab
                        P[('I','J')] = ab
                        P[('I','K')] = 0.4-a
                        P[('I','L')] = 0.6-b
                        for s in S[12:]:
                            P[('I',s)] = ab
                            
                        ### line 10
                        P[('J','A')] = ab
                        P[('J','B')] = ab
                        P[('J','C')] = ab
                        P[('J','D')] = ab
                        P[('J','E')] = ab
                        P[('J','F')] = ab
                        P[('J','G')] = ab
                        P[('J','H')] = ab
                        P[('J','I')] = ab
                        P[('J','J')] = ab
                        P[('J','K')] = 0.5-a
                        P[('J','L')] = 0.5-b
                        for s in S[12:]:
                            P[('J',s)] = ab
                            
                        ### line 11
                        P[('K','A')] = delta_a
                        P[('K','B')] = delta_a
                        P[('K','C')] = delta_a
                        P[('K','D')] = delta_a
                        P[('K','E')] = delta_a
                        P[('K','F')] = delta_a
                        P[('K','G')] = delta_a
                        P[('K','H')] = delta_a
                        P[('K','I')] = delta_a
                        P[('K','J')] = delta_a
                        P[('K','K')] = delta_a
                        P[('K','L')] = delta_a
                        P[('K','M')] = 1-a
                        P[('K','N')] = delta_a
                        P[('K','O')] = delta_a
                        P[('K','P')] = delta_a
                            
                        ### line 12
                        P[('L','A')] = delta_a
                        P[('L','B')] = delta_a
                        P[('L','C')] = delta_a
                        P[('L','D')] = delta_a
                        P[('L','E')] = delta_a
                        P[('L','F')] = delta_a
                        P[('L','G')] = delta_a
                        P[('L','H')] = delta_a
                        P[('L','I')] = delta_a
                        P[('L','J')] = delta_a
                        P[('L','K')] = delta_a
                        P[('L','L')] = delta_a
                        P[('L','M')] = 1-a
                        P[('L','N')] = delta_a
                        P[('L','O')] = delta_a
                        P[('L','P')] = delta_a
                        
                        ### line 13
                        P[('M','A')] = ab
                        P[('M','B')] = ab
                        P[('M','C')] = ab
                        P[('M','D')] = ab
                        P[('M','E')] = ab
                        P[('M','F')] = ab
                        P[('M','G')] = ab
                        P[('M','H')] = ab
                        P[('M','I')] = 0.4-a
                        P[('M','J')] = 0.6-b
                        P[('M','K')] = ab
                        P[('M','L')] = ab
                        P[('M','M')] = ab
                        P[('M','N')] = ab
                        P[('M','O')] = ab
                        P[('M','P')] = ab
                        
                        ### line 14
                        P[('N','A')] = ab
                        P[('N','B')] = ab
                        P[('N','C')] = ab
                        P[('N','D')] = ab
                        P[('N','E')] = ab
                        P[('N','F')] = ab
                        P[('N','G')] = ab
                        P[('N','H')] = ab
                        P[('N','I')] = ab
                        P[('N','J')] = ab
                        P[('N','K')] = ab
                        P[('N','L')] = ab
                        P[('N','M')] = ab
                        P[('N','N')] = 0.4-a
                        P[('N','O')] = 0.6-b
                        P[('N','P')] = ab
                        
                        ### line 15
                        P[('O','A')] = ab
                        P[('O','B')] = ab
                        P[('O','C')] = ab
                        P[('O','D')] = ab
                        P[('O','E')] = ab
                        P[('O','F')] = ab
                        P[('O','G')] = ab
                        P[('O','H')] = ab
                        P[('O','I')] = ab
                        P[('O','J')] = ab
                        P[('O','K')] = ab
                        P[('O','L')] = ab
                        P[('O','M')] = ab
                        P[('O','N')] = 0.8-a
                        P[('O','O')] = ab
                        P[('O','P')] = 0.2-b
                        
                        ### line 16
                        P[('P','A')] = ab
                        P[('P','B')] = ab
                        P[('P','C')] = ab
                        P[('P','D')] = ab
                        P[('P','E')] = ab
                        P[('P','F')] = ab
                        P[('P','G')] = ab
                        P[('P','H')] = ab
                        P[('P','I')] = ab
                        P[('P','J')] = ab
                        P[('P','K')] = ab
                        P[('P','L')] = ab
                        P[('P','M')] = ab
                        P[('P','N')] = 0.5-a
                        P[('P','O')] = ab
                        P[('P','P')] = 0.5-b
                        
                        G = nx.DiGraph(list(P.keys()), directed=True)
                        assert(nx.is_strongly_connected(G) and nx.is_aperiodic(G))
                    
                        ### write in fn file
                        with open(fn,'w') as f:
                            for s1 in S:
                                for s2 in S:
                                    f.write("%s %s %f\n"%(s1,s2,P[s1,s2]))