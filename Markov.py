# -*- coding: utf-8 -*-

import numpy as np

from nbPartitions import calculSbyGordon
from Partition import getPartitionsList, findCombinations

class Markov():

    def __init__(self,p):
        self._p = p

        self._states_as_number = np.array(range(len(self._p)))
        self._m = None

        self._state_to_label = None
        self._steady = None
        self._nb_partitions = np.zeros(len(self._states_as_number), dtype=int)
        self._combinations = [c for c in findCombinations(len(self._states_as_number))]

    def Steady(self):
        if self._steady:
            import operator
            a=np.append(np.transpose(self._p)-np.identity(3),[[1,1,1]],axis=0)
            b=np.transpose(np.array([0,0,0,1]))
            c=np.transpose(a)
            self._steady = np.linalg.solve(c.dot(a), c.dot(b))
            d = [x[1] for x in sorted(self._state_to_label.items(), key=operator.itemgetter(1),reverse=True)]
            self._labeled_steady = dict(zip(d, self._steady))
        return self._steady

    def LabeledSteady(self):
        if not self._steady:
            self._steady = True
            self.Steady()
        return self._labeled_steady

    def SetState(self,states=[],labels=None):
        
        if labels:
            ### association to be more readable (optional)
            self._state_to_label = dict(zip(states,labels))
            self._states_as_label = labels

        self._states_as_number = states

    def GetStateAsNumber(self):
        return self._states_as_number

    def GetStateAsLabel(self):
        if self._states_as_label:
            return self._states_as_label

    def GetLabels(self):
        return self._state_to_label

    def GetProbabilityRow(self,r):
        for (k, v) in self.GetLabels().items():
            if v == r:
                return self._p[k]

    def GetProbabilityCol(self,c):
        for (k, v) in self.GetLabels().items():
            if v == c:
                return self._p[k]

    def GetProbabiltyValue(self,r,c):
        a,b=None,None
        for (k, v) in self.GetLabels().items():
            if v == r:
                a = k
                if b:break
            elif v == c:
                b = k
                if a:break

        return self._p[a][b]

    def GetPartition(self,k):

        ### the total number possible partition (by M.Gordon) 
        if self._nb_partitions[k] == 0.0:
            self._nb_partitions[k] = calculSbyGordon(len(self._states_as_number),k)            

        ### list of partitions
        self._m = np.array(['']*self._nb_partitions[k])
        self._m = self._m.astype('object')
    
        ### change self._m in place
        getPartitionsList(self._m, self._states_as_number, self.GetCombination(k), self._nb_partitions[k])

        if self._state_to_label:
            ### 
            self._labeled_partition = {}
            ### prepare the labeled partition
            for i,partition in enumerate(self._m):
                ### 0, ((0, 2), (1,))
                self._labeled_partition[i] = {}
                for j,p in enumerate(eval(partition)):
                    ### 0, (0, 2)
                    self._labeled_partition[i][''.join(['NS',str(j)])]=str([self._state_to_label[k] for k in p])
            
        return self._m

    def GetLabeledPartition(self,k=None):
        """ Return the labeled partition depending on k
            for example k=0 gives {'NS0': "['S', 'R']", 'NS1': "['C']"}
        """
        if not self._state_to_label:
            self._state_to_label = True
            self.GetPartition(k)
        return self._labeled_partition[k]

    def GetCombination(self, k):
        return [c for c in self.GetCombinations() if len(c) == k]

    def GetCombinations(self):
        return self._combinations