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
### Name: main_loop.py
### Author: L. Capocchi
### Version: 1.0
### Description: script to compute the loop lumping from the matrix given in argv[1]. See the __main__ for the usage. 
### Dependencies: pykov, networkx, numpy, matplotlib
### Python version: 3.9.17
### Date: 18/11/2023
#######################################################################

import os, sys, time
import pykov
import networkx as nx
import statistics

from more_itertools import locate
from sys import platform
from multiprocessing import freeze_support

from Partition import Partition
from Lifting import Lump, KL, Lifting

PLOT = False
WRITE_FILE = False
STAT = False
EXPORT_GRAPH = True

if PLOT:
    import matplotlib.pyplot as plt

def export_graph_to_graphml(P, filename):

    G = nx.Graph()
    
    for k,v in P.items():
        G.add_edge(k[0].replace('/',' \n '), k[1].replace('/',' \n '), weight=float(v))

    nx.write_graphml(G, filename)

def get_edges_with_weights(P:pykov.Chain)->tuple:
    """
    """  
    for s1 in P.states():
        for s2,v in P.mfpt_to(s1).items():
            yield (s1,s2,v)

__getMFTP_cache = {}
def getMFTP(P:pykov.Chain,s):
    if s not in __getMFTP_cache:
        __getMFTP_cache[s] = P.mfpt_to(s).items()
    return __getMFTP_cache[s]

def getMFTPs(P:pykov.Chain,p)->float:
    L = [v for s in p for s2,v in getMFTP(P,s) if s2 not in p]
    
    
#    for s in p:
#        for s2,v in sorted(getMFTP(P,s)):
#            print(s2,p,v)
#            sys.exit()           
    return min([abs(a-b) for a,b in zip(L[::2], L[1::2])])
    
def getMFTPs2(P:pykov.Chain,p)->float:
    for s in p:
        a = None
        b = None
        r = 0.0
        min = 100000000
        for s2,v in getMFTP(P,s):
            if s2 not in p:
                if a:
                    b=v
                else:
                    a=v
                if a and b:
                    r = abs(a-b)
                    if r < min:
                        min = r
                    a = None
                    b = None
    return min

    
def get_ordered_partitions_from_mftp(S:[str],P:pykov.Chain)->tuple:
    """ Get orderded list of partitions.
    """ 
    n = len(S)
    partitionObject = Partition(n)
    partitionObject.AddStateLabels(S)
    
    dd = {}
    
    for c in partitionObject.GetLabeled(k=n-1):
        listOfLengths = map(len,c)
        ### for all couple of states p / for k=n-1 the lenght of the list is 1
        ### consider only the partition with 2 states
        ### the aggregation of 2 states is the best choice... (only for these kind of matrix ?)
        
        for i in locate(listOfLengths, lambda a: a == 2):
            p = c[i]
            mfpt = getMFTPs(P,p)
            if p in dd:
               if mfpt < dd[d][-1]:
                   dd[p]=(c,mfpt)
            else:
                dd[p]=(c,mfpt)
              
    ### heristic is based on the mean of mftp values
    mean = statistics.fmean([v[1] for v in dd.values()])

    for _,v in dd.items():
        if v[1] <= mean:
            yield v[0]

def getMFTPAnalysis(S:[str],P:pykov.Chain)->tuple:
    Pi = P.steady()
    ### result varaible - kl = 1.0 is the max value; we find the min.
    result = {'kl':1000.0,'partition':None, 'Q':None}

    ### loop on partitions to find the best from the mftp analysis
    for p in get_ordered_partitions_from_mftp(S,P):
        #partition = {"/".join(a):a for a in p}
        partition = {''.join(['NS',str(i)]):a for i,a in enumerate(p)}
         
        ### compute the kl divergence rate
        Q = Lump(partition, Pi, P)
        p = [(v, k1) for k1,v1 in partition.items() for v in v1]
        Q_mu = Lifting(Q, Pi, S, p)
        kl = KL(S, P, Pi, Q_mu)
    
        ### store the best kl and partition
        if kl < result['kl']:
            result['kl']=kl
            result['partition']=p
            result['Q']=Q
    
    return result

def getMFTPAnalysis2(S:[str],P:pykov.Chain)->tuple:
    
    Pi = P.steady()

    ### result varaible - kl = 1.0 is the max value; we find the min.
    result = []
   
    ### loop on partitions to find the best from the mftp analysis
    for p in get_ordered_partitions_from_mftp(S,P):
        #partition = {"/".join(a):a for a in p}
        partition = {''.join(['NS',str(i)]):a for i,a in enumerate(p)}
         
        ### compute the kl divergence rate
        Q = Lump(partition, Pi, P)
        pp = [(v, k1) for k1,v1 in partition.items() for v in v1]
        Q_mu = Lifting(Q, Pi, S, pp)
        kl = KL(S, P, Pi, Q_mu)
         
        result.append((kl,pp,Q))
    
    for n in sorted(result,key=lambda x: x[0]):
        yield n

def getMFTPAnalysis3(S:[str],P:pykov.Chain)->tuple:
    
    Pi = P.steady()

    ### result varaible - kl = 1.0 is the max value; we find the min.
    result = []
   
    ### loop on partitions to find the best from the mftp analysis
    for p in get_ordered_partitions_from_mftp(S,P):
        #partition = {"/".join(a):a for a in p}
        partition = {''.join(['NS',str(i)]):a for i,a in enumerate(p)}
         
        ### compute the kl divergence rate
        Q = Lump(partition, Pi, P)
        pp = [(v, k1) for k1,v1 in partition.items() for v in v1]
        # Q_mu = Lifting(Q, Pi, S, pp)
        # kl = KL(S, P, Pi, Q_mu)
         
        # result.append((kl,pp,Q))
        result.append((pp,Q))

    for n in sorted(result,key=lambda x: x[0]):
        yield n

def display_tree(P, root):
  """Displays a tree structure based on the provided transition matrix P and root node.

  Args:
      P (dict): Transition matrix, where keys are state pairs (source, target) and
                 values are transition probabilities.
      root (str): The starting node of the tree.
  """

  # Create a directed acyclic graph (DAG) to ensure tree structure
  G = nx.DiGraph()

  # Breadth-First Search (BFS) traversal to efficiently add nodes and edges
  queue = [root]
  visited = set()

  current_node = queue.pop(0)
  if current_node not in visited:
      visited.add(current_node)
      G.add_node(current_node)

        # Find outgoing edges (children) from the current node
      outgoing_edges = [(source_target[1],w) for source_target, w in P.items() if source_target[0] == current_node]
    #   for source, target in P.items():
        #   print(source,target, current_node)
      for child,w in outgoing_edges:
        G.add_edge(current_node, child, weight=float(w))
        queue.append(child)

  # Separate edges based on transition probability (optional customization)
  elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] > 0.5]
  esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= 0.5]

  # Node positions (consider using a more layout algorithm for complex trees)
  pos = nx.spring_layout(G, seed=7)
  
  # Nodes
  nx.draw_networkx_nodes(G, pos, node_color='#DCDCDC', node_size=700)

  # Edges
  nx.draw_networkx_edges(G, pos, edgelist=elarge, width=3, alpha=1.0)
  nx.draw_networkx_edges(
      G, pos, edgelist=esmall, width=2, alpha=0.7, edge_color="b", style="dashed"
  )

  # Labels
  nx.draw_networkx_labels(G, pos, font_size=10, font_family="sans-serif")

  # edge weight labels
  edge_labels = nx.get_edge_attributes(G, "weight")
  
  for source_target, w in P.items():
      if source_target in edge_labels.keys():
        edge_labels[source_target] = f"prob={str(edge_labels[source_target])}"

  nx.draw_networkx_edge_labels(G, pos, edge_labels)

  # Customize plot appearance (optional)
  ax = plt.gca()
  ax.margins(0.08)
  plt.axis("off")
  plt.tight_layout()
  plt.title(f"Tree Structure (Root: {root})")  # Add a title with root node

  plt.show()

def displayGraph(P):
    """Display the state graph of P
        https://networkx.org/documentation/stable/auto_examples/drawing/plot_weighted_graph.html

    Args:
        P (_type_): Matrix
    """

    G = nx.Graph()
    
    for k,v in P.items():
        G.add_edge(k[0].replace('/',' \n '), k[1].replace('/',' \n '), weight=float(v))
        
    elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] > 0.5]
    esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= 0.5]

    pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility
    
    # nodes
    nx.draw_networkx_nodes(G, pos, node_size=700,  node_color='#DCDCDC')
    
    # edges
    nx.draw_networkx_edges(G, pos, edgelist=elarge, width=3)
    nx.draw_networkx_edges(
        G, pos, edgelist=esmall, width=3, alpha=0.5, edge_color="b", style="dashed"
    )

    # labels
    nx.draw_networkx_labels(G, pos, font_size=10, font_family="sans-serif")
    #edge_labels = nx.get_edge_attributes(G, "weight")
    #nx.draw_networkx_edge_labels(G, pos, edge_labels)

    ax = plt.gca()
    ax.margins(0.08)
    plt.axis("off")
    plt.tight_layout()
    plt.show()

    #if WRITE_FILE:
    #nx.write_graphml_lxml(G, f"{len(P)}.graphml")
        
def trace(n:int, kl:float, p:Partition)->None:
    """Print results

    Args:
        n (int): size of the aggregated current Markov chain
        kl (float): KL
        p (Partition): best partition
    """
    
    print(f"----------------------------------------------------------{n}x{n}")
    print(f"Best KL:{kl}")

    ### transform states in order to avoid losing the fusion of states
    lst = list(map(lambda a: a[-1],p))
    new_state_for_aggregation = max(lst,key=lst.count)
    d = dict(p)
    state_to_aggregate = [k for k,v in d.items() if v == new_state_for_aggregation]
    d = {value: key for key, value in d.items()}
    d[new_state_for_aggregation] = "/".join(state_to_aggregate)

    print(f"Aggregate states: {d[new_state_for_aggregation]}")
    
    return d, d[new_state_for_aggregation]

if __name__ == '__main__':

    # for Windows support of tqdm
    if platform == "win32":
        freeze_support()

    ### filename of the P Matrix
    try:
        fn = sys.argv[1]
    except:
        sys.exit()
    else:
        # starting time
        start1 = time.time()
        
        ### nom de fichier
        filename = os.path.splitext(fn)[0]

        ### Répertoire de base
        base_dir = os.path.dirname(os.path.abspath(__file__))

        ### list of aggregate states during the process
        AS = []

        ### Markov matrix
        #try:
        #    P = pykov.readmat(fn)
        #except AttributeError as info:
        P = pykov.Chain()
        with open(fn,'r') as f:
            for s in f.readlines():
                try:
                    s1,s2,p = s.split(' ')
                    P[(s1,s2)]=float(p.strip())
                except Exception as e:
                    pass
        ### extract the matrix dim from the filename passed as input of the script
        #n = int(fn.split('x')[-1].split('_')[0].split('.dat')[0])
        
        S = tuple(sorted(list(P.states())))
        kl,p,Q = next(getMFTPAnalysis2(S,P))
        n = len(S)
        
        # starting time
        end1 = time.time()

        # total time taken
        print(f"\nRuntime of the first phase of the algorithm BESTA is {end1 - start1}s")

        if PLOT:
            X = []
            K_L = []
        
        if STAT:    
            STEADY = {n:P.steady()}
            if not PLOT:
                K_L = []
                
        if PLOT:
            X.append(n)
            K_L.append(kl)

        ### condition for table 2
        #cond = "n>2"

        ### condition for table3
        kl=new_kl=diff=0.0
        # cond = "new_kl <= kl*(1+0.5) or kl==0.0"
    
        cond = "kl==0.0"

        # N=len(S)/2
        # cond = 'n>N'
        
        ### stopping condition
        while(eval(cond)):
            
            # starting time
            start2 = time.time()

            ### P must be ergotic i.e. the transition matrix must be irreducible and acyclic.            
            G = nx.DiGraph(list(P.keys()), directed=True)
            nx.strongly_connected_components(G)
            assert nx.is_strongly_connected(G) and nx.is_aperiodic(G), f"Matrix is not ergodic!"
            
            ### just for stdout
            d,aggregate_state = trace(n,kl,p)
            AS.append(aggregate_state)
            #ranks_pr = nx.pagerank(G)
            #print(statistics.mean(ranks_pr.values()))

            #print(f"PageRank : {ranks_pr}")
            #nx.draw_circular(G, node_size=400, with_labels=True)
            #plt.show()
            #exit()

            ### P transition matrix of the Markoc chain to lump
            P = pykov.Chain()
            for k,v in Q.items():
                P[(d[k[0]],d[k[1]])] = v
            
            ### set of states
            S = tuple(sorted(P.states()))

            if "n>2" not in cond:
                kl = new_kl
                
            ### mfpt analisys (REDUCED function that call FOP)
            new_kl,p,Q = next(getMFTPAnalysis2(S,P))
            
            ### update variables
            n = len(S)
           
            if "n>2" in cond:
                kl = new_kl

            # print(new_kl, kl*(1+0.5))
            
            # starting time
            end2 = time.time()

            # total time taken
            print(f"\nRuntime of the 'while' phase of the algorithm BESTA for n={n} is {end2 - start2}s")

            if PLOT:
                X.append(n)
                K_L.append(kl)
                displayGraph(dict(P))
                root = next((s for s in S if '000' in s), None)
                display_tree(dict(P),root)
            
            if STAT:
                STEADY[n]=P.steady()
                if not PLOT:
                    K_L.append(kl)
                
            if WRITE_FILE:
                
                # Nom du fichier transformé                
                new_fn = f"{filename}_{n}x{n}.dat"

                # Répertoire cible
                target_dir = os.path.join(base_dir, '..', 'Matrix')
                
                # Chemin complet du fichier
                fn = os.path.join(target_dir, new_fn)

                # fn = os.path.join(os.pardir,'Matrix',f"{fn.split('.dat')[0]}_{n}x{n}.dat")
                if os.path.exists(fn):
                    os.remove(fn)

                with open(fn, 'w') as f:
                    for k,v in dict(P).items():
                        f.write(f"{k[0]} {k[1]} {v} \n")
         
                print(f"Write file in {fn}")

        ### just for stdout
        _, aggregate_state = trace(n,kl,p)
        #ranks_pr = nx.pagerank(G)
        #print(statistics.mean(ranks_pr.values()))
        #print(f"PageRank : {ranks_pr}")
          
        AS.append(aggregate_state)

        # end time
        end3 = time.time()

        # total time taken
        print(f"\nRuntime of the program is {end3 - start1}s")
     
        if STAT:    
            import pandas as pd
            import pprint
            s = pd.Series(K_L)
            print(s.describe())
            pprint.pprint(STEADY)
            
        if PLOT:
            X = list(map(lambda a : f"{a} x {a}", map(str,X)))
            plt.plot(X, K_L, label="KL",marker="o")

            # Ajouter les valeurs à la verticale au-dessus de chaque point
            for i, (x, y) in enumerate(zip(X, K_L)):
                plt.text(i, y+0.001, f"{AS[i]}", ha='center', va='bottom', rotation=90, fontsize=8)

            if STAT:
                plt.plot(X,[s.describe()['std']]*len(X), label='std', linestyle="-")
                
            # show a legend on the plot
            #plt.legend()
            #plt.axis([max(X),min(X),min(K_L),max(K_L)])
            plt.xticks(rotation=45)
            plt.grid()
            
            plt.ylabel("KL")
            
            # function to show the plot
            plt.show()
        
        if EXPORT_GRAPH:
            filename = "graph"
            export_graph_to_graphml(dict(P),f"{filename}_{n}x{n}.graphml")
            print(f"{filename}_{n}x{n}.graphml exported!")

            ### test
            # G = import_graph_from_graphml(f"{filename}_{n}x{n}.graphml")

            # Vous pouvez maintenant utiliser G comme un objet Graph networkx
            # print(G.nodes())  # Affiche les nœuds du graphe
            # print(G.edges())  # Affiche les arêtes du graphe