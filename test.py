 
import pykov
import networkx as nx
import pprint 
import sys


#S = ['S', 'R', 'C']

fn = sys.argv[-1]
P = pykov.readmat(fn)

#print(P.mfpt_to('A'))
#print(P.mfpt_to('B'))
#print(P.mfpt_to('C'))
#print(P.mfpt_to('D'))
#print(P.mfpt_to('E'))
#print(P.mfpt_to('F'))
#print(P.mfpt_to('G'))
#print(P.mfpt_to('F'))

#G = nx.DiGraph(list(P.keys()), directed=True)
#assert(nx.is_strongly_connected(G) and nx.is_aperiodic(G))

edges_with_weights = [ (s1,s2,v) for s1 in P.states() for s2,v in P.mfpt_to(s1).items()]

pprint.pprint(edges_with_weights)

import numpy as np

def get_mfpt(P):
    t = { (s1,s2):v for s1 in P.states() for s2,v in P.mfpt_to(s1).items()}
    r = {}
    for k,v in t.items():
        a = tuple(sorted(k))
        if a not in r.keys():
            r[a]=v
        else:
            r[a]+=v

    return [(k[0],k[1],v) for k,v in r.items()]

pprint.pprint(get_mfpt(P))
#pprint.pprint(sorted(get_mfpt(P), key=lambda tup: tup[-1]))

G = nx.Graph()
G.add_weighted_edges_from(get_mfpt(P))
#G.add_weighted_edges_from(edges_with_weights, weight='weight')

#print(nx.clustering(G))
print(nx.generalized_degree(G))

#import community
#print(community.best_partition(G))

#from networkx.algorithms.community import coverage, is_partition
#print(coverage(G,('S','CR')), is_partition(G,('SRC')))

#pprint.pprint(get_mfpt(P))
#pprint.pprint(len(set([a[-1] for a in get_mfpt(P)])))

d={}
for c in edges_with_weights:
    if tuple(sorted(c[:-1])) not in d:
        d.update({tuple(sorted(c[:-1])):c[-1]})
    else:
        d[tuple(sorted(c[:-1]))] = d[tuple(sorted(c[:-1]))]-c[-1] if d[tuple(sorted(c[:-1]))]>c[-1] else c[-1] - d[tuple(sorted(c[:-1]))]

                    

print("\nGeneralized_degree:")
print(set([b for v in nx.generalized_degree(G).values() for b in v[1]]))

#pprint.pprint(d)


#            from networkx.algorithms import community
#            print(list(community.girvan_newman(G)))


#print(edges_with_weights)

#G = nx.DiGraph()
#G.add_weighted_edges_from(get_mfpt(P), weight='weight')

#print(nx.generalized_degree(G))

#import community

#dendo = community.generate_dendrogram(G, resolution=2)
#for level in range(len(dendo)):
#    print("partition at level", level,"is", community.partition_at_level(dendo, level))

#import collections

#partition = community.best_partition(G, weight='weight')
#print(partition)
#print(community.modularity(partition, G))
#values = [partition.get(node) for node in G.nodes()]
#counter=collections.Counter(values)
#print(counter)

##########################################################


import matplotlib.pyplot as plt
#sp = nx.spring_layout(G)
#nx.draw_networkx(G, pos=sp, with_labels=True)
            # plt.axes('off')
#plt.show()

#pos = nx.spring_layout(G)

#nx.draw(G, pos,with_labels=True)
#nx.draw_networkx_edge_labels(G,pos=pos)
#plt.show()

#preds = nx.jaccard_coefficient(G, G.edges)
#for u, v, j in preds:
#    print(u,v,j)

#print(edges_with_weights)

    #P = pykov.Chain({(S[0],S[0]):0.97,
    #                 (S[0],S[1]):0.01,
    #                 (S[0],S[2]):0.02,
    #                 (S[1],S[0]):0.02,
    #                 (S[1],S[1]):0.48,
    #                 (S[1],S[2]):0.50,
    #                 (S[2],S[0]):0.01,
    #                 (S[2],S[1]):0.75,
    #                 (S[2],S[2]):0.24
    #                 })
