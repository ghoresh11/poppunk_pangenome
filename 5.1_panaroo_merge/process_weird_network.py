import sys
import networkx as nx
from numpy import mean, std
import operator
from math import sqrt

print("Reading network...")
G = nx.read_gml(sys.argv[1])
print("Calculating betweenness centrality...")
print(len(G))
k = int(sqrt(len(G)))
betweenness_centrality = nx.betweenness_centrality(G, k = k)
betweenness_centrality = sorted(betweenness_centrality.items(), key=operator.itemgetter(1), reverse = True)
values = [x[1] for x in betweenness_centrality]
min_val = mean(values) + std(values)

remove = []
for item in betweenness_centrality:
    if item[1] > min_val:
        remove.append(item[0])
        continue
    break


print(remove)
print("Removing nodes..")
G.remove_nodes_from(remove)

print("Getting connected_components")
ccs = nx.connected_components(G)
for c in ccs:
	print(len(c))

nx.write_gml(G, path = "test.gml")
