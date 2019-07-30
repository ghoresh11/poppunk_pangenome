import sys
import networkx as nx
from numpy import mean, std
import operator


print("Reading network...")
G = nx.read_gml(sys.argv[1])
print("Calculating betweenness centrality...")
betweenness_centrality = nx.betweenness_centrality(G)
betweenness_centrality = sorted(betweenness_centrality.items(), key=operator.itemgetter(1), reverse = True)
values = [x[1] for x in betweenness_centrality]
min = mean(values) + std(values)

remove = []
for item in betweenness_centrality:
    print(item[1])
    if item[1] > min:
        remove.append(item[0])
        continue
    break


print(remove)
print("Removing nodes..")
G.remove_nodes_from(remove)

print("Getting connected_components")
ccs = nx.connected_components(G)
print(len(list(ccs)))

nx.write_gml(G, path = "test.gml")
