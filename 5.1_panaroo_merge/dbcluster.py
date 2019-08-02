import networkx as nx
import igraph
import numpy as np
from joblib import Parallel, delayed
from sklearn.cluster import dbscan
from sklearn import metrics
import matplotlib.pyplot as plt
import sys

G = nx.read_gml(sys.argv[1])
G = nx.convert_node_labels_to_integers(G, label_attribute = "Name")

# use igraph to calculate Jaccard distances quickly
edges = zip(*nx.to_edgelist(G))
G1 = igraph.Graph(len(G), zip(*edges[:2]))
X = 1 - np.array(G1.similarity_jaccard(loops=False)) ## X is an NxN matrix with the pairwise Jaccard distance between every two nodes

# DBSCAN is much faster with metric='precomputed'
db = dbscan(X, metric='precomputed', eps=0.7, min_samples=4)

labels = db[1]
unique_labels = set(labels)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

core_samples_mask = np.zeros_like(labels, dtype=bool)
core_samples_mask[np.argwhere(labels > -1)] = True
n_clusters_ = len(unique_labels) - 1

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = (labels == k)
    print(class_member_mask)
    xy = X[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=14)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
            markeredgecolor='k', markersize=6)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()

label_per_node = {}
nodes = list(G.nodes())
for i in range(0, len(labels)):
	label_per_node[nodes[i]] = labels[i]

nx.set_node_attributes(G, label_per_node, 'dbscan')
## fix the noise -> count how many neighbours from each cluster the noisy ones have and 
## give the correct label
for n in G.nodes(data=True):
	if n[1]["dbscan"] != -1:
		continue
	curr_neighbours = G.neighbors(n[0])
	neighbour_clusters = []
	for k in curr_neighbours:
		neighbour_clusters.append(G.nodes[k]['dbscan'])
	curr_cluster = max(set(neighbour_clusters), key = neighbour_clusters.count) 
	n[1]['dbscan'] = curr_cluster


nx.write_gml(G, path = sys.argv[1] + ".gml")

