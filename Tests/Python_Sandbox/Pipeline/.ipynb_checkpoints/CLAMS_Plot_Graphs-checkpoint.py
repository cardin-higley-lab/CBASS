import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import networkx as nx
import community as community_louvain

def PlotGraphs (db2Adj, labels=None):
    fig,ax1=plt.subplots(figsize=(10,10))
    G=nx.from_numpy_matrix(db2Adj)
    partition = community_louvain.best_partition(G, resolution=1,random_state=0)
    pos = nx.spring_layout(G)
    # color the nodes according to their partition
    cmap = cm.get_cmap('viridis', max(labels) + 1)
    
    if labels is None:
        print('using the labels from Louvain')
        labels = (partition.values())
        
    nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=600,
                           cmap=cmap, node_color=list(labels), ax=ax1)
    nx.draw_networkx_edges(G, pos, alpha=0.5, ax=ax1)

    nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif", font_color='r', ax=ax1)
    all_weights = []

    #4 a. Iterate through the graph nodes to gather all the weights
    for (node1,node2,data) in G.edges(data=True):
        all_weights.append(data['weight']) #we'll use this when determining edge thickness

    #4 b. Get unique weights
    unique_weights = list(set(all_weights))

    # Plot the edges
    for weight in unique_weights:
        weighted_edges = [(node1,node2) for (node1,node2,edge_attr) in G.edges(data=True) if edge_attr['weight']==weight]
        width = weight*len(labels)*10/sum(all_weights)
        nx.draw_networkx_edges(G,pos,edgelist=weighted_edges,width=width, ax=ax1)
    plt.show()