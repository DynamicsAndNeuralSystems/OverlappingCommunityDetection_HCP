import scipy
import pandas as pd
import numpy as np
from scipy import io
import matplotlib.pyplot as plt
import networkx as nx
from fastconsensus.core import construct_sparse_consensus_matrix,threshold_matrix,triadic_closure,check_convergence,get_algorithm
import igraph as ig

# Load in connectivity data
RH_connectivity_data = pd.DataFrame(scipy.io.loadmat('../data/RH.mat')['RH'])

# Define Louvain + consensus clustering parameters
n_partitions=100
threshold=0.2
algorithm='louvain'
max_iterations = 100
convergence_threshold = 0.05
community_algorithm = get_algorithm(algorithm)

# Use a gamma value of 1.3 for Louvain
gamma_val = 1.3
# Convert RH_connectivity_data to ig.Graph object
Glasser_RH_graph = ig.Graph.Weighted_Adjacency(RH_connectivity_data.values.tolist(), mode="undirected")
consensus_matrix = Glasser_RH_graph.copy()

# Initialize weights if they don't exist
if 'weight' not in consensus_matrix.edge_attributes():
    print("Initializing edge weights to 1.0")
    consensus_matrix.es['weight'] = 1.0

# For each iteration, apply Louvain algorithm and update consensus matrix
for iteration in range(max_iterations):
        print(iteration, end = '\r')
        max_triads = consensus_matrix.ecount()

        # Generate n_partitions Louvain clustering solutions
        partitions = []
        for i in range(n_partitions):
            partition = Glasser_RH_graph.community_multilevel(weights='weight', return_levels=True, resolution=gamma_val)[0]
            partitions.append({v: partition.membership[v] for v in range(Glasser_RH_graph.vcount())})
        
        # Construct sparse consensus matrix
        consensus_matrix = construct_sparse_consensus_matrix(consensus_matrix, partitions)
        
        # Thresholding
        consensus_matrix = threshold_matrix(consensus_matrix, threshold)

        # Triadic closure
        consensus_matrix = triadic_closure(consensus_matrix, max_triads, partitions)
        
        # Check for convergence
        if check_convergence(consensus_matrix, convergence_threshold):
            break

# Get final partition based on consensus matrix
final_partition = consensus_matrix.community_multilevel(weights='weight', return_levels=True, resolution=gamma_val)[0]
final_partition = {v: final_partition.membership[v] for v in range(Glasser_RH_graph.vcount())}

# Convert final_partition to a pandas DataFrame
louvain_consensus_modules = pd.DataFrame(final_partition.items(), columns=["Glasser_RH_ROI", "module"])

# Write to CSV
louvain_consensus_modules.to_csv('../data/Louvain_Glasser_RH_SC_consensus_modules_gamma1.3.csv', index=False)
