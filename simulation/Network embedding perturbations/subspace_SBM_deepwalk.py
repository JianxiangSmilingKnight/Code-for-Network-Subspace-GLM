

import math
import networkx as nx
from node2vec import Node2Vec
import numpy as np
import pandas as pd
import random
# Define the community memberships for each node
# Assign community labels to nodes
np.random.seed(1)
# Define the probabilities within and between communities
# For a 2-block SBM, you'll have a 2x2 probability matrix
# Example: 0.9 within-community probability, 0.1 between-community probability
n = 2000  # Total number of nodes
blocksize=666
p = np.array([[1, 0.3, 0.3],
     [0.3, 1, 0.3],
     [0.3, 0.3, 1]]) # Probability matrix
current_degree = np.sum(p, axis=1)
p0_factor = 2*math.log(n)*3*pow(n,-1)/ np.mean(current_degree)  # Calculate the scaling factor
p0 = p * p0_factor
p1_factor = pow(n,1/2)*3*pow(n,-1)/ np.mean(current_degree)  # Calculate the scaling factor
p1 = p * p1_factor
p2_factor = pow(n,2/3)*3*pow(n,-1)/ np.mean(current_degree)  # Calculate the scaling factor
p2 = p * p2_factor
p3_factor = pow(n,1)/6*3*pow(n,-1)/ np.mean(current_degree)  # Calculate the scaling factor
p3 = p * p3_factor



# Number of different networks to generate
B =100

# Generate B different networks based on the probability matrix and store their embeddings

G = nx.generators.community.stochastic_block_model([blocksize, blocksize+1, blocksize+1], p=p0)

embedding_results = []
for _ in range(B):
    # Generate a stochastic block model graph based on the probability matrix   
    # Perform node2vec embedding
    node2vec = Node2Vec(G, dimensions=3, walk_length=80, num_walks=10)
    model = node2vec.fit(window=10, min_count=1, batch_words=4)
    # Record the learned embeddings in the list
    embeddings = np.array([model.wv[str(node)] for node in G.nodes()])
    embedding_results.append(embeddings)

merged_matrix_2logn = np.concatenate(embedding_results, axis=1)
df = pd.DataFrame(merged_matrix_2logn)  
df.to_csv('80_2logn_embedding_results.csv', index=False)

G = nx.generators.community.stochastic_block_model([blocksize, blocksize+1, blocksize+1], p=p1)

embedding_results = []
for _ in range(B):
    # Generate a stochastic block model graph based on the probability matrix   
    # Perform node2vec embedding
    node2vec = Node2Vec(G, dimensions=3, walk_length=80, num_walks=10)
    model = node2vec.fit(window=10, min_count=1, batch_words=4)
    # Record the learned embeddings in the list
    embeddings = np.array([model.wv[str(node)] for node in G.nodes()])
    embedding_results.append(embeddings)

merged_matrix_n12 = np.concatenate(embedding_results, axis=1)
df = pd.DataFrame(merged_matrix_n12)  
df.to_csv('80_n12_embedding_results.csv', index=False)


G = nx.generators.community.stochastic_block_model([blocksize, blocksize+1, blocksize+1], p=p2)
embedding_results = []
for _ in range(B):
    # Generate a stochastic block model graph based on the probability matrix   
    # Perform node2vec embedding
    node2vec = Node2Vec(G, dimensions=3, walk_length=80, num_walks=10)
    model = node2vec.fit(window=10, min_count=1, batch_words=4)
    # Record the learned embeddings in the list
    embeddings = np.array([model.wv[str(node)] for node in G.nodes()])
    embedding_results.append(embeddings)

merged_matrix_n23 = np.concatenate(embedding_results, axis=1)
df = pd.DataFrame(merged_matrix_n23)  
df.to_csv('80_n23_embedding_results.csv', index=False)
