# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 22:07:32 2024

@author: User
"""

import math
import networkx as nx
from node2vec import Node2Vec
import numpy as np
import pandas as pd
import random
# Define the community memberships for each node
# Assign community labels to nodes

# Define the probabilities within and between communities
# For a 2-block SBM, you'll have a 2x2 probability matrix
# Example: 0.9 within-community probability, 0.1 between-community probability
np.random.seed(1)

N = 4000  # network size
U = 0.8 * np.array([range(1, N+1)]) / (N + 1)
V = 0.8 * np.array([range(1, N+1)]) / (N + 1)
UU = np.transpose(U) @ np.ones((1, N))
VV = np.ones((N, 1)) @ V
W = 15 * (np.abs(UU - VV))**(0.8) - 0.1
W = 1 - 1 / (1 + 0.75 * np.exp(-W))
neighborhood_matrix = W
# Generate a permutation vector and permute the rows and columns of the neighborhood matrix
permutation_vector = np.random.permutation(N)
big_P = W[permutation_vector[:, None], permutation_vector]

# Further operations
n = 2000
P = big_P[0:n, 0:n]
p0 = P / np.sum(P) * 2 * np.log(n) * n
p0=np.array(p0)
p1 = p0 * (n**(1/2) / 2 / np.log(n))
p2 = p0 * (n**(2/3) / 2 / np.log(n))
p3 = p0 * (n / 6 / 2 / np.log(n))

p0_net = np.zeros_like(p0, dtype=int)

# Generate binary data based on the probabilities in the input matrix
for i in range(n):
    for j in range(i+1,n):
        p0_net[i][j] = np.random.choice([0, 1], p=[1 - p0[i][j], p0[i][j]])

for i in range(n):
    for j in range(i):
        p0_net[i][j] = p0_net[j][i]

p1_net = np.zeros_like(p1, dtype=int)

# Generate binary data based on the probabilities in the input matrix
for i in range(n):
    for j in range(i+1,n):
        p1_net[i][j] = np.random.choice([0, 1], p=[1 - p1[i][j], p1[i][j]])

for i in range(n):
    for j in range(i):
        p1_net[i][j] = p1_net[j][i]

p2_net = np.zeros_like(p2, dtype=int)

# Generate binary data based on the probabilities in the input matrix
for i in range(n):
    for j in range(i+1,n):
        p2_net[i][j] = np.random.choice([0, 1], p=[1 - p2[i][j], p2[i][j]])

for i in range(n):
    for j in range(i):
        p2_net[i][j] = p2_net[j][i]


# Number of different networks to generate
B =100

# Generate B different networks based on the probability matrix and store their embeddings
G = nx.Graph(p0_net)
embedding_results = []
for _ in range(B):
    # Generate a stochastic block model graph based on the probability matrix   
    # Perform node2vec embedding
    node2vec = Node2Vec(G, dimensions=3, walk_length=80, num_walks=10,p=1,q=0.5)
    model = node2vec.fit(window=10, min_count=1, batch_words=4)
    # Record the learned embeddings in the list
    embeddings = np.array([model.wv[str(node)] for node in G.nodes()])
    embedding_results.append(embeddings)

merged_matrix_2logn = np.concatenate(embedding_results, axis=1)
df = pd.DataFrame(merged_matrix_2logn)  
df.to_csv('80_2logn_embedding_results.csv', index=False)

  
G = nx.Graph(p1_net)
embedding_results = []
for _ in range(B):
    # Generate a stochastic block model graph based on the probability matrix   
    # Perform node2vec embedding
    node2vec = Node2Vec(G, dimensions=3, walk_length=80, num_walks=10,p=1,q=0.5)
    model = node2vec.fit(window=10, min_count=1, batch_words=4)
    # Record the learned embeddings in the list
    embeddings = np.array([model.wv[str(node)] for node in G.nodes()])
    embedding_results.append(embeddings)

merged_matrix_n12 = np.concatenate(embedding_results, axis=1)
df = pd.DataFrame(merged_matrix_n12)  
df.to_csv('80_n12_embedding_results.csv', index=False)


G = nx.Graph(p2_net)
embedding_results = []
for _ in range(B):
    # Generate a stochastic block model graph based on the probability matrix   
    # Perform node2vec embedding
    node2vec = Node2Vec(G, dimensions=3, walk_length=80, num_walks=10,p=1,q=0.5)
    model = node2vec.fit(window=10, min_count=1, batch_words=4)
    # Record the learned embeddings in the list
    embeddings = np.array([model.wv[str(node)] for node in G.nodes()])
    embedding_results.append(embeddings)

merged_matrix_n23 = np.concatenate(embedding_results, axis=1)
df = pd.DataFrame(merged_matrix_n23)  
df.to_csv('80_n23_embedding_results.csv', index=False)

