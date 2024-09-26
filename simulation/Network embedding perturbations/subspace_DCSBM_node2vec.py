# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:55:18 2024

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
np.random.seed(2)
n=2000
p0 = pd.read_csv("C:/Users/User/.spyder-py3/P_0_2000.csv")
p0=p0.values
p1 = pd.read_csv("C:/Users/User/.spyder-py3/P_1_2000.csv")
p1=p1.values
p2 = pd.read_csv("C:/Users/User/.spyder-py3/P_2_2000.csv")
p2=p2.values

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
