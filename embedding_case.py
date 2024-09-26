# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 20:02:17 2024

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
np.random.seed(1)

p0 = pd.read_csv("G.csv")
p0=p0.values




# Generate B different networks based on the probability matrix and store their embeddings

G = nx.Graph(p0)
B=1
    # Generate a stochastic block model graph based on the probability matrix   
    # Perform node2vec embedding
node2vec = Node2Vec(G, dimensions=34, walk_length=80, num_walks=10,p=1,q=1)
model = node2vec.fit(window=10, min_count=1, batch_words=4)
# Record the learned embeddings in the list
embeddings = np.array([model.wv[str(node)] for node in G.nodes()])
df = pd.DataFrame(embeddings)  
df.to_csv('embedding_case_study_G_node2vec.csv', index=False)


