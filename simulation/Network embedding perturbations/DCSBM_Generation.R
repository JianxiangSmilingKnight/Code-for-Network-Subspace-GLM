library(randnet)
library(RSpectra)
library(foreach)
library(doParallel)
library(glmnet)
source("simulation_function.R")

# This is to generate the DCBM network for network embedding
set.seed(1)
N <-2000 ## network size
big.model <- BlockModel.Gen(lambda=2*log(N),n=N,beta=0.3,K=3,rho = 0.8,power = FALSE,simple = FALSE)
W <- big.model$P
neighborhood_matrix <- W

# Generate a permutation vector
permutation_vector <- sample(nrow(neighborhood_matrix))

# Permute the rows and columns of the neighborhood matrix
big.P <- neighborhood_matrix[permutation_vector, permutation_vector]

n <- 2000
P <- big.P[1:n,1:n]
P_0<-P/sum(P)*2*log(n)*n
P_1<-P_0*(n^{1/2}/2/log(n))
P_2<-P_0*(n^{2/3}/2/log(n))
P_3<-P_0*(n/6/2/log(n))
write.csv(P_0, "P_0_2000.csv", row.names=FALSE)
write.csv(P_1, "P_1_2000.csv", row.names=FALSE)
write.csv(P_2, "P_2_2000.csv", row.names=FALSE)
