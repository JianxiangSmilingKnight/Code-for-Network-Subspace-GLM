library(randnet)
library(RSpectra)
library(foreach)
library(doParallel)
registerDoParallel(cores=90)
source("simulation_function.R")

set.seed(1)
N <-4000 ## network size
big.model <- BlockModel.Gen(lambda=2*log(N),n=N,beta=0.3,K=3,rho = 0.8,power = FALSE,simple = FALSE)
W <- big.model$P
neighborhood_matrix <- W

# Generate a permutation vector
permutation_vector <- sample(nrow(neighborhood_matrix))

# Permute the rows and columns of the neighborhood matrix
big.P <- neighborhood_matrix[permutation_vector, permutation_vector]

n <- 500
P <- big.P[1:n,1:n]
P_0<-P/sum(P)*2*log(n)*n
P_1<-P_0*(n^{1/2}/2/log(n))
P_2<-P_0*(n^{2/3}/2/log(n))
P_3<-P_0*(n/6/2/log(n))

M=500
B=100

eigen.P <- eigs_sym(A=P_0,k=4,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <- (1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_500_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-c(mean(tmp[,1]))
}
save(record_500_2logn,file="record_500_2logn.Rda")


record_500_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-c(mean(tmp[,1]))
}

n <- 1000
P <- big.P[1:n,1:n]
P_0<-P/sum(P)*2*log(n)*n
P_1<-P_0*(n^{1/2}/2/log(n))
P_2<-P_0*(n^{2/3}/2/log(n))
P_3<-P_0*(n/6/2/log(n))

M=1000
B=100

eigen.P <- eigs_sym(A=P_0,k=4,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <- (1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_1000_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_1000_2logn[i,]<-c(mean(tmp[,1]))
}
save(record_1000_2logn,file="record_1000_2logn.Rda")

set.seed(1)
N <-4000 ## network size
big.model <- BlockModel.Gen(lambda=2*log(N),n=N,beta=0.3,K=3,rho = 0.2)
W <- big.model$P
neighborhood_matrix <- W

# Generate a permutation vector
permutation_vector <- sample(nrow(neighborhood_matrix))

# Permute the rows and columns of the neighborhood matrix
big.P <- neighborhood_matrix[permutation_vector, permutation_vector]

n <- 500
P <- big.P[1:n,1:n]
P_0<-P/sum(P)*2*log(n)*n
P_1<-P_0*(n^{1/2}/2/log(n))
P_2<-P_0*(n^{2/3}/2/log(n))
P_3<-P_0*(n/6/2/log(n))

M=1000
B=100

eigen.P <- eigs_sym(A=P_0,k=4,which = "LA")

eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <- (1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_500_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-c(mean(tmp[,1]))
}
save(record_500_2logn,file="record_500_2logn.Rda")


record_500_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-c(mean(tmp[,1]))
}

n <- 1000
P <- big.P[1:n,1:n]
P_0<-P/sum(P)*2*log(n)*n
P_1<-P_0*(n^{1/2}/2/log(n))
P_2<-P_0*(n^{2/3}/2/log(n))
P_3<-P_0*(n/6/2/log(n))

M=1000
B=100

eigen.P <- eigs_sym(A=P_0,k=4,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <- (1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_1000_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_1000_2logn[i,]<-c(mean(tmp[,1]))
}
save(record_1000_2logn,file="record_1000_2logn.Rda")



set.seed(1)
N <-4000 ## network size
U = 0.8*matrix(1:N,nrow=1) / (N+1)
V = 0.8*matrix(1:N,nrow=1) / (N+1)
UU = t(U)%*%matrix(1,nrow=1,ncol=N)
VV = matrix(1,N,1)%*%V
W = 15*(abs(UU-VV))^(4/5)-0.1
W = 1-1/(1+0.75*exp(-W))
neighborhood_matrix <- W

# Generate a permutation vector
permutation_vector <- sample(nrow(neighborhood_matrix))

# Permute the rows and columns of the neighborhood matrix
big.P <- neighborhood_matrix[permutation_vector, permutation_vector]

n <- 500
P <- big.P[1:n,1:n]
P_0<-P/sum(P)*2*log(n)*n
P_1<-P_0*(n^{1/2}/2/log(n))
P_2<-P_0*(n^{2/3}/2/log(n))
P_3<-P_0*(n/6/2/log(n))

M=1000
B=100

eigen.P <- eigs_sym(A=P_0,k=4,which = "LA")
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <- (1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_500_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-c(mean(tmp[,1]))
}
save(record_500_2logn,file="record_500_2logn.Rda")


record_500_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-c(mean(tmp[,1]))
}

n <- 1000
P <- big.P[1:n,1:n]
P_0<-P/sum(P)*2*log(n)*n
P_1<-P_0*(n^{1/2}/2/log(n))
P_2<-P_0*(n^{2/3}/2/log(n))
P_3<-P_0*(n/6/2/log(n))

M=1000
B=100

eigen.P <- eigs_sym(A=P_0,k=4,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <- (1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_1000_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <-glm(Y~X.true,family = "binomial")
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  record_1000_2logn[i,]<-c(mean(tmp[,1]))
}
save(record_1000_2logn,file="record_1000_2logn.Rda")