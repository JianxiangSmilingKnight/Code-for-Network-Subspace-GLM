.libPaths("/panfs/jay/groups/31/tianxili/wan01965/R/x86_64-pc-linux-gnu-library/4.1")
library(randnet)
library(netcoh)
library(RSpectra)
library(foreach)
library(doParallel)
registerDoParallel(cores=50)
source("/panfs/jay/groups/31/tianxili/wan01965/simulation_function.R")

set.seed(1)
N <-4000 ## network size
big.model <- BlockModel.Gen(lambda=2*log(N),n=N,beta=0.3,K=3)
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

M=500
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
EY <-(1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_500_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- as.matrix(rbinom(n,1,EY))
    lambda = seq(2,1,length.out=10)
    perm.rnc.cv <- rncregpath(A=as.matrix(A_0),X=X.true,Y=Y,lambdaseq=lambda,model="logistic",cv=5,cv.seed=999)
    perm.rnc.lambda <- lambda[perm.rnc.cv$cv.min.index]
    perm.rnc.fit <- rncreg(A=as.matrix(A_0),X=X.true,Y=Y,model="logistic",lambda=perm.rnc.lambda)
    fitted_value<-exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha)/(1+exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha))
    SPE<-as.matrix(t(fitted_value-EY)%*%(fitted_value-EY)/n)
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-mean(tmp)
}
save(record_500_2logn,file="SBM_L_record_500_2logn.Rda")

eigen.P <- eigs_sym(A=P_1,k=4,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
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
EY <-(1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_500_n1.2<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- as.matrix(rbinom(n,1,EY))
    lambda = seq(2,1,length.out=10)
    perm.rnc.cv <- rncregpath(A=as.matrix(A_1),X=X.true,Y=Y,lambdaseq=lambda,model="logistic",cv=5,cv.seed=999)
    perm.rnc.lambda <- lambda[perm.rnc.cv$cv.min.index]
    perm.rnc.fit <- rncreg(A=as.matrix(A_1),X=X.true,Y=Y,model="logistic",lambda=perm.rnc.lambda)
    fitted_value<-exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha)/(1+exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha))
    SPE<-as.matrix(t(fitted_value-EY)%*%(fitted_value-EY)/n)
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_n1.2[i,]<-mean(tmp)
}
save(record_500_n1.2,file="SBM_L_record_500_n1.2.Rda")


eigen.P <- eigs_sym(A=P_2,k=4,which = "LA")
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
EY <-(1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_500_n2.3<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- as.matrix(rbinom(n,1,EY))
    lambda = seq(2,1,length.out=10)
    perm.rnc.cv <- rncregpath(A=as.matrix(A_2),X=X.true,Y=Y,lambdaseq=lambda,model="logistic",cv=5,cv.seed=999)
    perm.rnc.lambda <- lambda[perm.rnc.cv$cv.min.index]
    perm.rnc.fit <- rncreg(A=as.matrix(A_2),X=X.true,Y=Y,model="logistic",lambda=perm.rnc.lambda)
    fitted_value<-exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha)/(1+exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha))
    SPE<-as.matrix(t(fitted_value-EY)%*%(fitted_value-EY)/n)
    c(SPE)
  }
  tmp<-t(tmp)
  record_500_n2.3[i,]<-mean(tmp)
}

save(record_500_n2.3,file="SBM_L_record_500_n2.3.Rda")

n <- 1000
P <- big.P[1:n,1:n]
P_0<-P/sum(P)*2*log(n)*n
P_1<-P_0*(n^{1/2}/2/log(n))
P_2<-P_0*(n^{2/3}/2/log(n))

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
EY <-(1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_1000_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- as.matrix(rbinom(n,1,EY))
    lambda = seq(2,1,length.out=10)
    perm.rnc.cv <- rncregpath(A=as.matrix(A_0),X=X.true,Y=Y,lambdaseq=lambda,model="logistic",cv=5,cv.seed=999)
    perm.rnc.lambda <- lambda[perm.rnc.cv$cv.min.index]
    perm.rnc.fit <- rncreg(A=as.matrix(A_0),X=X.true,Y=Y,model="logistic",lambda=perm.rnc.lambda)
    fitted_value<-exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha)/(1+exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha))
    SPE<-as.matrix(t(fitted_value-EY)%*%(fitted_value-EY)/n)
    c(SPE)
  }
  tmp<-t(tmp)
  record_1000_2logn[i,]<-mean(tmp)
}
save(record_1000_2logn,file="SBM_L_record_1000_2logn.Rda")

eigen.P <- eigs_sym(A=P_1,k=4,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
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
EY <-(1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_1000_n1.2<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- as.matrix(rbinom(n,1,EY))
    lambda = exp(seq(log(1000),log(0.001),length.out=10))
    perm.rnc.cv <- rncregpath(A=as.matrix(A_1),X=X.true,Y=Y,lambdaseq=lambda,model="logistic",cv=5,cv.seed=999)
    perm.rnc.lambda <- lambda[perm.rnc.cv$cv.min.index]
    perm.rnc.fit <- rncreg(A=as.matrix(A_1),X=X.true,Y=Y,model="logistic",lambda=perm.rnc.lambda)
    fitted_value<-exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha)/(1+exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha))
    SPE<-as.matrix(t(fitted_value-EY)%*%(fitted_value-EY)/n)
    c(SPE)
  }
  tmp<-t(tmp)
  record_1000_n1.2[i,]<-mean(tmp)
}
save(record_1000_n1.2,file="SBM_L_record_1000_n1.2.Rda")


eigen.P <- eigs_sym(A=P_2,k=4,which = "LA")
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
EY <-(1+exp(-Xtheta-Xbeta-alpha))^(-1)
record_1000_n2.3<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- as.matrix(rbinom(n,1,EY))
    lambda = exp(seq(log(1000),log(0.001),length.out=10))
    perm.rnc.cv <- rncregpath(A=as.matrix(A_2),X=X.true,Y=Y,lambdaseq=lambda,model="logistic",cv=5,cv.seed=999)
    perm.rnc.lambda <- lambda[perm.rnc.cv$cv.min.index]
    perm.rnc.fit <- rncreg(A=as.matrix(A_2),X=X.true,Y=Y,model="logistic",lambda=perm.rnc.lambda)
    fitted_value<-exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha)/(1+exp(X.true%*%perm.rnc.fit$beta+perm.rnc.fit$alpha))
    SPE<-as.matrix(t(fitted_value-EY)%*%(fitted_value-EY)/n)
    c(SPE)
  }
  tmp<-t(tmp)
  record_1000_n2.3[i,]<-mean(tmp)
}

save(record_1000_n2.3,file="SBM_L_record_1000_n2.3.Rda")