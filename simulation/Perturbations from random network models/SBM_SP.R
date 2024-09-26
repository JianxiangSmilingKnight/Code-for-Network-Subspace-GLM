# Subspace Regression for SBM

library(randnet)
library(RSpectra)
library(foreach)
library(doParallel)
registerDoParallel(cores=10)
source("simulation_function.R")

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
P_3<-P_0*(n/6/2/log(n))

M=500
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
record_500_2logn<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_0,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_500_2logn,file="Lrecord_500_2logn.Rda")
record_500_n1.2<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_1,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_500_n1.2[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_500_n1.2,file="Lrecord_500_n1.2.Rda")


record_500_n2.3<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_2,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_500_n2.3[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}

save(record_500_n2.3,file="Lrecord_500_n2.3.Rda")


record_500_n.6<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_3 <- net.gen.from.P(P_3)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_3,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_500_n.6[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_500_n.6,file="Lrecord_500_n.6.Rda")


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
record_1000_2logn<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_0,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_1000_2logn[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_1000_2logn,file="Lrecord_1000_2logn.Rda")
record_1000_n1.2<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_1,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_1000_n1.2[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_1000_n1.2,file="Lrecord_1000_n1.2.Rda")


record_1000_n2.3<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_2,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_1000_n2.3[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}

save(record_1000_n2.3,file="Lrecord_1000_n2.3.Rda")


record_1000_n.6<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_3 <- net.gen.from.P(P_3)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_3,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_1000_n.6[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_1000_n.6,file="Lrecord_1000_n.6.Rda")


n <- 2000
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
record_2000_2logn<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_0,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_2000_2logn[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_2000_2logn,file="Lrecord_2000_2logn.Rda")
record_2000_n1.2<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_1,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_2000_n1.2[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_2000_n1.2,file="Lrecord_2000_n1.2.Rda")


record_2000_n2.3<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_2,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_2000_n2.3[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}

save(record_2000_n2.3,file="Lrecord_2000_n2.3.Rda")


record_2000_n.6<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_3 <- net.gen.from.P(P_3)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_3,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_2000_n.6[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_2000_n.6,file="Lrecord_2000_n.6.Rda")

n <- 4000
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
record_4000_2logn<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_0,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_4000_2logn[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_4000_2logn,file="Lrecord_4000_2logn.Rda")
record_4000_n1.2<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_1,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_4000_n1.2[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_4000_n1.2,file="Lrecord_4000_n1.2.Rda")


record_4000_n2.3<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_2,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_4000_n2.3[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}

save(record_4000_n2.3,file="Lrecord_4000_n2.3.Rda")


record_4000_n.6<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_3 <- net.gen.from.P(P_3)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    fit <- SP.Inf(X.true,Y,A_3,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_4000_n.6[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_4000_n.6,file="Lrecord_4000_n.6.Rda")

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
P_3<-P_0*(n/6/2/log(n))

M=500
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
EY <- exp(Xtheta+Xbeta+alpha)
record_500_2logn<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_0,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_500_2logn[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_500_2logn,file="Precord_500_2logn.Rda")
record_500_n1.2<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_1,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_500_n1.2[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_500_n1.2,file="Precord_500_n1.2.Rda")


record_500_n2.3<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_2,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_500_n2.3[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}

save(record_500_n2.3,file="Precord_500_n2.3.Rda")


record_500_n.6<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_3 <- net.gen.from.P(P_3)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_3,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_500_n.6[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_500_n.6,file="Precord_500_n.6.Rda")


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
EY <- exp(Xtheta+Xbeta+alpha)
record_1000_2logn<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_0,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_1000_2logn[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_1000_2logn,file="Precord_1000_2logn.Rda")
record_1000_n1.2<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_1,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_1000_n1.2[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_1000_n1.2,file="Precord_1000_n1.2.Rda")


record_1000_n2.3<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_2,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_1000_n2.3[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}

save(record_1000_n2.3,file="Precord_1000_n2.3.Rda")


record_1000_n.6<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_3 <- net.gen.from.P(P_3)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_3,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_1000_n.6[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_1000_n.6,file="Precord_1000_n.6.Rda")


n <- 2000
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
EY <- exp(Xtheta+Xbeta+alpha)
record_2000_2logn<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_0,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_2000_2logn[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_2000_2logn,file="Precord_2000_2logn.Rda")
record_2000_n1.2<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_1,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_2000_n1.2[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_2000_n1.2,file="Precord_2000_n1.2.Rda")


record_2000_n2.3<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_2,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_2000_n2.3[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}

save(record_2000_n2.3,file="Precord_2000_n2.3.Rda")


record_2000_n.6<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_3 <- net.gen.from.P(P_3)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_3,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_2000_n.6[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_2000_n.6,file="Precord_2000_n.6.Rda")

n <- 4000
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
EY <- exp(Xtheta+Xbeta+alpha)
record_4000_2logn<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_0 <- net.gen.from.P(P_0)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_0,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_4000_2logn[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_4000_2logn,file="Precord_4000_2logn.Rda")
record_4000_n1.2<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_1 <- net.gen.from.P(P_1)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_1,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_4000_n1.2[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_4000_n1.2,file="Precord_4000_n1.2.Rda")


record_4000_n2.3<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_2 <- net.gen.from.P(P_2)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_2,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_4000_n2.3[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}

save(record_4000_n2.3,file="Precord_4000_n2.3.Rda")


record_4000_n.6<-matrix(0,nrow = B,ncol = 6)
for (i in c(1:B)) {
  A_3 <- net.gen.from.P(P_3)
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <- SP.Inf.Poisson(X.true,Y,A_3,K=3,r=1,boot.thr=FALSE)
    bias=fit$beta[2,]-Xrho
    count.CI<-(Xrho>fit$CI.lower[2])*(Xrho<fit$CI.upper[2])
    SPE<-t(fit$fitted-EY)%*%(fit$fitted-EY)/n
    SEE<-fit$SEE[2]
    c(bias,SEE,SPE,count.CI)
  }
  tmp<-t(tmp)
  record_4000_n.6[i,]<-c(mean(tmp[,1]),sd(tmp[,2]),mean(tmp[,1]^2),mean(tmp[,2]),mean(tmp[,3]),mean(tmp[,4]))
}
save(record_4000_n.6,file="Precord_4000_n.6.Rda")