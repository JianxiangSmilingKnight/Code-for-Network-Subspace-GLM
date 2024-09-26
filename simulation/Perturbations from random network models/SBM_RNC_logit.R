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
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  D <- diag(rowSums(A_0))
  gamma<-0.05
  L_gamma <- D - A_0 + diag(rep(gamma,n))
  eigen_decomp <- eigen(L_gamma, symmetric = TRUE)
  tmp<-foreach(i=1:M,.packages="glmnet",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    Xtmp<-cbind(X.true,(eigen_decomp$vectors))
    penalty_factors <- c(rep(0, 2), eigen_decomp$values)
    cv_fit_perm <- cv.glmnet(Xtmp, Y, family = "binomial",nfolds = 10, penalty.factor = penalty_factors,lambda = seq(2,1,length.out=10), alpha = 0)
    optimal_lambda_perm <- cv_fit_perm$lambda.min
    fit_perm <- glmnet(Xtmp, Y, family = "binomial",thresh = 1e-03, alpha = 0,intercept = FALSE, penalty.factor =penalty_factors,lambda = optimal_lambda_perm)
    fitted_value<-exp(Xtmp%*%fit_perm$beta)
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
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  D <- diag(rowSums(A_1))
  gamma<-0.05
  L_gamma <- D - A_1 + diag(rep(gamma,n))
  eigen_decomp <- eigen(L_gamma, symmetric = TRUE)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    Xtmp<-cbind(X.true,eigen_decomp$vectors)
    penalty_factors <- c(rep(0, 2), eigen_decomp$values)
    cv_fit_perm <- cv.glmnet(Xtmp, Y, family = "binomial",nfolds = 10, penalty.factor = penalty_factors, alpha = 0,lambda = seq(2,1,length.out=10))
    optimal_lambda_perm <- cv_fit_perm$lambda.min
    fit_perm <- glmnet(Xtmp, Y,  family = "binomial",thresh = 1e-03, alpha = 0,intercept = FALSE, penalty.factor =penalty_factors,lambda = optimal_lambda_perm)
    fitted_value<-exp(Xtmp%*%fit_perm$beta)
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
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  D <- diag(rowSums(A_2))
  gamma<-0.05
  L_gamma <- D - A_2 + diag(rep(gamma,n))
  eigen_decomp <- eigen(L_gamma, symmetric = TRUE)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    Xtmp<-cbind(X.true,(eigen_decomp$vectors))
    penalty_factors <- c(rep(0, 2), eigen_decomp$values)
    cv_fit_perm <- cv.glmnet(Xtmp, Y, family = "binomial",nfolds = 10, penalty.factor = penalty_factors, alpha = 0,lambda = seq(2,1,length.out=10))
    optimal_lambda_perm <- cv_fit_perm$lambda.min
    fit_perm <- glmnet(Xtmp, Y, family = "binomial",thresh = 1e-03, alpha = 0,intercept = FALSE, penalty.factor =penalty_factors,lambda = optimal_lambda_perm)
    fitted_value<-exp(Xtmp%*%fit_perm$beta)
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
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  D <- diag(rowSums(A_0))
  gamma<-0.05
  L_gamma <- D - A_0 + diag(rep(gamma,n))
  eigen_decomp <- eigen(L_gamma, symmetric = TRUE)
  tmp<-foreach(i=1:M,.packages="glmnet",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    Xtmp<-cbind(X.true,(eigen_decomp$vectors))
    penalty_factors <- c(rep(0, 2), eigen_decomp$values)
    cv_fit_perm <- cv.glmnet(Xtmp, Y, family = "binomial",nfolds = 10, penalty.factor = penalty_factors,lambda = seq(2,1,length.out=10), alpha = 0)
    optimal_lambda_perm <- cv_fit_perm$lambda.min
    fit_perm <- glmnet(Xtmp, Y, family = "binomial",thresh = 1e-03, alpha = 0,intercept = FALSE, penalty.factor =penalty_factors,lambda = optimal_lambda_perm)
    fitted_value<-exp(Xtmp%*%fit_perm$beta)
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
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  D <- diag(rowSums(A_1))
  gamma<-0.05
  L_gamma <- D - A_1 + diag(rep(gamma,n))
  eigen_decomp <- eigen(L_gamma, symmetric = TRUE)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    Xtmp<-cbind(X.true,eigen_decomp$vectors)
    penalty_factors <- c(rep(0, 2), eigen_decomp$values)
    cv_fit_perm <- cv.glmnet(Xtmp, Y, family = "binomial",nfolds = 10, penalty.factor = penalty_factors, alpha = 0,lambda = seq(2,1,length.out=10))
    optimal_lambda_perm <- cv_fit_perm$lambda.min
    fit_perm <- glmnet(Xtmp, Y,  family = "binomial",thresh = 1e-03, alpha = 0,intercept = FALSE, penalty.factor =penalty_factors,lambda = optimal_lambda_perm)
    fitted_value<-exp(Xtmp%*%fit_perm$beta)
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
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  D <- diag(rowSums(A_2))
  gamma<-0.05
  L_gamma <- D - A_2 + diag(rep(gamma,n))
  eigen_decomp <- eigen(L_gamma, symmetric = TRUE)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rbinom(n,1,EY)
    Xtmp<-cbind(X.true,(eigen_decomp$vectors))
    penalty_factors <- c(rep(0, 2), eigen_decomp$values)
    cv_fit_perm <- cv.glmnet(Xtmp, Y, family = "binomial",nfolds = 10, penalty.factor = penalty_factors, alpha = 0,lambda = seq(2,1,length.out=10))
    optimal_lambda_perm <- cv_fit_perm$lambda.min
    fit_perm <- glmnet(Xtmp, Y, family = "binomial",thresh = 1e-03, alpha = 0,intercept = FALSE, penalty.factor =penalty_factors,lambda = optimal_lambda_perm)
    fitted_value<-exp(Xtmp%*%fit_perm$beta)
    SPE<-as.matrix(t(fitted_value-EY)%*%(fitted_value-EY)/n)
    c(SPE)
  }
  tmp<-t(tmp)
  record_1000_n2.3[i,]<-mean(tmp)
}

save(record_1000_n2.3,file="SBM_L_record_1000_n2.3.Rda")