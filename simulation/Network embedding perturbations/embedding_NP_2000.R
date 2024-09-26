# This is network embedding results with naive poisson regression

library(randnet)
library(RSpectra)
library(foreach)
library(doParallel)
registerDoParallel(cores=50)
source("simulation_function.R")
set.seed(1)


csv_record_2000_2logn<-read.csv("node2vec_A/SBM_2000_3_0.5/80_2logn_embedding_results.csv")
record_2000_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_2logn[[i]]<-csv_record_2000_2logn[,(3*i-2):(3*i)]
}
projection_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_2logn<-projection_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}

eigen.P <- eigs_sym(A=projection_2logn,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_2logn[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_2logn,file="Node2Vec_80_estimate_2000_2logn.Rda")


csv_record_2000_n12<-read.csv("node2vec_A/SBM_2000_3_0.5/80_n12_embedding_results.csv")
record_2000_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n12[[i]]<-csv_record_2000_n12[,(3*i-2):(3*i)]
}
projection_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n12<-projection_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_n12,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n12<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n12[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n12,file="Node2Vec_80_estimate_2000_n12.Rda")

csv_record_2000_n23<-read.csv("node2vec_A/SBM_2000_3_0.5/80_n23_embedding_results.csv")
record_2000_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n23[[i]]<-csv_record_2000_n23[,(3*i-2):(3*i)]
}
projection_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n23<-projection_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_n23,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n23<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n23[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n23,file="Node2Vec_80_estimate_2000_n23.Rda")


csv_record_2000_2logn<-read.csv("node2vec_A/SBM_2000_3_1/80_2logn_embedding_results.csv")
record_2000_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_2logn[[i]]<-csv_record_2000_2logn[,(3*i-2):(3*i)]
}
projection_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_2logn<-projection_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_2logn,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_2logn[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_2logn,file="DeepWalk_80_SBM_estimate_2000_2logn.Rda")


csv_record_2000_n12<-read.csv("node2vec_A/SBM_2000_3_1/80_n12_embedding_results.csv")
record_2000_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n12[[i]]<-csv_record_2000_n12[,(3*i-2):(3*i)]
}
projection_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n12<-projection_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_n12,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n12<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n12[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n12,file="DeepWalk_80_SBM_estimate_2000_n12.Rda")

csv_record_2000_n23<-read.csv("node2vec_A/SBM_2000_3_1/80_n23_embedding_results.csv")
record_2000_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n23[[i]]<-csv_record_2000_n23[,(3*i-2):(3*i)]
}
projection_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n23<-projection_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}



eigen.P <- eigs_sym(A=projection_n23,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n23<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n23[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n23,file="DeepWalk_80_SBM_estimate_2000_n23.Rda")


csv_record_2000_2logn<-read.csv("node2vec_A/DCBM_2000_3_0.5/80_2logn_embedding_results.csv")
record_2000_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_2logn[[i]]<-csv_record_2000_2logn[,(3*i-2):(3*i)]
}
projection_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_2logn<-projection_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_2logn,k=6,which = "LA")
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
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_2logn[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_2logn,file="Node2Vec_DCBM_80_estimate_2000_2logn.Rda")


csv_record_2000_n12<-read.csv("node2vec_A/DCBM_2000_3_0.5/80_n12_embedding_results.csv")
record_2000_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n12[[i]]<-csv_record_2000_n12[,(3*i-2):(3*i)]
}
projection_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n12<-projection_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}

eigen.P <- eigs_sym(A=projection_n12,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n12<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n12[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n12,file="Node2Vec_DCBM_80_estimate_2000_n12.Rda")

csv_record_2000_n23<-read.csv("node2vec_A/DCBM_2000_3_0.5/80_n23_embedding_results.csv")
record_2000_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n23[[i]]<-csv_record_2000_n23[,(3*i-2):(3*i)]
}
projection_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n23<-projection_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_n23,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n23<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n23[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n23,file="Node2Vec_DCBM_80_estimate_2000_n23.Rda")

csv_record_2000_2logn<-read.csv("node2vec_A/DCBM_2000_3_1/80_2logn_embedding_results.csv")
record_2000_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_2logn[[i]]<-csv_record_2000_2logn[,(3*i-2):(3*i)]
}
projection_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_2logn<-projection_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_2logn,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]

X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_2logn[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_2logn,file="DW_DCBM_80_estimate_2000_2logn.Rda")


csv_record_2000_n12<-read.csv("node2vec_A/DCBM_2000_3_1/80_n12_embedding_results.csv")
record_2000_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n12[[i]]<-csv_record_2000_n12[,(3*i-2):(3*i)]
}
projection_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n12<-projection_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_n12,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n12<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n12[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n12,file="DW_DCBM_80_estimate_2000_n12.Rda")

csv_record_2000_n23<-read.csv("node2vec_A/DCBM_2000_3_1/80_n23_embedding_results.csv")
record_2000_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n23[[i]]<-csv_record_2000_n23[,(3*i-2):(3*i)]
}
projection_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n23<-projection_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_n23,k=6,which = "LA")
eigen.P$vectors[,1]<--eigen.P$vectors[,1]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n23<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n23[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n23,file="DW_DCBM_80_estimate_2000_n23.Rda")

csv_record_2000_Diff2Vec_2logn<-read.csv("diff2vec_A/SBM_2000/20_Diff2Vec_2logn_embedding_results.csv")
record_2000_Diff2Vec_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_2logn[[i]]<-csv_record_2000_Diff2Vec_2logn[,(3*i-2):(3*i)]
}
projection_Diff2Vec_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_2logn<-projection_Diff2Vec_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_Diff2Vec_2logn,k=6,which = "LA")
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_2logn<-matrix(0,nrow = B,ncol = 1)
for (i in c(1:B)) {
  beta_hat_result<-matrix(0,nrow = M,ncol = 2)
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_2logn[i,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_2logn,file="diff2vec_SBM_20_estimate_2000_Diff2Vec_2logn.Rda")


csv_record_2000_Diff2Vec_n12<-read.csv("diff2vec_A/SBM_2000/20_Diff2Vec_n12_embedding_results.csv")
record_2000_Diff2Vec_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_n12[[i]]<-csv_record_2000_Diff2Vec_n12[,(3*i-2):(3*i)]
}
projection_Diff2Vec_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_n12<-projection_Diff2Vec_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_Diff2Vec_n12,k=6,which = "LA")
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_n12<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_n12[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_n12,file="diff2vec_SBM_20_estimate_2000_Diff2Vec_n12.Rda")

csv_record_2000_Diff2Vec_n23<-read.csv("diff2vec_A/SBM_2000/20_Diff2Vec_n23_embedding_results.csv")
record_2000_Diff2Vec_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_n23[[i]]<-csv_record_2000_Diff2Vec_n23[,(3*i-2):(3*i)]
}
projection_Diff2Vec_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_n23<-projection_Diff2Vec_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_Diff2Vec_n23,k=6,which = "LA")
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_n23<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_n23[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_n23,file="diff2vec_SBM_20_estimate_2000_Diff2Vec_n23.Rda")

csv_record_2000_Diff2Vec_2logn<-read.csv("diff2vec_A/DCBM_2000/20_Diff2Vec_2logn_embedding_results.csv")
record_2000_Diff2Vec_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_2logn[[i]]<-csv_record_2000_Diff2Vec_2logn[,(3*i-2):(3*i)]
}
projection_Diff2Vec_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_2logn<-projection_Diff2Vec_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}
svd_result <- svd(projection_Diff2Vec_2logn)
reconstructed_matrix <- svd_result$u[, 1:ncol(col_space)] %*% t(svd_result$v[, 1:ncol(col_space)])
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]-reconstructed_matrix),"2")
  }
  
}
print(mean(perturbation))
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]%*%svd_result$u[, 4]-reconstructed_matrix%*%svd_result$u[, 4]),"2")
  }
}
print(mean(perturbation))

eigen.P <- eigs_sym(A=projection_Diff2Vec_2logn,k=6,which = "LA")
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_2logn<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_2logn[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_2logn,file="Diff2Vec_DCBM_20_estimate_2000_Diff2Vec_2logn.Rda")

csv_record_2000_Diff2Vec_n12<-read.csv("diff2vec_A/DCBM_2000/20_Diff2Vec_n12_embedding_results.csv")
record_2000_Diff2Vec_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_n12[[i]]<-csv_record_2000_Diff2Vec_n12[,(3*i-2):(3*i)]
}
projection_Diff2Vec_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_n12<-projection_Diff2Vec_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}
svd_result <- svd(projection_Diff2Vec_n12)
reconstructed_matrix <- svd_result$u[, 1:ncol(col_space)] %*% t(svd_result$v[, 1:ncol(col_space)])
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]-reconstructed_matrix),"2")
  }
  
}
print(mean(perturbation))
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]%*%svd_result$u[, 4]-reconstructed_matrix%*%svd_result$u[, 4]),"2")
  }
}
print(mean(perturbation))

eigen.P <- eigs_sym(A=projection_Diff2Vec_n12,k=6,which = "LA")
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_n12<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_n12[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_n12,file="Diff2Vec_DCBM_20_estimate_2000_Diff2Vec_n12.Rda")

csv_record_2000_Diff2Vec_n23<-read.csv("diff2vec_A/DCBM_2000/20_Diff2Vec_n23_embedding_results.csv")
record_2000_Diff2Vec_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_n23[[i]]<-csv_record_2000_Diff2Vec_n23[,(3*i-2):(3*i)]
}
projection_Diff2Vec_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_n23<-projection_Diff2Vec_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}
svd_result <- svd(projection_Diff2Vec_n23)
reconstructed_matrix <- svd_result$u[, 1:ncol(col_space)] %*% t(svd_result$v[, 1:ncol(col_space)])
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]-reconstructed_matrix),"2")
  }
  
}
print(mean(perturbation))
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]%*%svd_result$u[, 4]-reconstructed_matrix%*%svd_result$u[, 4]),"2")
  }
}
print(mean(perturbation))


eigen.P <- eigs_sym(A=projection_Diff2Vec_n23,k=6,which = "LA")
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_n23<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_n23[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_n23,file="Diff2Vec_DCBM_20_estimate_2000_Diff2Vec_n23.Rda")

csv_record_2000_2logn<-read.csv("node2vec_A/Diag_2000_1/80_2logn_embedding_results.csv")
record_2000_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_2logn[[i]]<-csv_record_2000_2logn[,(3*i-2):(3*i)]
}
projection_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_2logn<-projection_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}
svd_result <- svd(projection_2logn)
reconstructed_matrix <- svd_result$u[, 1:ncol(col_space)] %*% t(svd_result$v[, 1:ncol(col_space)])
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]-reconstructed_matrix),"2")
  }
  
}
print(mean(perturbation))
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]%*%svd_result$u[, 4]-reconstructed_matrix%*%svd_result$u[, 4]),"2")
  }
  
}
print(mean(perturbation))

eigen.P <- eigs_sym(A=projection_2logn,k=6,which = "LA")
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_2logn<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_2logn[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_2logn,file="DeepWalk_Diag_80_estimate_2000_2logn.Rda")


csv_record_2000_n12<-read.csv("node2vec_A/Diag_2000_1/80_n12_embedding_results.csv")
record_2000_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n12[[i]]<-csv_record_2000_n12[,(3*i-2):(3*i)]
}
projection_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n12<-projection_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}
svd_result <- svd(projection_n12)
reconstructed_matrix <- svd_result$u[, 1:ncol(col_space)] %*% t(svd_result$v[, 1:ncol(col_space)])
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]-reconstructed_matrix),"2")
  }
  
}
print(mean(perturbation))
perturbation<-rep(0,100)
for (i in 1:100) {
  tmp<-matrix_list[[i]]
  if(is.null(tmp)==FALSE){
    perturbation[i]<-norm((matrix_list[[i]]%*%svd_result$u[, 4]-reconstructed_matrix%*%svd_result$u[, 4]),"2")
  }
}
print(mean(perturbation))

eigen.P <- eigs_sym(A=projection_n12,k=6,which = "LA")
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n12<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n12[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n12,file="DeepWalk_Diag_80_estimate_2000_n12.Rda")

csv_record_2000_n23<-read.csv("node2vec_A/Diag_2000_1/80_n23_embedding_results.csv")
record_2000_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n23[[i]]<-csv_record_2000_n23[,(3*i-2):(3*i)]
}
projection_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n23<-projection_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_n23,k=6,which = "LA")
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n23<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n23[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n23,file="DeepWalk_Diag_80_estimate_2000_n23.Rda")

csv_record_2000_2logn<-read.csv("node2vec_A/Diag_2000_0.5_tmp/80_2logn_embedding_results.csv")
record_2000_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_2logn[[i]]<-csv_record_2000_2logn[,(3*i-2):(3*i)]
}
projection_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_2logn<-projection_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}

eigen.P <- eigs_sym(A=projection_2logn,k=6,which = "LA")
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_2logn<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_2logn[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_2logn,file="Node2Vec_Diag_80_estimate_2000_2logn.Rda")

csv_record_2000_2logn<-read.csv("node2vec_A/Diag_2000_0.5/80_2logn_embedding_results.csv")
record_2000_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_2logn[[i]]<-csv_record_2000_2logn[,(3*i-2):(3*i)]
}
projection_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_2logn<-projection_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_2logn,k=6,which = "LA")
eigen.P$vectors[,2]<--eigen.P$vectors[,2]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_2logn<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_2logn[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_2logn,file="Node2Vec_Diag_80_estimate_2000_2logn.Rda")


csv_record_2000_n12<-read.csv("node2vec_A/Diag_2000_0.5/80_n12_embedding_results.csv")
record_2000_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n12[[i]]<-csv_record_2000_n12[,(3*i-2):(3*i)]
}
projection_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n12<-projection_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}

eigen.P <- eigs_sym(A=projection_n12,k=6,which = "LA")
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n12<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n12[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n12,file="Node2Vec_Diag_80_estimate_2000_n12.Rda")

csv_record_2000_n23<-read.csv("node2vec_A/Diag_2000_0.5/80_n23_embedding_results.csv")
record_2000_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_n23[[i]]<-csv_record_2000_n23[,(3*i-2):(3*i)]
}
projection_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_n23<-projection_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}


eigen.P <- eigs_sym(A=projection_n23,k=6,which = "LA")
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_n23<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_n23[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_n23,file="Node2Vec_Diag_80_estimate_2000_n23.Rda")

csv_record_2000_Diff2Vec_2logn<-read.csv("diff2vec_A/Diag2000/20_Diff2Vec_2logn_embedding_results.csv")
record_2000_Diff2Vec_2logn <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_2logn[[i]]<-csv_record_2000_Diff2Vec_2logn[,(3*i-2):(3*i)]
}
projection_Diff2Vec_2logn<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_2logn[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_2logn<-projection_Diff2Vec_2logn+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}

eigen.P <- eigs_sym(A=projection_Diff2Vec_2logn,k=6,which = "LA")
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_2logn<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_2logn[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_2logn,file="Diff2Vec_Diag_20_estimate_2000_Diff2Vec_2logn.Rda")


csv_record_2000_Diff2Vec_n12<-read.csv("diff2vec_A/Diag2000/20_Diff2Vec_n12_embedding_results.csv")
record_2000_Diff2Vec_n12 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_n12[[i]]<-csv_record_2000_Diff2Vec_n12[,(3*i-2):(3*i)]
}
projection_Diff2Vec_n12<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_n12[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_n12<-projection_Diff2Vec_n12+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}

eigen.P <- eigs_sym(A=projection_Diff2Vec_n12,k=6,which = "LA")
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_n12<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_n12[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_n12,file="Diff2Vec_Diag_20_estimate_2000_Diff2Vec_n12.Rda")

csv_record_2000_Diff2Vec_n23<-read.csv("diff2vec_A/Diag2000/20_Diff2Vec_n23_embedding_results.csv")
record_2000_Diff2Vec_n23 <- list()
B=100
n=2000
M=1000
for (i in 1:B) {
  record_2000_Diff2Vec_n23[[i]]<-csv_record_2000_Diff2Vec_n23[,(3*i-2):(3*i)]
}
projection_Diff2Vec_n23<-matrix(0,n,n)
matrix_list <- list()
for (i in 1:B) {
  tmp<-record_2000_Diff2Vec_n23[[i]]
  if(is.null(tmp)==FALSE){
    tmp <- as.matrix(tmp)
    projection_Diff2Vec_n23<-projection_Diff2Vec_n23+tmp%*%t(tmp)/B
    svd_result <- svd(tmp)
    col_space <- svd_result$u[, 1:length(svd_result$d)]
    matrix_list[[i]]<-col_space%*%t(col_space)
  }
}

eigen.P <- eigs_sym(A=projection_Diff2Vec_n23,k=6,which = "LA")
eigen.P$vectors[,3]<--eigen.P$vectors[,3]
eigen.P$vectors[,4]<--eigen.P$vectors[,4]
X.true1<-1/sqrt(25)*eigen.P$vectors[,2]+sqrt(24)/sqrt(25)*eigen.P$vectors[,4]
X.true <-sqrt(n)*cbind(eigen.P$vectors[,1],X.true1)
Xrho<-0.5
theta<- matrix(c(0.5,0),ncol=1)
beta<- matrix(c(0,Xrho),ncol=1)
Xtheta <- X.true%*%theta
Xbeta <- X.true%*%beta
W<-sqrt(n)*eigen.P$vectors[,1:3]
rho<-0.5
alpha.coef <-rho* matrix(c(0,1,1),ncol=1)
alpha <- W%*%alpha.coef
EY <-exp(Xtheta+Xbeta+alpha)
estimate_2000_Diff2Vec_n23<-matrix(0,nrow = B,ncol = 1)
for (j in c(1:B)) {
  
  tmp<-foreach(i=1:M,.packages="RSpectra",.combine = cbind) %dopar% {
    Y <- rpois(n,EY)
    fit <-glm(Y~X.true,family = poisson(link = "log"))
    SPE<-t(fit$fitted.values-EY)%*%(fit$fitted.values-EY)/n
    c(SPE)
  }
  tmp<-t(tmp)
  estimate_2000_Diff2Vec_n23[j,]<-c(mean(tmp[,1]))
}
save(estimate_2000_Diff2Vec_n23,file="Diff2Vec_Diag_20_estimate_2000_Diff2Vec_n23.Rda")