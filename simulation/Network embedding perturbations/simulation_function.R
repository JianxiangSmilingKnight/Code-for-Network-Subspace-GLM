### It is suggested to have X standardized before model fitting
library(RSpectra)
library(Matrix)

net.gen.from.P <- function(P,mode="undirected"){
  n <- nrow(P)
  if(mode=="undirected"){
    upper.index <- which(upper.tri(P))
    upper.p <- P[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    
    
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    
  }else{
    A <- matrix(0,n,n)
    R <- matrix(runif(n^2),n,n)
    A[R<P] <- 1
    diag(A) <- 0
  }
  return(A)
}

node2vec.SVD<-function(A,tL,tU,L){
  n<-nrow(A)
  gamma<-(L-tL-tU)*(tU-tL+1)/2
  A_sum<-sum(A)
  sqrt_DA_inverse<-diag(1/sqrt(apply(A, 1, sum)))
  W<-sqrt_DA_inverse%*%A%*%sqrt_DA_inverse
  W[is.na(W)] <- 0
  svd_result <- svd(W)
  tmp<-rep(0,n)
  for(i in tL:tU)
  {
    tmp<-tmp+(L-i)*svd_result$d^i
  }
  tmp<-2*A_sum/gamma*sqrt_DA_inverse%*%svd_result$u %*%diag(tmp)%*% t(svd_result$v)%*%sqrt_DA_inverse
  M_0<-log(tmp)
  return(M_0)
}

SP.Inf <- function(X,Y,A,K,r=NULL,thr=NULL,boot.thr=TRUE,boot.n= 50){
  n <- nrow(X)
  p <- ncol(X)
  eigen.A <- eigs_sym(A=A,k=K,which = "LA")
  X.svd <- svd(X)
  x.proj <- X.svd$v%*%(t(X.svd$u)/X.svd$d)
  R <- X.svd$u
  Uhat <- matrix(eigen.A$vectors[,1:K],ncol=K)
  hat.SVD <- svd(t(R)%*%Uhat,nv=K,nu=p)
  Vhat <- hat.SVD$v
  U.curl <- sqrt(n)*Uhat%*%Vhat
  R.curl <- sqrt(n)*R%*%hat.SVD$u
  if(is.null(r)){
  if(is.null(thr)){
    if(boot.thr){
      P0 <- eigen.A$vectors %*%t(eigen.A$vectors*eigen.A$values)
      P0[P0>1] <- 1
      P0[P0<0] <- 0
      P0 <- P0/mean(rowSums(P0))*mean(rowSums(A))
      eig.P0 <- eigs_sym(P0,k=K)
      fake.svd <- svd(t(R)%*%matrix(eig.P0$vectors[,1:K],ncol=K))
      print("Use bootstrapping to find r-threshold....")
      sim.s <- matrix(0,boot.n,min(p,K))
      for(II in 1:boot.n){
        A.sim <- net.gen.from.P(P0)
        eig.A.sim <- eigs_sym(A.sim,k=K)
        sim.svd <- svd(t(R)%*%matrix(eig.A.sim$vectors[,1:K],ncol=K))
        sim.s[II,] <- sim.svd$d
      }
      
      thr <- 1-max(abs(t(sim.s)-fake.svd$d))
      print(paste("Select r by threshold",thr))
      
    }else{
      dhat <- max(rowSums(A))
      thr <- 1- 4*sqrt(p*K*log(n))/dhat
      if(thr < 0.05){
        thr <- 0.05
      }
      print(paste("Select r by asymptotic threshold",thr))
    }
  }
    print("Observed singular values are ")
    print(hat.SVD$d)
    r <- sum(hat.SVD$d>=thr)
    print(paste("Select r =",r))
  }
  if(r >0){
  Z <- cbind(R.curl,U.curl[,-(1:r)])
  fit<-glm(Y~0+Z,family = "binomial")
  theta<-as.matrix(fit$coefficients[1:r])
  beta<-as.matrix(fit$coefficients[(r+1):p])
  VCOV<-vcov(fit)
  cov_theta<-VCOV[1:r,1:r]
  cov_beta<-VCOV[(r+1):p,(r+1):p]
  Xtheta<-R.curl[,(1:r)]%*%theta
  Xbeta<-R.curl[,(r+1):p]%*%beta
  theta.hat<-x.proj%*%Xtheta
  cov_theta<-x.proj%*%R.curl[,(1:r)]%*%cov_theta%*%t(x.proj%*%R.curl[,(1:r)])
  }else{
    Z <- cbind(R.curl,U.curl)
    fit<-glm(Y~0+Z,family = "binomial")
    beta<-as.matrix(fit$coefficients[1:p])
    VCOV<-vcov(fit)
    cov_beta<-VCOV[1:p,1:p]
    Xbeta<-R.curl[,1:p]%*%beta
  }
  beta.hat <- x.proj%*%Xbeta
  cov_beta<-x.proj%*%R.curl[,(r+1):p]%*%cov_beta%*%t(x.proj%*%R.curl[,(r+1):p])
  diag.sd<-sqrt(diag(cov_beta))
  CI.lower <- beta.hat - qnorm(1-0.05/2)*diag.sd
  CI.upper <- beta.hat + qnorm(1-0.05/2)*diag.sd
  SEE<-diag(cov_beta)
  cov_network<-VCOV[(p+1):(p+K-r),(p+1):(p+K-r)]
  network_effect_test_statistics<-t(as.matrix(fit$coefficients[(p+1):(p+K-r)]))%*%solve(cov_network)%*%as.matrix(fit$coefficients[(p+1):(p+K-r)])
  network<-(network_effect_test_statistics>qchisq(1-0.05,K-r))
  chisq.pval<-pchisq(network_effect_test_statistics,df=K-r,lower.tail=FALSE)
  return(list(beta=beta.hat,theta=theta.hat,
              fitted=fit$fitted.values,CI.lower=CI.lower,CI.upper=CI.upper,exist_network=network,SEE=SEE,model=fit,chisq.value= chisq.pval))
}

SP.Inf.Poisson <- function(X,Y,A,K,r=NULL,thr=NULL,boot.thr=TRUE,boot.n= 50){
  n <- nrow(X)
  p <- ncol(X)
  eigen.A <- eigs_sym(A=A,k=K,which = "LA")
  X.svd <- svd(X)
  x.proj <- X.svd$v%*%(t(X.svd$u)/X.svd$d)
  R <- X.svd$u
  Uhat <- matrix(eigen.A$vectors[,1:K],ncol=K)
  hat.SVD <- svd(t(R)%*%Uhat,nv=K,nu=p)
  Vhat <- hat.SVD$v
  U.curl <- sqrt(n)*Uhat%*%Vhat
  R.curl <- sqrt(n)*R%*%hat.SVD$u
  if(is.null(r)){
    if(is.null(thr)){
      if(boot.thr){
        P0 <- eigen.A$vectors %*%t(eigen.A$vectors*eigen.A$values)
        P0[P0>1] <- 1
        P0[P0<0] <- 0
        P0 <- P0/mean(rowSums(P0))*mean(rowSums(A))
        eig.P0 <- eigs_sym(P0,k=K)
        fake.svd <- svd(t(R)%*%matrix(eig.P0$vectors[,1:K],ncol=K))
        print("Use bootstrapping to find r-threshold....")
        sim.s <- matrix(0,boot.n,min(p,K))
        for(II in 1:boot.n){
          A.sim <- net.gen.from.P(P0)
          eig.A.sim <- eigs_sym(A.sim,k=K)
          sim.svd <- svd(t(R)%*%matrix(eig.A.sim$vectors[,1:K],ncol=K))
          sim.s[II,] <- sim.svd$d
        }
        
        thr <- 1-max(abs(t(sim.s)-fake.svd$d))
        print(paste("Select r by threshold",thr))
        
      }else{
        dhat <- max(rowSums(A))
        thr <- 1- 4*sqrt(p*K*log(n))/dhat
        if(thr < 0.05){
          thr <- 0.05
        }
        print(paste("Select r by asymptotic threshold",thr))
      }
    }
    print("Observed singular values are ")
    print(hat.SVD$d)
    r <- sum(hat.SVD$d>=thr)
    print(paste("Select r =",r))
  }
  if(r >0){
    Z <- cbind(R.curl,U.curl[,-(1:r)])
    fit<-glm(Y~0+Z,family = poisson(link = "log"))
    theta<-as.matrix(fit$coefficients[1:r])
    beta<-as.matrix(fit$coefficients[(r+1):p])
    VCOV<-vcov(fit)
    cov_theta<-VCOV[1:r,1:r]
    cov_beta<-VCOV[(r+1):p,(r+1):p]
    Xtheta<-R.curl[,(1:r)]%*%theta
    Xbeta<-R.curl[,(r+1):p]%*%beta
    theta.hat<-x.proj%*%Xtheta
    cov_theta<-x.proj%*%R.curl[,(1:r)]%*%cov_theta%*%t(x.proj%*%R.curl[,(1:r)])
  }else{
    Z <- cbind(R.curl,U.curl)
    fit<-glm(Y~0+Z,family = poisson(link = "log"))
    beta<-as.matrix(fit$coefficients[1:p])
    VCOV<-vcov(fit)
    cov_beta<-VCOV[1:p,1:p]
    Xbeta<-R.curl[,1:p]%*%beta
  }
  beta.hat <- x.proj%*%Xbeta
  cov_beta<-x.proj%*%R.curl[,(r+1):p]%*%cov_beta%*%t(x.proj%*%R.curl[,(r+1):p])
  diag.sd<-sqrt(diag(cov_beta))
  CI.lower <- beta.hat - qnorm(1-0.05/2)*diag.sd
  CI.upper <- beta.hat + qnorm(1-0.05/2)*diag.sd
  SEE<-diag(cov_beta)
  cov_network<-VCOV[(p+1):(p+K-r),(p+1):(p+K-r)]
  network_effect_test_statistics<-t(as.matrix(fit$coefficients[(p+1):(p+K-r)]))%*%solve(cov_network)%*%as.matrix(fit$coefficients[(p+1):(p+K-r)])
  network<-(network_effect_test_statistics>qchisq(1-0.05,K-r))
  chisq.pval<-pchisq(network_effect_test_statistics,df=K-r,lower.tail=FALSE)
  return(list(beta=beta.hat,theta=theta.hat,
              fitted=fit$fitted.values,CI.lower=CI.lower,CI.upper=CI.upper,exist_network=network,SEE=SEE,model=fit,chisq.value= chisq.pval))
}

SP.Inf.P <- function(X,Y,A,K,r=NULL,thr=NULL,boot.thr=TRUE,boot.n= 50){
  n <- nrow(X)
  p <- ncol(X)
  eigen.A <- eigs_sym(A=A,k=K,which = "LA")
  X.svd <- svd(X)
  x.proj <- X.svd$v%*%(t(X.svd$u)/X.svd$d)
  R <- X.svd$u
  Uhat <- matrix(eigen.A$vectors[,1:K],ncol=K)
  hat.SVD <- svd(t(R)%*%Uhat,nv=K,nu=p)
  Vhat <- hat.SVD$v
  U.curl <- sqrt(n)*Uhat%*%Vhat
  R.curl <- sqrt(n)*R%*%hat.SVD$u
  if(is.null(r)){
    if(is.null(thr)){
      if(boot.thr){
        P0 <- eigen.A$vectors %*%t(eigen.A$vectors*eigen.A$values)
        P0[P0>1] <- 1
        P0[P0<0] <- 0
        P0 <- P0/mean(rowSums(P0))*mean(rowSums(A))
        eig.P0 <- eigs_sym(P0,k=K)
        fake.svd <- svd(t(R)%*%matrix(eig.P0$vectors[,1:K],ncol=K))
        print("Use bootstrapping to find r-threshold....")
        sim.s <- matrix(0,boot.n,min(p,K))
        for(II in 1:boot.n){
          A.sim <- net.gen.from.P(P0)
          eig.A.sim <- eigs_sym(A.sim,k=K)
          sim.svd <- svd(t(R)%*%matrix(eig.A.sim$vectors[,1:K],ncol=K))
          sim.s[II,] <- sim.svd$d
        }
        
        thr <- 1-max(abs(t(sim.s)-fake.svd$d))
        print(paste("Select r by threshold",thr))
        
      }else{
        dhat <- max(rowSums(A))
        thr <- 1- 4*sqrt(p*K*log(n))/dhat
        if(thr < 0.05){
          thr <- 0.05
        }
        print(paste("Select r by asymptotic threshold",thr))
      }
    }
    print("Observed singular values are ")
    print(hat.SVD$d)
    r <- sum(hat.SVD$d>=thr)
    print(paste("Select r =",r))
  }
  if(r >0){
    Z <- cbind(R.curl,U.curl[,-(1:r)])
    fit<-glm(Y~0+Z,family = poisson(link = "log"))
    theta<-as.matrix(fit$coefficients[1:r])
    beta<-as.matrix(fit$coefficients[(r+1):p])
    VCOV<-vcov(fit)
    cov_theta<-VCOV[1:r,1:r]
    cov_beta<-VCOV[(r+1):p,(r+1):p]
    Xtheta<-R.curl[,(1:r)]%*%theta
    Xbeta<-R.curl[,(r+1):p]%*%beta
    theta.hat<-x.proj%*%Xtheta
    cov_theta<-x.proj%*%R.curl[,(1:r)]%*%cov_theta%*%t(x.proj%*%R.curl[,(1:r)])
  }else{
    Z <- cbind(R.curl,U.curl)
    fit<-glm(Y~0+Z,family = poisson(link = "log"))
    beta<-as.matrix(fit$coefficients[1:p])
    VCOV<-vcov(fit)
    cov_beta<-VCOV[1:p,1:p]
    Xbeta<-R.curl[,1:p]%*%beta
    theta.hat<-NULL
  }
  beta.hat <- x.proj%*%Xbeta
  cov_beta<-x.proj%*%R.curl[,(r+1):p]%*%cov_beta%*%t(x.proj%*%R.curl[,(r+1):p])
  diag.sd<-sqrt(diag(cov_beta))
  CI.lower <- beta.hat - qnorm(1-0.05/2)*diag.sd
  CI.upper <- beta.hat + qnorm(1-0.05/2)*diag.sd
  SEE<-diag(cov_beta)
  cov_network<-VCOV[(p+1):(p+K-r),(p+1):(p+K-r)]
  network_effect_test_statistics<-t(as.matrix(fit$coefficients[(p+1):(p+K-r)]))%*%solve(cov_network)%*%as.matrix(fit$coefficients[(p+1):(p+K-r)])
  network<-(network_effect_test_statistics>qchisq(1-0.05,K-r))
  chisq.pval<-pchisq(network_effect_test_statistics,df=K-r,lower.tail=FALSE)
  alpha<-U.curl[,(r+1):K]%*%as.matrix(fit$coefficients[(p+1):(p+K-r)])
  abs.t.val <- abs(beta.hat/diag.sd)
  t.pval <- 1-pnorm(abs.t.val)
  coef.mat <- cbind(beta.hat,CI.lower,CI.upper,t.pval)
  colnames(coef.mat) <- c("coefficient","CI-lower","CI-upper","p-val")
  return(list(beta=beta.hat,theta=theta.hat,alpha=alpha,
              fitted=fit$fitted.values,coef.mat=coef.mat,diag.sd=diag.sd,CI.lower=CI.lower,CI.upper=CI.upper,exist_network=network,chisq.value= chisq.pval,SEE=SEE,model=fit))
}


