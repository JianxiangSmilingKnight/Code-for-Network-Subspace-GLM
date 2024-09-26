library(latex2exp)
library(irlba)
library(glmnet)
library(RSpectra)
library(doParallel)
library(igraph)
library(randnet)
library(testit)
library(foreach)

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
    fit<-glm(Y~0+Z,family = "binomial",control = list(maxit = 100))
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
    fit<-glm(Y~0+Z,family = "binomial",control = list(maxit = 100))
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

SP.semi.Inf <- function(X,Y,A,K,r=NULL,thr=NULL,boot.thr=TRUE,boot.n= 50){
  n <- nrow(X)
  p <- ncol(X)
  N <- length(Y)
  if(N>= n){
    print("Not correct dimension for semi-supervised version....")
    return(NA)
  }
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
    fit<-glm(Y~0+Z[1:N,],family = "binomial",control = list(maxit = 100))
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
    fit<-glm(Y~0+Z[1:N,],family = "binomial",control = list(maxit = 100))
    beta<-as.matrix(fit$coefficients[1:p])
    VCOV<-vcov(fit)
    cov_beta<-VCOV[1:p,1:p]
    Xbeta<-R.curl%*%beta
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
  alpha.hat <- alpha[1:N,]
  full.fitted <-X[(N+1):n,]%*%beta.hat + U.curl[(1+N):n,(r+1):K]%*%as.matrix(fit$coefficients[(p+1):(p+K-r)])
  predict_alpha= U.curl[(1+N):n,(r+1):K]%*%as.matrix(fit$coefficients[(p+1):(p+K-r)])
  full.fitted<-exp(full.fitted)/(1+exp(full.fitted))
  return(list(beta=beta.hat,theta=theta.hat,alpha=alpha.hat,
              fitted=fit$fitted.values,CI.lower=CI.lower,CI.upper=CI.upper,exist_network=network,chisq.value= chisq.pval,full.fitted=full.fitted,predict_alpha=predict_alpha,SEE=SEE,model=fit))
}

logistic.back <- function(X,Y,pval.thr=0.1,adj=TRUE,keep.index){
  n <- nrow(X)
  p <- ncol(X)
  X.keep <- X[,keep.index]
  X.select <- X[,-keep.index]
  p.keep <- length(keep.index)
  new.X <- cbind(X.keep,X.select)
  select.terminate <- FALSE
  while(!select.terminate){
    p.current <- ncol(new.X)
    select.index <- setdiff(1:p.current,1:p.keep)
    
    fit.old <- glm(Y~new.X-1,family = "binomial",control = list(maxit = 100))
    ss <- summary(fit.old)
    orig.p.val <- ss$coefficients[select.index,4]
    if(adj){
      select.p.val <- p.adjust(orig.p.val,method ="bonferroni")
    }else{
      select.p.val <- orig.p.val
    }
    elim.index <- select.index[which.max(select.p.val)]
    if(max(select.p.val)<=pval.thr){
      select.terminate <- TRUE
    }else{
      print(paste("Remove",colnames(new.X)[elim.index],"..."))
      new.X <- new.X[,-elim.index]
      if(ncol(new.X) == p.keep){
        select.terminate <- TRUE
      }
    }
  }
  fit.final <- glm(Y~new.X-1,family = "binomial",control = list(maxit = 100))
  return(list(model=fit.final,X=new.X,Y=Y))
}

logistic.back.sch <- function(X,Y,pval.thr=0.1,adj=TRUE,keep.index){
  n <- nrow(X)
  p <- ncol(X)
  X.keep <- X[,keep.index]
  X.select <- X[,-keep.index]
  p.keep <- length(keep.index)
  new.X <- cbind(X.keep,X.select)
  select.terminate <- FALSE
  while(!select.terminate){
    p.current <- ncol(new.X)
    select.index <- setdiff(1:p.current,1:p.keep)
    
    fit.old <- glm(Y~new.X-1,family = "binomial",control = list(maxit = 100))
    ss <- summary(fit.old)
    orig.p.val <- ss$coefficients[select.index,4]
    if(adj){
      select.p.val <- p.adjust(orig.p.val,method ="bonferroni")
    }else{
      select.p.val <- orig.p.val
    }
    elim.index <- select.index[which.max(select.p.val)]
    if(max(select.p.val)<=pval.thr){
      select.terminate <- TRUE
    }else{
      print(paste("Remove",colnames(new.X)[elim.index],"..."))
      new.X <- new.X[,-elim.index]
      if(ncol(new.X) == p.keep){
        select.terminate <- TRUE
      }
    }
  }
  fit.final <- glm(Y~new.X-1,family = "binomial",control = list(maxit = 100))
  return(list(model=fit.final,X=new.X,Y=Y))
}


SP.back <- function(X,Y,A,K,thr = 0.95,pval.thr=0.1,keep.index,adj=TRUE){
  n <- nrow(X)
  p <- ncol(X)
  X.keep <- X[,keep.index]
  X.select <- X[,-keep.index]
  p.keep <- length(keep.index)
  new.X <- cbind(X.keep,X.select)
  select.terminate <- FALSE
  while(!select.terminate){
    p.current <- ncol(new.X)
    select.index <- setdiff(1:p.current,1:p.keep)
    fit.old <- SP.Inf(X=new.X,Y=Y,A=A,K=K,thr = thr)
    orig.p.val <- fit.old$coef.mat[select.index,4]
    if(adj){
      select.p.val <- p.adjust(orig.p.val,method ="bonferroni")
    }else{
      select.p.val <- orig.p.val
    }
    elim.index <- select.index[which.max(select.p.val)]
    if(max(select.p.val)<=pval.thr){
      select.terminate <- TRUE
    }else{
      print(paste("Remove",colnames(new.X)[elim.index],"..."))
      new.X <- new.X[,-elim.index]
      if(ncol(new.X) == p.keep){
        select.terminate <- TRUE
      }
    }
  }
  fit.final <- SP.Inf(X=new.X,Y=Y,A=A,K=K,thr = thr)
  return(list(model=fit.final,X=new.X,Y=Y,A=A,K=K,thr=thr))
}

USVT.orig <- function (A)
{
  n <- nrow(A)
  K <- ceiling(n^(1/2))
  SVD <- irlba(A, nv = K, nu = K)
  Km <- sum(SVD$d>(2.01*sqrt(n*max(rowMeans(A)))))
  print(Km)
  Ahat <- SVD$u[,1:Km] %*% (t(SVD$v[,1:Km]) * SVD$d[1:Km])
  Ahat[Ahat > 1] <- 1
  Ahat[Ahat < 0] <- 0
  return(list(Ahat=Ahat,Km=Km))
}

USVT.embed <- function (A)
{
  n <- nrow(A)
  K <- ceiling(n^(1/2))
  eigens <- eigs(A, k = K)
  Km <- sum(eigens$values^2>(2.01*sqrt(n*max(rowMeans(A)))))
  print(Km)
  return(list(Km=Km))
}

