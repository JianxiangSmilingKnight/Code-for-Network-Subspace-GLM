source("case_study_function.R")
registerDoParallel(cores=10)

load("averge_network_embed.Rda")
load("Y.Rda")
load("X_true_embed.Rda")
load("X_Naive_embed.Rda")
load("X_SP_embed.Rda")
load("cv_order.Rda")


h=200
n<-nrow(X.true)
cv_order_index<-1:n
A_new<-G

###Our method
usvt.reg <- USVT.orig(A_new)
Y<-as.matrix(Y_D)
X<-X_SP

prediction<-matrix(0,nrow = n,ncol = 1)
prediction.result <-for (II in 1:h) {
  index <- c(cv_order_index[cv_order!=II],cv_order_index[cv_order==II])
  M <- sum(cv_order==II)
  #index <- sample(1:n,n,replace=FALSE)
  
  perm.X <- X[index,]
  perm.Y <- Y[index,]
  perm.A.orig <- A_new[index,index]
  tmp.result <- NULL
  SP.perm <- SP.semi.Inf(X=perm.X,Y = matrix(perm.Y[1:(n-M)],ncol=1),A=perm.A.orig,K = usvt.reg$Km,r=0,thr=0.95)
  prediction[cv_order==II,]=SP.perm$full.fitted
}

A_new<-G
Y<-as.matrix(Y_D)
X<-X_Naive

###Naive Logistic Reg
prediction_NL<-matrix(0,nrow = n,ncol = 1)
for (II in 1:h) {
  index <- c(cv_order_index[cv_order!=II],cv_order_index[cv_order==II])
  M <- sum(cv_order==II)
  
  perm.X <- X[index,]
  perm.Y <- Y[index,]
  perm.A.orig <- A_new[index,index]
  tmp.result <- NULL
  NL.perm <- glm(matrix(perm.Y[1:(n-M)],ncol=1)~0+perm.X[1:(n-M),],family = "binomial",control = list(maxit = 100))
  prediction_NL[cv_order==II,]<-perm.X[(n-M+1):n,]%*%NL.perm$coefficients
  print(II)
  
}
prediction_NL<-exp(prediction_NL)/(1+exp(prediction_NL))


Y<-as.matrix(Y_D)
X<-X.true[,2:12]
n=nrow(X.true)
prediction_RNC<-matrix(0,nrow = n,ncol = 1)

###RNC, with glmnet()

predictionRNC<-foreach(II=1:h,.packages="RSpectra",.combine = cbind) %dopar% {
  index <- c(cv_order_index[cv_order!=II],cv_order_index[cv_order==II])
  M <- sum(cv_order==II)
  #index <- sample(1:n,n,replace=FALSE)
  
  perm.X <- X[index,]
  perm.Y <- Y[index,]
  perm.A.orig <- A_new[index,index]
  perm.A.sub <- perm.A.orig[1:(n-M),1:(n-M)]
  gamma<-0.05
  perm.X.sub<-perm.X[1:(n-M),]
  D <- diag(rowSums(perm.A.sub))
  L_gamma <- D - perm.A.sub + diag(rep(gamma,n-M))
  eigen_decomp <- eigen(L_gamma, symmetric = TRUE)
  Xtmp<-cbind(perm.X.sub,(eigen_decomp$vectors))
  penalty_factors <- c(rep(0, 11), eigen_decomp$values)
  cv_fit_perm <- cv.glmnet(Xtmp, perm.Y[1:(n-M)], family = "binomial",nfolds = 5, penalty.factor = penalty_factors, alpha = 0,lambda =  exp(seq(log(1000),log(0.001),length.out=10)))
  optimal_lambda_perm <- cv_fit_perm$lambda.min
  fit_perm <- glmnet(Xtmp, perm.Y[1:(n-M)], family = "binomial", alpha = 0,intercept = FALSE, penalty.factor =penalty_factors,lambda = optimal_lambda_perm)
  alpha<-(eigen_decomp$vectors)%*%as.matrix(fit_perm$beta[(ncol(X)+1):(ncol(Xtmp))])
  Xbeta<-perm.X[(n-M+1):n,]%*%as.matrix(fit_perm$beta[1:ncol(X)])
  D <- diag(rowSums(perm.A.orig))
  L <- D - perm.A.orig + diag(rep(gamma,n))
  L11 <- Matrix(L[(n-M+1):n,(n-M+1):n],sparse=TRUE)
  L12 <- -1*L[(n-M+1):n,1:(n-M)]
  alpha_pred <- solve(L11,L12%*%alpha,sparse=TRUE)
  c(exp(as.matrix(alpha_pred+Xbeta))/(1+exp(as.matrix(alpha_pred+Xbeta))))
}
prediction_RNC<-rep(0,n)
for (II in 1:h) {
  prediction_RNC[cv_order==II]=predictionRNC[,II]
}

save(prediction,file="embed_prediction.Rda")
save(prediction_NL,file="embed_prediction_NL.Rda")
save(prediction_RNC,file="embed_prediction_RNC.Rda")