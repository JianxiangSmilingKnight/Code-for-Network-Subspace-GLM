source("case_study_function.R")
registerDoParallel(cores=10)

load("average_network_embed.Rda")
load("Y.Rda")
load("X_true_embed.Rda")
load("cv_order.Rda")


h=200
n<-nrow(X.true)
cv_order_index<-1:n
A_new<-G

###Our method
usvt.reg <- USVT.embed(A_new)
Y<-as.matrix(Y)
X<-X.true

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

prediction_embed<-prediction


save(prediction_embed,file="predictionN2V.Rda")