### run this to load all needed helpers for testing

#import functions
library(glmnet)
library(Matrix)
library(ggplot2)
library(lubridate)
library(caTools)
library(forecast)
library(zoo)
library(ranger)
library(missRanger)
library(GPfit)
library(lattice)
library(RColorBrewer)
library(segmented)

######################################   model code   ##################################################### nolint

#get order and location of all coefficient function nodes
get_index_mat <- function(M = 4,max_lag=3) {
  K <- (M+1)*max_lag # total number of nodes
  what<-expand.grid(0:(max_lag-1),0:M)
  data.frame("Node" = 1:K, "Y_coord" = what$Var1, "X_coord" = what$Var2)
}

#build smoothness regularization matrices
get_D <- function(M,max_lag) {
  ind_mat <- get_index_mat(M = M,max_lag = max_lag)
  K <- nrow(ind_mat)
  num_pairs <- (M+1)*(max_lag-1)
  
  # First we build the D_h matrix
  D_v <- Matrix(0, nrow = num_pairs, ncol = K)
  ind_diff_v <- which(diff(ind_mat[,3]) == 0)
  
  D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v])] <- -1
  D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v]+1)] <- 1
  D_v<-rbind(D_v,sparseMatrix(i=1:(M+1),j=which(ind_mat$Y_coord == max(ind_mat$Y_coord)), x=rep(1,(M+1)),dims = c(M+1,K) )) # add vanishing boundary to top
  
  # Secondly we have the D_v matrix
  num_pairs <- (M)*(max_lag)
  D_h <- Matrix(0, nrow = num_pairs, ncol = K)
  ind_mat2 <- ind_mat[order(ind_mat[,2]),]
  ind_diff_h<- which(diff(ind_mat2[,2]) == 0)
  D_h[cbind(1:num_pairs, sort(ind_mat2[ind_diff_h,1]))] <- -1
  D_h[cbind(1:num_pairs,sort(ind_mat2[ind_diff_h+1,1]))] <- 1
  D_h<-rbind(D_h,sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == min(ind_mat$X_coord)), x=rep(1,max_lag),dims = c(max_lag,K) ) + sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == max(ind_mat$X_coord)), x=rep(-1,max_lag),dims = c(max_lag,K) )) # add periodicity condition
  
  return(list(D_h,D_v))
}

get_grp_normsH <- function(bvec, M,max_lag) {
  grps <- as.factor(rep(0:M, rep(max_lag,M+1)))
  norms <- tapply(bvec, grps, FUN = function(x){
    cumsum(as.numeric(x[length(x):1])^2)[length(x):1]
  })
  
  as.numeric(do.call("c", norms))
}

#plotting beta(s,t)
plot_bvec_better<-function(bvec,M,max_lag){
  index_dat<-get_index_mat(M,max_lag)
  delta<-compute_delta(bvec,M,max_lag = max_lag)
  
  mat<-matrix(0,ncol = M+1,nrow = max_lag)
  mat[cbind(index_dat$Y_coord+1,index_dat$X_coord+1)]<-bvec
  mat<-matrix(mat[1:min(ceiling(max(delta)/10)*10,max_lag),],ncol = M+1)
  range_p<-max(abs(mat[mat>0]),2e-15)
  range_n<-max(abs(mat[mat<0]),2e-15)
  max_range<-max(abs(mat))
  tot_range<-range_p+ range_n
  totalcols<-35
  rangepercol<-tot_range/totalcols
  numblues<-ceiling(range_p/rangepercol)
  numreds<-totalcols- numblues +1
  
  coul <- colorRampPalette(brewer.pal(8, "Blues"))(max(numblues+2,ceiling(totalcols/2)))
  coul2 <- colorRampPalette(brewer.pal(8, "Reds"))(max(numreds+2,ceiling(totalcols/2)))
  colorvec<-c(coul2[(numreds-1+2):2],"#FFFFFF",coul[3:(numblues-1+3)])
  
  rownames(mat)<-0:(nrow(mat)-1)
  mat<-setNames(reshape2::melt(t(mat)),c("x","y","z"))
  
  levelplot(z~x*y,mat,col.regions = colorvec,at=sort(c(seq(from=-1e-15,to=-range_n,length.out=numreds+1 ),seq(from=1e-15,to=range_p,length.out=numblues+1))),
            ylab=list("Lag (Days)",cex=1.4),xlim = c(0.5,M+1+0.5),aspect = "xy",xlab=list("Day of the year" ,cex=1.4),colorkey=list(axis.text=list(cex=1.2)),
            scales=list(x=list(cex=1.2),y=list(cex=1.2)))
}

# form matrix of input var
get_x<-function(x_ts,M,max_lag,first_day){
  K=(M+1)*max_lag
  nrowX<-(length(x_ts)-max_lag+1)
  i_vec<-rep(0, max_lag*nrowX)
  j_vec<-rep(0, max_lag*nrowX)
  x_vec<-rep(0, max_lag*nrowX)
  
  indexCounter<-0
  for(i in 1:nrowX){
    i_vec[(indexCounter+1):(indexCounter+max_lag)]<-rep(i,max_lag)
    j_vec[(indexCounter+1):(indexCounter+max_lag)]<-(((indexCounter + (first_day-1)*max_lag) %% K):((indexCounter+max_lag-1 + (first_day-1)*max_lag) %% K))+1
    x_vec[(indexCounter+1):(indexCounter+max_lag)]<-x_ts[(i+max_lag-1):(i)]
    indexCounter<-indexCounter+max_lag
  }
  
  sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(nrowX, K))
}

compute_delta<-function(bvec,M,max_lag){
  deltavec<-rep(0,M+1)
  
  index_dat<-get_index_mat(M,max_lag)
  
  for(i in 0:M){
    curobs<-bvec[index_dat$X_coord==i]
    curobs<-curobs[length(curobs):1]
    deltavec[i+1]<-max_lag-max(which(cumsum(curobs)==0),0)
  }
  deltavec
}

#optimize w_h and w_v hyperparameters with bayesian optimization
optimize_hp<-function(hp_range,XTX,hTh,vTv,XTy,Xval,ytrain,yval,n_hp,To_Zero){
  n_init<-30 #get starting points
  R2<-rep(0,round(n_init))
  
  max_w_h<-hp_range[[1]][2]
  max_w_v<-hp_range[[2]][2]
  
  min_w_h<-hp_range[[1]][1]
  min_w_v<-hp_range[[2]][1]
  
  w_h<-runif(length(R2),min=min_w_h,max=max_w_h)
  w_v<-runif(length(R2),min=min_w_v,max=max_w_v)
  
  for(h in 1:length(R2)){
    pred_beta<-rep(0,ncol(Xval))
    pred_beta[!To_Zero] <- Matrix::solve(XTX+ exp(w_h[h])*hTh + exp(w_v[h])*vTv, XTy, sparse = TRUE)
    Y_hat_val=Xval %*% pred_beta
    R2[h]<-1-sum((yval- Y_hat_val)^2)/sum((mean(yval)-yval)^2)
    R2[h]<-max(R2[h],-1)
    print(paste0("This is R2[h]: ",R2[h], " , and h is: ",h," and w_h is: ",w_h[h], " and w_v is: ", w_v[h]))
  }
  
  while(h<n_hp){ #start gaussian process
    
    w_h_scaled<-(w_h- min_w_h)/(max_w_h- min_w_h)
    w_v_scaled<-(w_v- min_w_v)/(max_w_v- min_w_v)

    mod<-GP_fit(X=data.frame(w_h=w_h_scaled, w_v=w_v_scaled),Y=R2)
    
    w_h_long<-runif(1000,min=min_w_h,max=max_w_h)
    w_v_long<-runif(1000,min=min_w_v,max=max_w_v)
    
    w_h_long_scaled<-(w_h_long- min_w_h)/(max_w_h- min_w_h)
    w_v_long_scaled<-(w_v_long- min_w_v)/(max_w_v- min_w_v)
    
    pred <- predict.GP(mod, xnew = data.frame(w_h=w_h_long_scaled,w_v=w_v_long_scaled))
    mu <- pred$Y_hat
    sigma <- sqrt(pred$MSE)
    Z <- (mu - max(pred$Y_hat))/sigma
    expected_imp <- sigma*(Z  * pnorm(Z) + dnorm(Z))
    expected_imp[is.na(expected_imp)]<-0
    
    nextone<-which.max(expected_imp)
    w_h<-c(w_h,w_h_long[nextone])
    w_v<-c(w_v,w_v_long[nextone])
    
    h=h+1
    pred_beta<-rep(0,ncol(Xval))
    pred_beta[!To_Zero] <- Matrix::solve(XTX+ exp(w_h[h])*hTh + exp(w_v[h])*vTv, XTy, sparse = TRUE)
    Y_hat_val=Xval %*% pred_beta
    curR2<- 1-sum((yval- Y_hat_val)^2)/sum((mean(yval)-yval)^2)
    curR2<-max(curR2,-1)
    R2<-c(R2,curR2)
    print(paste0("Gaussian process: This is R2[h]: ",R2[h], " , and h is: ",h," and w_h is: ",w_h[h], " and w_v is: ", w_v[h]))
  }
  list(w_h,w_v,R2)
}

#find elbow point of monotone decreasing plot
get_elbow<-function(qVec,R2Vec){
  qVecScaled<-(qVec-min(qVec))/(max(qVec)-min(qVec))
  R2VecScaled<- (R2Vec-min(R2Vec))/(max(R2Vec)-min(R2Vec))
  kneedledist<-R2VecScaled - (1-qVecScaled) 
  kneedledist
}

compute_beta_R2<-function(est,true){
  1-sum((est-true)^2)/sum((true-mean(true))^2)
}

optimize_q<-function(grpnorms,Nq=500,Xmat_train,Xmat_val,y_train,y_val,hMat,vMat,w_h_best,w_v_best){
  q_ngn<-seq(from=0,to=1-365/dim(Xmat_train)[2],length.out=100) #important to select this carefully
  R2<-rep(0,length(q_ngn))
  for(h in 1:length(R2)){
    To_Zero<-grpnorms<quantile(grpnorms,q_ngn[h])
    XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
    XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
    hProd_sparse<-crossprod(hMat[,!To_Zero])
    vProd_sparse<-crossprod(vMat[,!To_Zero])
    
    pred_beta<-rep(0,ncol(Xmat_val))
    pred_beta[!To_Zero]<- Matrix::solve(XXProd_sparse+ w_h_best*hProd_sparse + w_v_best*vProd_sparse, XyProd_sparse, sparse = TRUE)
    Y_hat_val=Xmat_val %*% pred_beta
    R2[h]<-1-sum((y_val- Y_hat_val)^2)/sum((y_val-mean(y_val))^2)
    print(paste0("This is sparse R2[h]: ",R2[h], " , and h is: ",h))
  }
  
  whichbestq<-which.max(get_elbow(q_ngn,R2))
  new_q<-q_ngn[q_ngn<= q_ngn[whichbestq]]
  new_R2<-R2[q_ngn<= q_ngn[whichbestq]]
  
  whichbestq2<-which.max(get_elbow(new_q,new_R2))
  finalq2<- new_q[whichbestq2]
  finalq1<- q_ngn[whichbestq]
  searchRange<-c(max(0,finalq2-0.05),min(1,finalq1+0.05))
  
  q_ngn2<-runif(Nq- length(R2),min=searchRange[1],max=searchRange[2]) #important to select this carefully
  R22<-rep(0,length(q_ngn2))
  for(h in 1:length(R22)){
    To_Zero<-grpnorms<quantile(grpnorms,q_ngn2[h])

    XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
    XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
    hProd_sparse<-crossprod(hMat[,!To_Zero])
    vProd_sparse<-crossprod(vMat[,!To_Zero])
    
    pred_beta<-rep(0,length(pred_beta))
    pred_beta[!To_Zero]<- Matrix::solve(XXProd_sparse+ w_h_best*hProd_sparse + w_v_best*vProd_sparse, XyProd_sparse, sparse = TRUE)
    Y_hat_val=Xmat_val %*% pred_beta
    R22[h]<-1-sum((y_val- Y_hat_val)^2)/sum((y_val-mean(y_val))^2)
    if(h%%50==0) print(h/length(R22))
  }
  q_all<-c(q_ngn,q_ngn2)
  R2_all<-c(R2,R22)
  
  return(list(q_all,R2_all))
}

get_x_kir<-function(x_ts,y_ts,M,max_lag,first_day,ars){
  differ<-length(x_ts)- length(y_ts)
  max_lag=max_lag+ars
  first_day= (first_day + ars) %% (M+1)
  nrowX<- length(x_ts)-max(max_lag-1,ars)
  
  K=(M+1)*(max_lag)

  i_vec<-rep(0, max_lag*nrowX)
  j_vec<-rep(0, max_lag*nrowX)
  x_vec<-rep(0, max_lag*nrowX)
  
  indexCounter<-0
  for(i in 1:nrowX){
    i_vec[(indexCounter+1):(indexCounter+max_lag)]<-rep(i,max_lag)
    j_vec[(indexCounter+1):(indexCounter+max_lag)]<-(((indexCounter + first_day*max_lag) %% K):((indexCounter+max_lag-1 + first_day*max_lag) %% K))+1
    x_vec[(indexCounter+1):(indexCounter+max_lag)]<-x_ts[(i+max_lag-1):i]
    indexCounter<-indexCounter+max_lag
  }
  
  AllX=sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(nrowX, K))
  
  i_vec<-rep(0, ars*nrowX)
  j_vec<-rep(0, ars*nrowX)
  x_vec<-rep(0, ars*nrowX)
  
  indexCounter<-0
  for(i in 1:nrowX){
    i_vec[(indexCounter+1):(indexCounter+ars)]<-rep(i,ars)
    j_vec[(indexCounter+1):(indexCounter+ars)]<-1:ars
    x_vec[(indexCounter+1):(indexCounter+ars)]<-y_ts[(i+ars-1):(i)]
    indexCounter<-indexCounter+ars
  }
  
  AllY=sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(nrowX, ars))
  
  Allboth<-cbind2(AllX,AllY)
  Allboth
}

get_truebeta<-function(pred_beta,indexmat,ars){
  y_cords<-unique(indexmat$Y_coord)
  big_indexmat<-get_index_mat(M=max(indexmat$X_coord),max_lag = max(y_cords+1)+ars)
  big_indexmat$val<-pred_beta[1:(length(pred_beta)-ars)]
  indexmat<-merge(indexmat,big_indexmat,by=c("X_coord","Y_coord"),all.x = T)
  indexmat<-indexmat[order(indexmat$Node.x),]
  arcs<-pred_beta[(length(pred_beta)-ars+1):length(pred_beta)]
  real_beta<-indexmat$val
  
  for(y_cord in y_cords[2:length(y_cords)]){
    aligned_beta<-real_beta[which(indexmat$Y_coord==y_cord-1)]
    aligned_beta<-aligned_beta[c(length(aligned_beta),1:(length(aligned_beta)-1))]
    real_beta[indexmat$Y_coord==y_cord]<-real_beta[indexmat$Y_coord==y_cord]+ arcs[1]* aligned_beta
    if(y_cord>1 & length(arcs)>1){
      aligned_beta<-real_beta[which(indexmat$Y_coord==y_cord-2)]
      aligned_beta<-aligned_beta[c(length(aligned_beta)-1,length(aligned_beta),1:(length(aligned_beta)-2))]
      real_beta[indexmat$Y_coord==y_cord]<-real_beta[indexmat$Y_coord==y_cord]+ arcs[2]* aligned_beta
    }
  }
  real_beta
}

#build smoothness regularization matrices
get_D_kir <- function(M,max_lag,est_arcoefs) {
  ars=length(est_arcoefs)
  max_lag=max_lag+ars
  ind_mat <- get_index_mat(M = M,max_lag = max_lag)
  K <- nrow(ind_mat)
  num_pairs <- (M+1)*(max_lag-1)
  
  D_v<-bandSparse(K,k=c(0,1),diagonals =list(rep(-1,K), rep(c(rep(1,max_lag-1),0),M+1)[-K]))
  #D_v<-bandSparse(K,k=c(-1,0,1),diagonals =list(rep(est_arcoefs[1]^2+est_arcoefs[2]-est_arcoefs[1],K-1),rep(est_arcoefs[1]-1,K), rep(c(rep(1,max_lag-1),0),M+1)[-K]))
  extracols<-Matrix(0,nrow = nrow(D_v),ncol = ars)
  D_v<-cbind2(D_v,extracols)
  extrarows<-Matrix(0,nrow = ars,ncol = ncol(D_v))
  D_v<-rbind2(D_v,extrarows)
  
  
  #D_v[1,1]<- est_arcoefs[1]-1
  #D_v[2,1]<- est_arcoefs[1]^2+ est_arcoefs[2]- est_arcoefs[1]
  #D_v[2,2]<- est_arcoefs[1]-1
  
  # max_lag=max_lag+ars
  # ind_mat <- get_index_mat(M = M,max_lag = max_lag)
  # K <- nrow(ind_mat)
  # Secondly we have the D_v matrix
  if(M>0){
    num_pairs <- (M)*(max_lag)
    D_h <- Matrix(0, nrow = num_pairs, ncol = K)
    ind_mat2 <- ind_mat[order(ind_mat[,2]),]
    ind_diff_h<- which(diff(ind_mat2[,2]) == 0)
    D_h[cbind(1:num_pairs, sort(ind_mat2[ind_diff_h,1]))] <- -1
    D_h[cbind(1:num_pairs,sort(ind_mat2[ind_diff_h+1,1]))] <- 1
    D_h<-rbind(D_h,sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == min(ind_mat$X_coord)), x=rep(1,max_lag),dims = c(max_lag,K) ) + sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == max(ind_mat$X_coord)), x=rep(-1,max_lag),dims = c(max_lag,K) )) # add periodicity condition
  } else{
    D_h= sparseMatrix(i=1,j=1,x=0,dims = c(1,K))
  }
  extracols<-Matrix(0,nrow = nrow(D_h),ncol = ars)
  D_h<-cbind2(D_h,extracols)
  extrarows<-Matrix(0,nrow = ars,ncol = ncol(D_h))
  D_h<-rbind2(D_h,extrarows)
  
  return(list(D_h,D_v))
  
}

Zero_out_pred_beta<- function(b,d){
  newb<-b
  M<-length(d)-1
  max_lag<-length(b)/(M+1)
  ind_mat<-get_index_mat(M,max_lag)
  for(m in 1:(M+1)){
    newb[ind_mat$Y_coord>(d[m]-1) & ind_mat$X_coord==(m-1)]<-0
  }
  newb
}

