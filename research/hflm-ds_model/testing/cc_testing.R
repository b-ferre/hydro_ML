rm(list = ls())

source("./testing_funs.R")

#############   single catchment beta estimation func (stolen from Joe's code)   ##########################

slapdash_test <- function(gridcode) {
    #### get data
    data_dir_path <- "./data/"

    simu_x<-readRDS(file=paste0(data_dir_path, "grid",gridcode,"_rain.rds"))
    Y<-readRDS(file = paste0(data_dir_path, "grid",gridcode,"_streamflow.rds"))

    ## interpolate any missing values
    ## TODO: do this in a smarter fashion
    for (day in 1:365) {
        simu_x[is.na(simu_x[,day]), day] = mean(as.numeric(simu_x[,day]), na.rm = TRUE)
        Y[is.na(Y[,day]), day] = mean(as.numeric(Y[,day]), na.rm = TRUE)
    }

    M <- ncol(simu_x)-1
    max_lag<-150
    K <- (M+1)*max_lag
    n_hp<-65
    ars=2

    #regularization
    RegMat<-get_D(M,max_lag)
    hProd<-crossprod(RegMat[[1]])
    vProd<-crossprod(RegMat[[2]])

    x<-as.numeric(t(simu_x))
    y<-as.numeric(t(Y))
    y<-y[max_lag:length(y)]

    trainfrac<-0.6
    valfrac<-0.2
    testfrac<-0.2

    y_train<-y[1:round(length(y)*trainfrac)]
    y_val<-y[(round(length(y)*trainfrac)+1):round(length(y)*(trainfrac+valfrac))]
    y_test<-y[(round(length(y)*(trainfrac+valfrac))+1):length(y)]
    y_full<-y[1:round(length(y)*(trainfrac+valfrac))]

    x_train<-x[1:(max_lag+length(y_train)-1)]
    x_val<-x[(length(y_train)+1):(max_lag+length(y_train)+length(y_val)-1)]
    x_test<-x[(length(y_train)+length(y_val)+1):length(x)]
    x_full<-x[1:(max_lag+length(y_train)+length(y_val)-1)]

    first_day<-max_lag
    first_day_train<-max_lag
    first_day_val<-(length(y_train)+max_lag) %% 365
    first_day_test<-(length(y_train)+length(y_val)+max_lag) %% 365
    first_day_full<-max_lag

    Xmat_train<-get_x(x_train,M,max_lag,first_day_train)
    Xmat_val<-get_x(x_val,M,max_lag,first_day = first_day_val)
    Xmat_test<-get_x(x_test,M,max_lag,first_day = first_day_test)
    Xmat_full<-get_x(x_full,M,max_lag,first_day = first_day_full)
    Xmat<-get_x(x,M,max_lag,first_day = first_day)

    XXProd_train<-crossprod(Xmat_train)
    XXProd_val<-crossprod(Xmat_val)
    XXProd_test<-crossprod(Xmat_test)
    XXProd_full<-crossprod(Xmat_full)
    XXProd<-crossprod(Xmat)

    XyProd_train<-crossprod(Xmat_train,matrix(y_train,ncol = 1))
    XyProd_val<-crossprod(Xmat_val,matrix(y_val,ncol = 1))
    XyProd_full<-crossprod(Xmat_full,matrix(y_full,ncol = 1))
    XyProd_test<-crossprod(Xmat_test,matrix(y_test,ncol = 1))
    XyProd<-crossprod(Xmat,matrix(y,ncol = 1))
    ####

    To_Zero<-rep(F,K)
    result<-optimize_hp(hp_range=list(c(10,20),c(-5,15)),XTX=XXProd_train,hTh=hProd,vTv=vProd,XTy=XyProd_train,Xval=Xmat_val
                        ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)


    best<-which.max(result[[3]])
    w_h<-exp(result[[1]])
    w_v<-exp(result[[2]])

    pred_beta <- as.numeric(Matrix::solve(XXProd_full+ w_h[best]*hProd + w_v[best]*vProd, XyProd_full, sparse = TRUE))

    #detect where low betas are
    grpnorms<-get_grp_normsH(pred_beta,M,max_lag)
    # plot_bvec_better(grpnorms,M,max_lag)

    #optimize quantile for thresholding
    q_result<-optimize_q(grpnorms = grpnorms,Nq=500,Xmat_train = Xmat_full,Xmat_val = Xmat_full,
                        y_train = y_full,y_val = y_full,hMat = RegMat[[1]],vMat = RegMat[[2]],w_h_best = w_h[best],w_v_best = w_v[best])
    q_ngn<-q_result[[1]]
    R2<-q_result[[2]]
    distances1<-get_elbow(q_ngn,R2)
    whichbestq<-which.max(distances1)
    finalq<-q_ngn[whichbestq]
    To_Zero<-grpnorms<quantile(grpnorms,finalq)


    big_indexmat<-get_index_mat(M=M,max_lag = ars+ max_lag)
    small_indexmat<-get_index_mat(M=M,max_lag = max_lag)
    To_Zero_big<-rep(F,nrow(big_indexmat))
    for(m in 1:(M+1)){
    tot_in_col<-sum(To_Zero[small_indexmat$X_coord==(m-1)])
    To_Zero_big[big_indexmat$Y_coord> max(big_indexmat$Y_coord)-tot_in_col & big_indexmat$X_coord==(m-1)]<- T
    }
    To_Zero<-c(To_Zero_big,rep(F,ars))

    Xmat<-get_x_kir(x,y,M,max_lag,first_day,ars = ars)
    Xmat_train<- Xmat[1:(nrow(Xmat_train)),]
    Xmat_val<-Xmat[(nrow(Xmat_train)+1):nrow(Xmat_full),]
    Xmat_full<-rbind2(Xmat_train,Xmat_val)

    y_train<-y[(ars+1):(nrow(Xmat_train)+ars)]
    y_val<-y[(nrow(Xmat_train)+ars+1):(nrow(Xmat_full)+ars)]
    y_full<-c(y_train,y_val)
    y<-y[(ars+1):length(y)]

    RegMat<-get_D_kir(M,max_lag,est_arcoefs = rep(0,ars))
    XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
    XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
    hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
    vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])

    result<-optimize_hp(hp_range=list(c(10,20),c(-5,15)),XTX=XXProd_sparse,hTh=hProd_sparse,vTv=vProd_sparse,XTy=XyProd_sparse,Xval=Xmat_val
                        ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)

    best<-which.max(result[[3]])
    w_h<-exp(result[[1]])
    w_v<-exp(result[[2]])
    #refit model and update beta
    XyProd_sparse<-crossprod(Xmat_full[,!To_Zero],matrix(y_full,ncol = 1))
    XXProd_sparse<-crossprod(Xmat_full[,!To_Zero])

    pred_beta<-rep(0,ncol(Xmat))
    pred_beta[!To_Zero]<-as.numeric(Matrix::solve(XXProd_sparse+ w_h[best]*hProd_sparse + w_v[best]*vProd_sparse, XyProd_sparse, sparse = TRUE))
    est_delta<-compute_delta(pred_beta[1:(length(pred_beta)-ars)],M,max_lag = max_lag+ars)-ars

    indexmat_new<-get_index_mat(M,max_lag)
    pred_beta<-get_truebeta(pred_beta,indexmat_new,ars=ars)
    pred_beta<-Zero_out_pred_beta(b=pred_beta,d=est_delta)

    Y_hat_test=Xmat_test %*% pred_beta
    R2 <- 1-sum((y_test- Y_hat_test)^2)/sum((y_test- mean(y_full))^2)
    return(R2)
}

###################   run single catchment func on nth thru (n + 5)th catchment   #########################

args = commandArgs(trailingOnly=TRUE)
load("./index.rda")

i0 <- as.integer(args[1])
for (i in i0:(i0 + 4)) {
  gridcode <- indexed_gridcodes[[i]]

  print(paste0("RUNNING SLAPDASH TEST ON CATCHMENT NUMBER ", gridcode))
  try({
    r2 <- slapdash_test(gridcode)
    save(r2, file = paste0("results/r2_scores/grid", gridcode, "_r2.rda"))
  }, silent = FALSE)
}
