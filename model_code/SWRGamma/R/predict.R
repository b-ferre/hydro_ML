#' @describeIn predict_target Convolve the time series with the window kernel
#' @param kernel a window kernel
#' @noRd
convolve_window <- function(ts, kernel){

  # compute lengths of ts and kernel
  ts_length <- length(ts)
  kernel_length <- length(kernel)

  if(ts_length < kernel_length){
    warning("kernel larger than ts, returning NA")
    return(rep(NA, ts_length))
  }

  # initialize prediction with 0
  res <- rep(0, ts_length)

  # convolve kernel
  for(i in kernel_length : ts_length){
    res[i] <- sum(kernel * ts[(i - kernel_length + 1) : i])
  }

  return(res)
}

#' @describeIn predict_target Build a **GAMMA** window kernel
#' @importFrom stats pnorm
#' @noRd
#' @export
build_gamma_kernel <- function(param){
  if(any(is.na(param)) | param[3] <= 0 | param[2] <= 0){
    return(0)
  }

  delta <- param[1]
  shape <- param[2]
  rate <- param[3]

  # define kernel rad. such that (0, r] covers 97.5% of the probability
  # mass of gamma dist, as w og paper
  kernel_rad <- qgamma(0.975, shape = shape, rate = rate)

  #print(param)
  #print(kernel_rad)

  # define kernel indexes
  calc_kernel_ind <- c(0, ceiling(kernel_rad))  # unshifted for integration
  acc_kernel_ind <- c(-ceiling(kernel_rad + delta), -ceiling(delta))

  # plot shape of gamma probability density on unshifted axis
  #ts <- seq(calc_kernel_ind[1], calc_kernel_ind[2], by = 0.1)
  #ys <- dgamma(ts, shape, rate = rate)
  #plot(ts, ys)
  # then on shifted access
  #acc_ts <- seq(acc_kernel_ind[1], acc_kernel_ind[2], by = 0.1)
  #plot(acc_ts, rev(ys))

  # initialize kernel values with 0 and define (UNADJ.) bins to make calc work
  x <- seq(acc_kernel_ind[1], 0)
  bins <- seq(calc_kernel_ind[1] - 0.5, calc_kernel_ind[2] + 0.5)

  # compute kernel values for window by integrating over Gamma dist. in each bin
  kernel <- diff(pgamma(bins, shape, rate = rate))
  kernel <- c(kernel, rep(0, length(x) - length(kernel)))
  names(kernel) <- x

  # scale kernel to account for (necessary) truncation
  kernel <- kernel / sum(kernel, na.rm = TRUE)

  #print(param)
  #print(kernel)

  return(kernel)
}

#' @title Predict target variable
#' @description predicts the target variable given a time series of inputs, and trained parameters
#' @param ts a vector or ts object of new model inputs to predict
#' @inheritParams createSWR
#' @param ... currently unused
#' @noRd
predict_target <- function(ts, mix, param, log = FALSE, ...){
  if(!is.vector(ts)){
    stop("Error in predict: ts must be a vector")
  }
  if(!is.vector(mix) || !is.matrix(param) || length(mix) != nrow(param)){
    stop("Error in predict: provided parameters are not consistent")
  }

  # compute offset and convolution for each window
  if(nrow(param) > 0){
    conv <- apply(param,
                  1,
                  function(x, ts){
                    kernel <- build_gamma_kernel(x)
                    return(convolve_window(ts, kernel))
                  },
                  ts = ts)
    if(log){
      conv <- log(conv)
    }
    res <- as.vector(conv %*% mix)

  }
  else{
    res <- rep(0, length(ts))
  }
  return(res)
}

#' @title Predict target variable
#' @description Predicts the output time series given an input time series using a trained `SWR` model object.
#' @param object an `SWR` model object created using \link{trainSWR}
#' @param newdata a vector or ts object of new model inputs to predict
#' @param ... currently unused
#' @importFrom stats predict
#' @importFrom methods is
#' @export
predict.SWR <- function(object, newdata,...){

  if(!is(object, "SWR")){
    stop("Wrong class of object")
  }

  return(predict_target(ts = newdata,
                 param = object$param,
                 mix = object$mix,
                 ...))
}