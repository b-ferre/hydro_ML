
rm(list=ls())  #clear environment
setwd("~/IRF paper")

library(data.table)
library(dplyr)

source("IRFnnhs.R")



# alpha=1 test against benchmark (Figure 2, left column)
set.seed(99)
alpha <- 1.0          # gamma shape factor: compare values of 1, 2, and 4 to see contrasting effects
tau <- 5             # mean of gamma distribution: compare values of 100, 50, and 20 to see contrasting effects

# # alpha=2 test against benchmark (Figure 2, middle column)
# set.seed(299)
# alpha <- 2.0          # gamma shape factor: compare values of 1, 2, and 4 to see contrasting effects
# tau <- 10             # mean of gamma distribution: compare values of 100, 50, and 20 to see contrasting effects

# # alpha=4 test against benchmark (Figure 2, right column)
# set.seed(599)
# alpha <- 4.0          # gamma shape factor: compare values of 1, 2, and 4 to see contrasting effects
# tau <- 20             # mean of gamma distribution: compare values of 100, 50, and 20 to see contrasting effects


n <- 1e4              # length of time series
m <- 60               # longest lag to be analyzed
nu <- 0.0             # regularization parameter

y_err_scale <- 1.0    # RMS amplitude of y errors, as a fraction of the standard deviation of the true y values
# y_err_model <- list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488))  # AR and MA coefficients of the ARMA model that will generate the noise
y_err_model <- list(ar = c(0.95), ma=c(-0.2, 0.2))  # AR and MA coefficients of the ARMA model that will generate the noise

x_rho <- 0.5         # first-order serial correlation coefficient of input time series

n_iter <- 100

t <- seq(n)                                              # time axis
lags <- 0:m                                              # lag time series


true_kernel <- dgamma(t-1, shape=alpha, scale=tau/alpha)   # true convolution kernel is gamma distribution

benchmk <- true_kernel[1:(m+1)]                          # this is the benchmark that we should be able to reproduce


RMSE <- rep(0, m+1)                                      # will be used to compute the root-mean-square error of the IRF relative to the benchmark
pooled_se <- rep(0, m+1)                                 # will be used to compute the pooled standard error of the IRF estimates
bias <- rep(0, m+1)                                      # will be used to compute any average bias in the IRF estimates

fileID <- paste0("alpha=", alpha, "_tau=", tau)

for (iter in 1:n_iter) {
  
  if (x_rho==0) x <- rnorm(n)
  else x <- as.vector(unname(arima.sim(n=n, model=list(ar=x_rho))))            # true input
  
  y_true <- convolve(x, true_kernel, conj=FALSE)           # true output
  
  y_err <- as.vector(unname(arima.sim(n=n, model=y_err_model)))      # error is ARMA noise
  y_err <- y_err *(y_err_scale * sd(y_true)/sd(y_err))
  
  y <- y_true + y_err                                      # measured output (true plus error)
  
  
  
  zz <- IRF(y, x, m=m, nu=nu)
  
  bb <- as.numeric(zz$IRF)
  se <- as.numeric(zz$se)
  
  
  RMSE <- RMSE + (bb - benchmk)^2
  bias <- bias + (bb - benchmk)
  pooled_se <- pooled_se + se^2
  
  if (iter==n_iter) cat("iteration", iter, "of", n_iter, "done        \n")
    else cat("iteration", iter, "of", n_iter, "done        \r")

} # next iter

RMSE <- sqrt(RMSE/n_iter)
bias <- bias/n_iter
pooled_se <- sqrt(pooled_se/n_iter)
bias_se <- pooled_se/sqrt(n_iter)                         # uncertainty in average IRF, and thus uncertainty in bias averaged over n_iter trials


# now plot time series

ysnip <- y[1001:2000]
ytruesnip <- y_true[1001:2000]
xsnip <- x[1001:2000]
tsnip <- 1001:2000

maxval <- max(xsnip, ysnip, ytruesnip, na.rm=TRUE) #generate upper limit for plotting
minval <- min(xsnip, ysnip, ytruesnip, na.rm=TRUE) #generate upper limit for plotting

# here is x, in gray
plot(1001:2000, xsnip, main="x, y, and ytrue", xlab="t", ylab="x and y", xlim=c(min(tsnip), max(tsnip)), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20, cex.main=0.6)

# here is y, in green
par(new=TRUE)
plot(1001:2000, ysnip,  xlab="", ylab="", xlim=c(min(tsnip), max(tsnip)), ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20)

# here is y, in blue
par(new=TRUE)
plot(1001:2000, ytruesnip,  xlab="", ylab="", xlim=c(min(tsnip), max(tsnip)), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20)



# now plot time series

ysnip <- y[1001:2000]
ytruesnip <- y_true[1001:2000]
xsnip <- x[1001:2000]
tsnip <- 1001:2000

maxval <- max(ysnip, ytruesnip, na.rm=TRUE) #generate upper limit for plotting
minval <- min(ysnip, ytruesnip, na.rm=TRUE) #generate upper limit for plotting

# here is y, in green
plot(1001:2000, ysnip, main="y and ytrue", xlab="t", ylab="y and y_true", xlim=c(min(tsnip), max(tsnip)), ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20, cex.main=0.6)

# here is y, in blue
par(new=TRUE)
plot(1001:2000, ytruesnip,  xlab="", ylab="", xlim=c(min(tsnip), max(tsnip)), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20)






# now plot results

maxval <- max(benchmk, bb, na.rm=TRUE) #generate upper limit for plotting
maxval <- 1.05*maxval  #and extend this a bit
minval <- -0.1*maxval  #use standardized scaling so plots can be compared

# here is the benchmark, in blue
plot(lags, benchmk, main="compare benchmark to IRF", xlab="lag", ylab="kernel", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20, cex.main=0.6)

# here is IRF, in green
par(new=TRUE)
plot(lags, bb, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20)

# and here are the uncertainty bounds (+/- 2 standard errors), in gray
par(new=TRUE)
plot(lags, bb+2*se, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20)

par(new=TRUE)
plot(lags, bb-2*se, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20)







# now plot results (RMSE vs pooled standard error)

maxval <- max(RMSE, pooled_se, na.rm=TRUE) #generate upper limit for plotting
maxval <- 1.05*maxval  #and extend this a bit
minval <- -0.1*maxval  #use standardized scaling so plots can be compared

# here is RMSE, in blue
plot(lags, RMSE, main="compare RMSE to pooled SE", xlab="lag", ylab="RMSE and se", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20, cex.main=0.6)

par(new=TRUE)
plot(lags, se, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20)



# now plot results (RMSE vs pooled standard error)

maxval <- max(bias, -bias, 2*bias_se, na.rm=TRUE) #generate upper limit for plotting
maxval <- 1.05*maxval  #and extend this a bit
minval <- -1.05*maxval  #use standardized scaling so plots can be compared

# here is bias, in blue
plot(lags, bias, main="compare bias to uncertainty bounds", xlab="lag", ylab="bias and 2*se", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20, cex.main=0.6)

par(new=TRUE)
plot(lags, 2*bias_se, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="gray", pch=20)

par(new=TRUE)
plot(lags, -2*bias_se, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="gray", pch=20)




print(cbind(lags, benchmk, bb, se, RMSE, pooled_se, bias, bias_se))

fwrite(data.table(cbind(t=1:n, x, y_true, y_err, y)), file=paste0("benchmark_time_series_", fileID, ".txt"), row.names=FALSE, sep="\t")

fwrite(data.table(cbind(lags, benchmk, bb, se, RMSE, pooled_se, bias, bias_se)), file=paste0("benchmark_test_IRF_", fileID, ".txt"), row.names=FALSE, sep="\t")

