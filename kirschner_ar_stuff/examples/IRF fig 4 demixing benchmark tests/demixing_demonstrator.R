setwd("~/IRF paper")    # change this to wherever the called scripts are


rm(list=ls())  #clear environment

library(data.table)
library(dplyr)


source("IRFnnhs.R")  # impulse response function script




n <- 1e4              # length of time series
m <- 40               # longest lag to be analyze


# # first benchmark test
# set.seed(456)
# alpha1 <- 2           # gamma shape factor for convolution 1
# alpha2 <- 6         # gamma shape factor for convolution 2
# alpha3 <- 12        # gamma shape factor for convolution 3
# tau1 <- 10            # mean response time for convolution 1
# tau2 <- 15            # mean response time for convolution 2
# tau3 <- 20            # mean response time for convolution 3
# 
# x_rho <- 0.0          # serial correlation coefficient for each of the input time series
# cross_rho <- 0.8      # correlation coefficient between input time series x1, x2, and x3
# 
# 
# y_err_scale <- 0.5    # RMS amplitude of y errors, as a fraction of the standard deviation of the true y values
# y_err_model <- list(ar=c(0.9), ma=c(-0.2, 0.2))
# fileID <- "demixing_benchmark_test1_Fig_4ace"


second benchmark test
set.seed(456)
alpha1 <- 2         # gamma shape factor for convolution 1
alpha2 <- 2        # gamma shape factor for convolution 2
alpha3 <- 2         # gamma shape factor for convolution 3
tau1 <- 5            # mean response time for convolution 1
tau2 <- 10            # mean response time for convolution 2
tau3 <- 20            # mean response time for convolution 3

x_rho <- 0.0           # serial correlation coefficient for each of the input time series
cross_rho <- 0.8            # correlation coefficient between input time series x1, x2, and x3


y_err_scale <- 0.5    # RMS amplitude of y errors, as a fraction of the standard deviation of the true y values
y_err_model <- list(ar=c(0.9), ma=c(-0.2, 0.2))
fileID <- "demixing_benchmark_test2_Fig_4bdf"








t <- seq(n)                                              # time axis
lags <- 0:m                                              # lag time series

verbose <- FALSE





true_kernel1 <- dgamma(t-1, shape=alpha1, scale=tau1/alpha1)   # true convolution kernel is gamma distribution
true_kernel2 <- dgamma(t-1, shape=alpha2, scale=tau2/alpha2)   # true convolution kernel is gamma distribution
true_kernel3 <- dgamma(t-1, shape=alpha3, scale=tau3/alpha3)   # true convolution kernel is gamma distribution


benchmk1 <- true_kernel1[1:(m+1)]                          # this is the benchmark that we should be able to reproduce
benchmk2 <- true_kernel2[1:(m+1)]                          # this is the benchmark that we should be able to reproduce
benchmk3 <- true_kernel3[1:(m+1)]                          # this is the benchmark that we should be able to reproduce



# first we need to calculate the correlated inputs
if (x_rho==0) {
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  randroot <- rnorm(n)
} else {
  x1 <- as.vector(unname(arima.sim(n=n, model=list(ar=x_rho))))           
  x2 <- as.vector(unname(arima.sim(n=n, model=list(ar=x_rho))))           
  x3 <- as.vector(unname(arima.sim(n=n, model=list(ar=x_rho))))           
  randroot <- as.vector(unname(arima.sim(n=n, model=list(ar=x_rho))))  
}

# then cross-correlate the inputs as required
x1 <- sqrt(cross_rho) * randroot + sqrt(1-abs(cross_rho))*x1       # using sqrt(rho) in place of the usual rho means that the correlation between TWO time series generated this way will be rho
x2 <- sqrt(cross_rho) * randroot + sqrt(1-abs(cross_rho))*x2       # using sqrt(rho) in place of the usual rho means that the correlation between TWO time series generated this way will be rho
x3 <- sqrt(cross_rho) * randroot + sqrt(1-abs(cross_rho))*x3       # using sqrt(rho) in place of the usual rho means that the correlation between TWO time series generated this way will be rho


# now convolve the inputs with their kernels and add them up to get the true y
y1 <- convolve(x1, true_kernel1, conj=FALSE)
y2 <- convolve(x2, true_kernel2, conj=FALSE)
y3 <- convolve(x3, true_kernel3, conj=FALSE)       # true output

y_true <- y1 + y2 + y3

# now add ARMA noise to y_true
if (is.null(y_err_model)) {
  y_err <- rnorm(n) 
} else y_err <- as.vector(unname(arima.sim(n=n, model=y_err_model)))      # error is ARMA noise

y_err <- y_err *(y_err_scale * sd(y_true)/sd(y_err))

y <- y_true + y_err                                      # measured output (true plus error)

zz <- IRF(y=y, x=cbind(x1,x2,x3), m=m, verbose=verbose)

bb1 <- as.numeric(zz$IRF[,1])
bb2 <- as.numeric(zz$IRF[,2])
bb3 <- as.numeric(zz$IRF[,3])

se1 <- as.numeric(zz$se[,1])
se2 <- as.numeric(zz$se[,2])
se3 <- as.numeric(zz$se[,3])




# now plot results


maxval <- max(benchmk1, bb1, na.rm=TRUE) #generate upper limit for plotting
maxval <- 1.05*maxval  #and extend this a bit
minval <- -0.1*maxval  #use standardized scaling so plots can be compared

# here is the benchmark, in blue
plot(lags, benchmk1, main="compare to benchmark1", xlab="lag", ylab="kernel", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20, cex.main=0.6)
# here is our result, in green
par(new=TRUE)
plot(lags, bb1, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20)
# and here are the uncertainty bounds (+/- one standard error), in gray
par(new=TRUE)
plot(lags, bb1+2*se1, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20)
par(new=TRUE)
plot(lags, bb1-2*se1, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20)

# here is the benchmark, in blue
plot(lags, benchmk2, main="compare to benchmark2", xlab="lag", ylab="kernel", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20, cex.main=0.6)
# here is our result, in green
par(new=TRUE)
plot(lags, bb2, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20)
# and here are the uncertainty bounds (+/- one standard error), in gray
par(new=TRUE)
plot(lags, bb2+2*se2, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20)
par(new=TRUE)
plot(lags, bb2-2*se2, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20)


# here is the benchmark, in blue
plot(lags, benchmk3, main="compare to benchmark3", xlab="lag", ylab="kernel", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20, cex.main=0.6)
# here is our result, in green
par(new=TRUE)
plot(lags, bb3, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20)
# and here are the uncertainty bounds (+/- one standard error), in gray
par(new=TRUE)
plot(lags, bb3+2*se3, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20)
par(new=TRUE)
plot(lags, bb3-2*se3, xlab="", ylab="", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=2, col="gray", pch=20)


# now calculate fraction of time that we are more than 2 standard errors away from the benchmark

# cat("fraction of misses in first kernel=",sum( abs(bb1-benchmk1)>(2*se1) )/(m+1), "\n")
# cat("fraction of misses in second kernel=",sum( abs(bb2-benchmk2)>(2*se2) )/(m+1), "\n")
# cat("fraction of misses in third kernel=",sum( abs(bb3-benchmk3)>(2*se3) )/(m+1), "\n")




maxval <- max(y, y_true, na.rm=TRUE)
minval <- min(y, y_true, na.rm=TRUE)

ll <- length(y)
ff <- max(ll-1000,0)
# here is y, in orange
plot(seq(from=ff,to=ll), y[ff:ll], main="true y and measured y", xlab="time step", ylab="y", ylim=c(minval, maxval), type="l", lwd=3, col="orange", pch=20, cex.main=0.6)
# here is y_true, in red
par(new=TRUE)
plot(seq(from=ff,to=ll), y_true[ff:ll], xlab="", ylab="",  ylim=c(minval, maxval), type="l", lwd=3, col="red", pch=20)
# here is y1, in green
par(new=TRUE)
plot(seq(from=ff,to=ll), y1[ff:ll], xlab="", ylab="",  ylim=c(minval, maxval), type="l", lwd=3, col="green", pch=20)
# here is y2, in blue
par(new=TRUE)
plot(seq(from=ff,to=ll), y2[ff:ll], xlab="", ylab="",  ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20)
# here is y1, in sea green
par(new=TRUE)
plot(seq(from=ff,to=ll), y3[ff:ll], xlab="", ylab="",  ylim=c(minval, maxval), type="l", lwd=3, col="seagreen", pch=20)


fwrite(data.table(cbind(t=1:n, x1, x2, x3, y1, y2, y3, y_true, y_err, y)), file=paste0(fileID, "_time_series.txt"), row.names=FALSE, sep="\t")

fwrite(data.table(cbind(lags, benchmk1, bb1, se1, benchmk2, bb2, se2, benchmk3, bb3, se3)), file=paste0(fileID, "_IRF.txt"), row.names=FALSE, sep="\t")



