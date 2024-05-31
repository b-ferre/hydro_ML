setwd("~/IRF paper")    # change this to wherever the called scripts are


rm(list=ls())  #clear environment

library(data.table)
library(dplyr)


source("IRFnnhs.R")  # impulse response function script



# benchmark test
set.seed(456)
n <- 1e4              # length of time series
m <- 40               # longest lag to be analyze
nx <- 3               # number of different input classes to be demixed

alpha <- rep(0, nx)
tau <- rep(0, nx)
alpha[1] <- 1
alpha[2] <- 2
alpha[3] <- 4
tau[1] <- 5
tau[2] <- 10
tau[3] <- 20

mean_interval_length <- 5

x_rho <- 0.0          # serial correlation coefficient for the input time series

y_err_scale <- 1.0    # RMS amplitude of y errors, as a fraction of the standard deviation of the true y values
y_err_model <- list(ar=c(0.9), ma=c(-0.2, 0.2))
fileID <- "nonstationarity_benchmark_test_Fig_6"














t <- seq(n)                                              # time axis
lags <- 0:m                                              # lag time series

verbose <- TRUE


# define convolution kernels
true_kernel <- matrix(data=0, nrow=n, ncol=nx)
# for (i in 1:nx) true_kernel[,i] <- dgamma(t-1, shape=alpha[i], scale=tau[i]/alpha[i])
# defining the kernel this way allows alpha values <1 (since it's the *average* of gamma over the interval)
for (i in 1:nx) true_kernel[,i] <- pgamma(t, shape=alpha[i], scale=tau[i]/alpha[i])-pgamma(t-1, shape=alpha[i], scale=tau[i]/alpha[i])


# save benchmarks (convolution kernels, truncated to the number of lags that we analyze)
benchmk <- matrix(data=0, nrow=m+1, ncol=nx)
for (i in 1:nx) benchmk[,i] <- true_kernel[1:(m+1),i]
colnames(benchmk) <- paste0("benchmk_", 1:nx)


# set up an index "xnum", taking on values between 1 and nx, that indicates which kernel is used at each time step
randx <- ceiling(runif(n)*nx)
xnum <- ifelse(runif(n)<1.0/(mean_interval_length*(nx-1)/nx), randx, 0)

xnum[1] <- randx[1]
for (i in 2:n) if (xnum[i]==0) xnum[i] <- xnum[i-1]


# first we need to calculate the correlated inputs
if (x_rho==0) {
  x <- rnorm(n)
} else {
  x <- as.vector(unname(arima.sim(n=n, model=list(ar=x_rho))))           
}

x <- x^2   # square x so inputs are always positive


# now split the input time series into a matrix of the x's for each of the types/classes
xs <- matrix(data=0, nrow=n, ncol=nx)
for (i in 1:nx) xs[,i] <- ifelse((xnum==i), x, 0)
colnames(xs) <- paste0("xs_", 1:nx)

# now generate the output time series by convolving each of the columns of xs with the corresponding kernel
ys <- matrix(data=0, nrow=n, ncol=nx)    # output time series corresponding to inputs in each type/class xs
ysum <- matrix(data=0, nrow=n, ncol=nx)  # sum of these output time series, in inverse order (the last one, then the sum of the last two, then the sum of the last 3, etc.)
y_true <- rep(0, n)

for (i in 1:nx) ys[,i] <- convolve(xs[,i], true_kernel[,i], conj=FALSE)
for (i in 1:nx) y_true <- y_true + ys[,i]
ysum[,nx] <- ys[,nx]
for (i in 1:(nx-1)) ysum[,(nx-i)] <- ysum[,(nx-i+1)] + ys[,(nx-i)]

colnames(xs) <- paste0("xs_", 1:nx)
colnames(ys) <- paste0("ys_", 1:nx)
colnames(ysum) <- paste0("ysum_", 1:nx)




# now add ARMA noise to y_true
if (is.null(y_err_model)) {
  y_err <- rnorm(n) 
} else y_err <- as.vector(unname(arima.sim(n=n, model=y_err_model)))      # error is ARMA noise

y_err <- y_err *(y_err_scale * sd(y_true)/sd(y_err))

y <- y_true + y_err                                      # measured output (true plus error)



# now analyze this using IRF
zz <- IRF(y=y, xs, m=m, verbose=verbose)


bb <- zz$IRF
se <- zz$se
colnames(bb) <- paste0("bb_", 1:nx)
colnames(se) <- paste0("se_", 1:nx)



# now plot results


maxval <- max(benchmk, bb, na.rm=TRUE) #generate upper limit for plotting
maxval <- 1.05*maxval  #and extend this a bit
minval <- -0.1*maxval  #use standardized scaling so plots can be compared


# now plot against all benchmarks
for (i in 1:nx) {
  # here is the benchmark, in blue
  plot(lags, benchmk[,i], main=paste0("compare to benchmark ", i), xlab="lag", ylab="kernel", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20, cex.main=0.6)
  # here is our result, in green
  lines(lags, bb[,i], type="l", lwd=3, col="green", pch=20)
  # and here are the uncertainty bounds (+/- one standard error), in gray
  lines(lags, bb[,i]+2*se[,i], type="l", lwd=2, col="gray", pch=20)
  lines(lags, bb[,i]-2*se[,i], type="l", lwd=2, col="gray", pch=20)
}





# now plot time series
maxval <- max(y, y_true, na.rm=TRUE)
minval <- min(y, y_true, na.rm=TRUE)

ll <- length(y)
ff <- max(ll-1000,0)
# here is y, in orange
plot(seq(from=ff,to=ll), y[ff:ll], main="true y and measured y", xlab="time step", ylab="y", ylim=c(minval, maxval), type="l", lwd=3, col="orange", pch=20, cex.main=0.6)
# here is y_true, in red
lines(seq(from=ff,to=ll), y_true[ff:ll], type="l", lwd=3, col="red", pch=20)

# and plot each of the ys components
# here is y, in orange
plot(seq(from=ff,to=ll), y[ff:ll], main="true y and measured y", xlab="time step", ylab="y", ylim=c(minval, maxval), type="l", lwd=3, col="orange", pch=20, cex.main=0.6)
# here is y_true, in red
for (i in 1:nx) lines(seq(from=ff,to=ll), ysum[ff:ll, i], type="l", lwd=3, col="cyan", pch=20)




fwrite(data.table(cbind(t=1:n, xs, x, ys, y_true, y_err, y)), file=paste0(fileID, "_time_series.txt"), row.names=FALSE, sep="\t")




fwrite(data.table(cbind(lags, benchmk, bb, se)), file=paste0(fileID, "_IRF.txt"), row.names=FALSE, sep="\t")









zz <- IRF(y=y, x=x, m=m, verbose=verbose)

bb <- as.numeric(zz$IRF[,1])

se <- as.numeric(zz$se[,1])

bench <- rowMeans(benchmk)

# now plot results


maxval <- max(benchmk, bb, na.rm=TRUE) #generate upper limit for plotting
maxval <- 1.05*maxval  #and extend this a bit
minval <- -0.1*maxval  #use standardized scaling so plots can be compared

# here is the benchmark, in blue
plot(lags, bench, main="compare aggregated fit to benchmark", xlab="lag", ylab="kernel", xlim=c(0, max(lags)+1), ylim=c(minval, maxval), type="l", lwd=3, col="blue", pch=20, cex.main=0.6)
for (i in 1:nx) lines(lags, benchmk[,i], col="cyan")
# here is our result, in green
lines(lags, bb,  type="l", lwd=3, col="green", pch=20)
# and here are the uncertainty bounds (+/- two standard errors), in gray
lines(lags, bb+2*se, type="l", lwd=2, col="gray", pch=20)
lines(lags, bb-2*se, type="l", lwd=2, col="gray", pch=20)







fwrite(data.table(cbind(lags, benchmk, bench, bb, se)), file=paste0(fileID, "_averageIRF from aggregated input.txt"), row.names=FALSE, sep="\t")





