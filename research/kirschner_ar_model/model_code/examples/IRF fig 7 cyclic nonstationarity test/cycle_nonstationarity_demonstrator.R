setwd("~/IRF paper")    # change this to wherever the called scripts are


rm(list=ls())  #clear environment

library(data.table)
library(dplyr)


source("IRFnnhs.R")  # impulse response function script



# first benchmark test
set.seed(456)
n <- 1e4              # length of time series
m <- 60               # longest lag to be analyze
nx <- 4               # number of different input classes to be demixed

alpha <- rep(0, 24)
tau <- rep(0, 24)
alpha_max <- 4
alpha_min <- 0.5
tau_max <- 40
tau_min <- 10



x_rho <- 0.0          # serial correlation coefficient for the input time series

y_err_scale <- 1.0   # RMS amplitude of y errors, as a fraction of the standard deviation of the true y values
y_err_model <- list(ar=c(0.9), ma=c(-0.2, 0.2))
fileID <- "nonstationary_cycle_Fig_7"




t <- seq(n)                                              # time axis
lags <- 0:m                                              # lag time series

verbose <- TRUE


# define convolution kernels
hourly_kernel <- matrix(data=0, nrow=n, ncol=24)
for (i in 1:24) {
  alpha <- alpha_max-(alpha_max-alpha_min)*(0.5+0.5*cos(2*pi*i/24))
  tau <- tau_min+(tau_max-tau_min)*(0.5+0.5*sin(2*pi*i/24))
  
  hourly_kernel[,i] <- pgamma(t, shape=alpha, scale=tau/alpha)-pgamma(t-1, shape=alpha, scale=tau/alpha)
  # defining the kernel this way allows alpha values <1 (since it's the *average* density of the gamma distribution over the interval)
}


# save benchmarks (convolution kernels, truncated to the number of lags that we analyze)
benchmk <- matrix(data=0, nrow=m+1, ncol=nx)
for (i in 1:nx) 
  for (j in 1:6) benchmk[,i] <- benchmk[,i] + hourly_kernel[1:(m+1),1+(6*i+j-5)%%24]/6
colnames(benchmk) <- paste0("benchmk_", 1:nx)


# first we need to calculate the correlated inputs
if (x_rho==0) {
  x <- rnorm(n)
} else {
  x <- as.vector(unname(arima.sim(n=n, model=list(ar=x_rho))))           
}

x <- x^2   # square x so inputs are always positive (and follow chi-square distribution with 1 degree of freedom)



# now generate the output time series by convolving each hour's x with the corresponding hourly kernel
y_true <- rep(0, n)
ys <- matrix(data=0, nrow=n, ncol=nx)
for (i in 1:24) y_true <- y_true + convolve(ifelse((((t-1)%%24)+1==i), x, 0), hourly_kernel[,i], conj=FALSE)
# the ifelse creates a time series with values of x for at hour i each day, otherwise zero

ys <- matrix(data=0, nrow=n, ncol=nx)
for (i in 1:nx) 
  for (j in 1:6) ys[,i] <- ys[,i] + convolve(ifelse((((t-1)%%24)+1==1+(6*i+j-5)%%24), x, 0), hourly_kernel[,i], conj=FALSE)

ysum <- matrix(data=0, nrow=n, ncol=nx)
ysum[,1] <- ys[,1]
for (i in 2:nx) ysum[,i] <- ysum[,(i-1)] + ys[,i]




# now split the input time series into a matrix of the x's for each of the types/classes
xs <- matrix(data=0, nrow=n, ncol=nx)
for (i in 1:nx) 
  for (j in 1:6) xs[,i] <- xs[,i] + ifelse((((t-1)%%24)+1==1+(6*i+j-5)%%24), x, 0)

colnames(xs) <- paste0("xs_", 1:nx)


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
ff <- max(ll-500,0)
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





