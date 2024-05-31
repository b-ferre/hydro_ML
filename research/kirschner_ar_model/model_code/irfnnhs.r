# IRFnnhs.R: Impulse Response Functions for nonlinear, nonstationary, and heterogeneous systems,
#                  estimated by deconvolution and demixing of time series,
#                  with corrections for contamination by ARMA and ARIMA noise
# 
#
#
# version 1.2  build 2022.04.24
# Author: James Kirchner, ETH Zurich
#
# Copyright (C) 2022 ETH Zurich and James Kirchner
# Public use of this script is permitted under GNU General Public License 3 (GPL3); 
# for details see <https://www.gnu.org/licenses/>
# 
#
# READ THIS CAREFULLY:
# ETH Zurich and James Kirchner make ABSOLUTELY NO WARRANTIES OF ANY KIND, including NO WARRANTIES, expressed or implied, that this software is
#     free of errors or is suitable for any particular purpose.  Users are solely responsible for determining the suitability and
#     reliability of this software for their own purposes.
#
# ALSO READ THIS:
# These scripts implement techniques presented in J.W. Kirchner, "Impulse response functions for nonlinear, nonstationary, and heterogeneous systems, 
#     estimated by deconvolution and de-mixing of noisy time series", Sensors, 2022 (hereafter, K2022)  Users should cite that paper.
#
# Equations denoted K2019 refer to J.W. Kirchner, Quantifying new water fractions and transit time distributions using ensemble 
#     hydrograph separation: theory and benchmark tests, Hydrology and Earth System Sciences, 23, 303-349, https://doi.org/10.5194/hess-23-303-2019, 2019.
#
# Equations denoted KK2020 refer to J.W. Kirchner and J.L.A. Knapp, Technical note: Calculation scripts for ensemble hydrograph 
#     separation, Hydrology and Earth System Sciences, 24, 5539-5558, https://doi.org/10.5194/hess-24-5539-2020, 2020.
#
# The equations may differ in detail.  For example, in these scripts, lag indices range from 1 to m+1 instead of 0 to m 
# due to the array indexing conventions of the R language. 



# SOME NOTES ABOUT PERFORMANCE:
#
# These scripts require solving large matrix problems when multiple input signals are deconvolved over many lags using long
# time series.  The order of difficulty for these matrix problems is approximately n * (m+1)^2 * nx^2, where n is the number 
# of time steps, nx is the number of input signals, and m+1 is the number of time lags over which their impulse response functions are estimated.  
#
# R version 4 (which is current as of early 2022) does not use particularly efficient BLAS (Basic Linear Algebra Subprograms) routines by default.
# Switching to a processor-optimized BLAS (such as the Intel MKL, in the case of Windows computers) can speed up these matrix operations by an 
# order of magnitude or more.  A straightforward way to make the Intel MKL available is to install Microsoft R Open (https://mran.microsoft.com/open).
# A web search will reveal other approaches for other operating systems and processors.
#
# Further (modest) performance gains *may* be possible by using the R command "options(matprod = "blas")".  This will invoke the BLAS without checking for NaN or Inf values, 
# which (if they are present), would lead to unpredictable results.  Thus this option should be used with caution.
#
# These large matrix problems could also lead to memory constraints.  The size of the design matrix X scales as n * (m+1) * nx (where, again, n is 
# the number of time steps, nx is the number of input time series, and m+1 is the number of time lags over which their IRFs are estimated).  The solvAR routine saves memory 
# by creating the X matrix in discrete chunks of rows (=time steps), calculating the necessary cross-products for each chunk, and then combining the cross-products.  
# This approach exploits the fact that the m+1 columns for each of the nx input signals are lags of one another, so they do not all need to be stored
# in memory.  This incurs a negligible performance penalty, and reduces the largest matrices that must be stored to order n * nx.
#
# For robust estimation it can be advantageous to create the design matrix all in one chunk, so that it can be re-used in each iteration rather than needing
# to be re-created.  However, if the design matrix is large enough to trigger memory paging, this advantage is lost.  The (approximate) upper limit on the size of the
# design matrix, above which it will be chunked, is determined by the parameter max.chunk.  Users who are experiencing unexpectedly slow performance can experiment
# with changing max.chunk to see whether that helps.
#
# For efficiency reasons, three time-consuming subroutines are coded in C++ under Rcpp.  






# if any of these packages are not already installed, do install.packages("<package name>")
# Windows users will also need to install rtools (see https://cran.r-project.org/bin/windows/Rtools/rtools40.html)
# for the necessary GNU C++ compiler; this should already be available for Mac and Linux implementations
library(dplyr)  
library(matrixStats)  
library(Rcpp)




TICK <- Sys.time()



#/////////////////////////////////////
#### create regularization matrix ####
#/////////////////////////////////////

tikh <- function(m) {     # creates an m x m Tikhonov - Phillips regularization matrix (equation 49 of K2019)
  
  v <- rep(0,m)           # start with a vector of zeroes
  v[1:3] <- c(6, -4, 1)    
  x <- toeplitz(v)        # make this into a toeplitz matrix
  
  #first two rows are special
  x[1, 1:2] <- c(1, -2)
  x[2, 1:3] <- c(-2, 5, -4)
  
  #last two rows are special
  x[m, (m-1):m] <- c(-2, 1)
  x[(m-1), (m-2):m] <- c(-4, 5, -2)
  
  return(x)
} # end tikh

#////////////////////////////////////
# END OF create regularization matrix
#////////////////////////////////////








#////////////////////////////////////////////////////////////////
#### Rcpp function to efficiently calculate predicted values ####
#////////////////////////////////////////////////////////////////
cppFunction('NumericVector calc_pred(      //THIS FUNCTION IS FOR USE IN SOLVAR, NOT IRF
               NumericVector x, 
               NumericVector y,
               NumericVector beta,
               int nx,
               int n,
               int nlag,
               int h) {

  int i, lag;
  
// x is an matrix of nx columns and n rows (which is handled here as if it were a vector, for efficiency and clarity)
// nlag is the number of lags in the xx matrix; h is the number of lags that we add from the y vector

  NumericVector pred(n-nlag+1, beta[beta.length()-1]);            //declare and initialize the target vector to constant term

  for (i=0; i<nx; i++) {             // step through columns of x
      for (lag=0; lag<nlag; lag++) {      // step through lags
          pred = pred + beta(i*nlag+lag)*x[Rcpp::Range((i*n+nlag-1-lag),(i*n+nlag-1-lag+pred.length()))];
      }
  }
  
  // now do lags of y for AR correction
  if (h>0) for (i=0; i<h; i++) {     // step through h
          pred = pred + beta(nx*nlag+i)*y[Rcpp::Range((nlag-1-(i+1)),(nlag-1-(i+1)+pred.length()))];
   }
  
  return pred;
  }
')





#////////////////////////////////////////////////////////////////
#### Rcpp function to efficiently calculate predicted values ####
#////////////////////////////////////////////////////////////////
cppFunction('NumericVector calc_ypred(      //THIS FUNCTION IS FOR USE IN IRF, NOT SOLVAR
               NumericVector x, 
               NumericVector beta,
               int nx,
               int n,
               int m,
               int h) {

  int i, lag, nlag;
  
  nlag = m+h+1;
  
// x is an matrix of nx columns and n rows (which is handled here as if it were a vector, for efficiency and clarity)
// nlag is the number of lags in the xx matrix; h is the number of lags that we add from the y vector

  NumericVector pred(n-nlag+1);            //declare and initialize the target vector

  for (i=0; i<nx; i++) {             // step through columns of x
      for (lag=0; lag<=m; lag++) {      // step through lags
          pred = pred + beta(i*(m+1)+lag)*x[Rcpp::Range((i*n+nlag-1-lag),(i*n+nlag-1-lag+pred.length()))];
      }
  }
  
  return pred;
  }
')










#//////////////////////////////////////////////////////
#### Rcpp function to efficiently create xx matrix ####
#//////////////////////////////////////////////////////
cppFunction('NumericVector buildxx(
               NumericVector x, 
               NumericVector y,
               int nx,
               int n,
               int nlag,
               int h,
               int chunktop,
               int chunkht) {


  NumericVector::iterator xroot = x.begin();
  NumericVector::iterator yroot = y.begin();
  NumericVector::iterator snip1;
  NumericVector::iterator snip2;
  NumericVector::iterator dest;

  int i, lag;
  
// x is an matrix of nx columns and n rows (which is handled here as if it were a vector, for efficiency and clarity)
// nlag is the number of lags in the x matrix; h is the number of lags that we add from the y vector

  NumericVector xx((nx*nlag+h+1)*chunkht);            //declare and initialize the target matrix (here shown only as a vector -- this will appear as if it were a matrix in the result)
  NumericVector::iterator xxroot = xx.begin();

  for (i=0; i<nx; i++) {             // step through columns of x
      for (lag=0; lag<nlag; lag++) {      // step through lags
          //this is the vector location where we need to start copying: move over to column i, down to chunktop (minus 1 because indexing starts at zero in C instead of 1 in R), 
          //and then back up by lag:
          snip1 = xroot + (i*n) + (chunktop-1) - lag; 
          snip2 = snip1 + chunkht;                           //this is the vector location where we should stop copying
          dest = xxroot + (i*chunkht*nlag) + lag*chunkht;    //this is where we should paste the vector in (each column of x becomes nlag columns here, each of length chunkht)
          std::copy(snip1, snip2, dest);
      }
  }
  
  // now do lags of y for AR correction (again we do this with iterators)
  if (h>0) for (i=0; i<h; i++) {     // step through h
          snip1 = yroot + (chunktop-1) - (i+1);              //step down to chunktop (with -1 due to indexing difference between the C and R languages), then back up by i+1 because i=0 corresponds to lag=1 here...
          snip2 = snip1 + chunkht;
          dest = xxroot + (nx*nlag*chunkht) + (i*chunkht);
          std::copy(snip1, snip2, dest);
  }
  
  // now add the constant column
  NumericVector unity(chunkht);
  std::fill(unity.begin(), unity.end(), 1.0);
  std::copy(unity.begin(), unity.end(), xxroot+(nx*nlag+h)*chunkht);
  return xx;
  }
')







#///////////////////////////////////////////////////////////////////////
#### calculate weighted autocorrelation and partial autocorrelation ####
#///////////////////////////////////////////////////////////////////////

wtd.pacf <-
  
  function(y,
           wt=rep(1, NROW(y)) , 
           maxlag=NULL )                    # maximum lag at which acf and pacf will be estimated
  {
    
    if (NROW(y)!=NROW(wt)) stop("error in wtd.pacf: data vector and weight vector must be same length")
    if (min(wt, na.rm=TRUE)<0) stop("error in wtd.pacf: weights cannot be less than zero")
    
    wt[is.na(wt)] <- 0                      # convert NA's to zero weight
    
    wt[is.na(y)] <- 0                       # set weight for any missing values to zero
    
    y[is.na(y)] <- 0                        # zero out NA values (but don't delete them: we need to preserve the temporal ordering!)
    
    n <- NROW(y)                            # length of input vector
    
    n.eff <- (sum(wt)^2)/sum(wt^2)          # effective sample size
    
    if (n.eff<30) stop("error in wtd.pacf: effective sample size of at least 30 is needed to estimate acf reliably")
    
    if (is.null(maxlag)) maxlag <- floor(min(n.eff/2, 10*log10(n.eff)))
    
    y <- y - sum(wt*y)/sum(wt)              # remove weighted mean from y
    
    acf <- rep(0, maxlag)                   # initialize acf vector
    pacf <- rep(0, maxlag)                  # initialize pacf vector
    
    #    y.sw <- y * sqrt(wt)                    # y weighted by square root of wt
    
    #    for (i in 1:maxlag) acf[i] <- cor(y.sw, dplyr::lag(y.sw, i), use="complete.obs")   # autocorrelation
    
    # this seems to behave more sensibly than the approach above with robustness weights (which makes sense because the full robustness weight is applied at each time step, not the square root) 
    for (i in 1:maxlag) acf[i] <- cor(y*wt, dplyr::lag(y, i), use="complete.obs")   # autocorrelation 
    
    vacf <- c(1, acf[-length(acf)])         # copy of autocorrelation, now for lags of zero to maxlag-1  (for toeplitz matrix that we use to calculate pacf)
    
    pacf[1] <- acf[1]/vacf[1]               # first-order partial autocorrelation equals autocorrelation at lag 1, divided by autocorrelation at lag zero (which =1)
    
    for (i in 2:maxlag) pacf[i] <- solve(toeplitz(vacf[1:i]), acf[1:i])[i]  # matrix solution for pacf
    # see, for example, https://stats.stackexchange.com/questions/129052/acf-and-pacf-formula
    
    return(list(
      acf=acf, 
      pacf=pacf))
    
  }










#////////////////////////////////////////////////////////////////////////////
#### calculate cross-products, (with chunking of xx matrix if necessary) ####
#////////////////////////////////////////////////////////////////////////////

calc_crossprods <-
  function(y,                               # any missing values must have been converted to zero
           x,                               # any missing values must have been converted to zero!
           wt=rep(1, NROW(x)) ,             # this must already have zeroes for any rows where x or its lags will be NA
           maxlag ,                         # this is n_lags in calling routine; does not include 1 more for lag zero
           h ,
           xx = NULL ,
           chunk.maxht ,
           verbose = FALSE )
  {
    
    nx <- NCOL(x)
    n <- NROW(x)
    nlags <- maxlag+1                                          # add one to account for lag zero
    wt[1:maxlag] <- 0                                          # just in case this hasn't been done already... the first maxlag rows will not be valid b/c one or more lags are missing
    sqwt <- as.vector(sqrt(wt))                                # sw is square root of weight (including robustness weights)
    
    if (!is.null(xx)) {                                        
      
      # if we already have a one-piece xx matrix, we can use it directly without compiling it again
      C <- crossprod(sqwt[(maxlag+1):n]*xx)                                        # use R's very fast cross-product routine (fast with a good BLAS, that is)
      xy <- crossprod(sqwt[(maxlag+1):n]*xx, sqwt[(maxlag+1):n]*y[(maxlag+1):n])   # note that because x and y are multiplied by sqrt(wt), these cross-products are weighted
      n.nz <- colSums(((xx*sqwt[(maxlag+1):n]) != 0))                              # count non-zero elements (with non-zero weight) in each x column
      
    } else {                                                   
      
      # if we don't already have an xx matrix, we'll have to build it, potentially with chunking
      # initialize outputs at zero because we'll be adding to them with each chunk
      C <- matrix(0, nrow=nx*nlags + h + 1, ncol=nx*nlags + h + 1)     # cross product matrix has an extra row and column for the constant term
      xy <- rep(0, nx*nlags + h + 1)                                   # x-y cross products (with an extra row for the constant term)
      n.nz <- rep(0, nx*nlags + h + 1)                                 # this vector tallies non-zero elements of x (with nonzero weight)
      
      
      chunkend <- maxlag                                           # initializing this at the last row that's guaranteed to be infeasible due to NA's in lags;
      # it'll change as soon as we enter the loop to build the xx matrix
      
      repeat{
        chunktop <- chunkend + 1                                   # new chunk starts on next step after the last one ended
        chunkend <- min(chunktop+chunk.maxht-1, n)                 # new chunk ends at chunk.top+chunk.maxht (or at the end of the data set)
        chunkht <- chunkend-chunktop + 1                           # height of this chunk
        
        if (sum(sqwt[chunktop:chunkend])>0) {                      # if all of the rows in this chunk have weight of zero then we can skip this one
          xx <- buildxx(x, y, nx, n, nlags, h, chunktop, chunkht)  # here is the Rcpp routine to efficiently build the xx matrix
          xx <- matrix(xx, nrow=chunkht)
          yy <- y[chunktop:chunkend]                               # take corresponding chunk of y
          sw <- sqwt[chunktop:chunkend]                            # and corresponding chunk of sqrt(weight)
          
          C <- C + crossprod(sw*xx)                                # use R's very fast cross-product routine (fast with a good BLAS, that is)
          xy <- xy + crossprod(sw*xx, sw*yy)                       # note that because x and yy are multiplied by sqrt(wt), these cross-products are weighted
          
          n.nz <- n.nz + colSums(((xx*sw) != 0))                   # count non-zero elements (with non-zero weight) in each x column and add to any tallies from previous chunks
          
        } # end if
        
        if (verbose) if (chunkend<n) cat(paste0(ceiling(100*chunkend/n), "% of cross-products done...: ", round(difftime(Sys.time(), TICK, units="secs"), 3), "                           \r"))
        
        if (chunkend >= n) {break}                                 # exit if this is the last chunk
      } # end repeat
      
    } # end else (i.e., we do't already have an xx matrix)
    
    
    if (chunk.maxht>=n) {
      return(list(                                                 # if we have a complete xx matrix, pass it along so it can be re-used...
        C = C ,
        xy = xy ,
        n.nz = n.nz ,
        xx = xx))
    } else {
      return(list(
        C = C ,
        xy = xy ,
        n.nz = n.nz ,
        xx = NULL))
    }
  }



















#////////////////////////////////////////////////////////////////////////////////////////////
#### Robust solution to matrix equations for IRF, with chunking of x matrix if necessary ####
#////////////////////////////////////////////////////////////////////////////////////////////

RobustSolvAR <-
  function(y, 
           x, 
           wt ,
           FD = FALSE ,
           h = 0 ,
           m ,
           nu = 0 ,
           robust = FALSE ,
           verbose = FALSE,
           chunk.maxht = chunk.maxht)
    
  {
    
    # takes as input:
    # y              vector of y values
    # x              matrix of x values
    # wt             vector of weight values
    # FD             flag for implementation of first differencing.  If FD==FALSE, no first differencing is applied.  
    #                          If FD==TRUE, y is first-differenced but x is not, and the coefficients and their covariance 
    #                          matrices are adjusted accordingly (as explained in K2022).
    # h              order of optional AR correction
    # m              integer value of maximum lag (first lag is zero, so total number of lags is m+1).  This doesn't count additional lags needed for AR correction or first differences  
    # nu             fractional weight, 0-1, to be given to Tikhonov-Phillips regularization (default is 0 = no regularization)
    # robust         flag controlling whether robust estimation by Iteratively Reweighted Least Squares (IRLS) will be used
    # verbose        controls whether progress reports are printed (TRUE) or suppressed (FALSE)
    # chunk.maxht    maximum height of design matrix that will not trigger chunking
    
    # returns a list as output, with the following objects:
    # b.hat          regression coefficients (matrix of nx columns and m+1 rows)
    # se             standard errors (matrix of nx columns and m+1 rows)
    # Kbb            stack of covariance matrices of coefficients for each x (array of m+1 by m+1 matrices)
    # phi            fitted AR coefficients (vector of length=h)
    # n              total number of time steps
    # n.eff          effective number of time steps, accounting for uneven weighting
    # n.nz           number of nonzero elements in the *weighted* x's (so values with zero weight don't count), stored as matrix of nx columns and m+1 rows
    # e              residuals
    # s2e            weighted residual variance
    # resid.acf      autocorrelation function of residuals
    # resid.pacf     partial autocorrelation function of residuals
    # rw             vector of robustness weights
    
    
    # If robust==TRUE, runtimes may be substantially faster if the robust solvAR can avoid compiling the cross-products in chunks, and instead create the
    # entire xx matrix all at once, so it can be re-used for each iteration of the robust estimation loop.  
    # This requires enough memory to hold this matrix (and potentially several copies of it) without triggering memory paging. 
    # On Windows, use memory.limit(size=###) to change the memory available to R, where ### is in MB.  
    
    
    maxlag <- m+h                   # total number of lags (not counting one for lag zero)
    if (FD) maxlag <- maxlag+1      # one extra lag if we are using first differences
    
    nx <- NCOL(x)                   # this is the number of different x variables
    
    mm <- (maxlag+1)*nx             # total number of b coefficients (and width of the x matrix, not counting the AR terms or the constant term)
    
    
    
    iter_limit <- 200               # iteration limit for robust estimation
    
    
    
    #///////////////////////////////////
    # zero out missing values of y and x
    
    orig_y <- y                                               # first make copies of y and x (for calculating residuals later)
    orig_x <- x
    orig_wt <- wt
    
    # note that all of these cases, and their lags, correspond to rows where wt has already been zeroed (in IRF), 
    # so replacing them with zeroes will not affect the values of the cross products (but will keep crossprod from choking on NA's)
    y[is.na(y)] <- 0                                          # replace missing values of y with zeroes
    x[is.na(x)] <- 0                                          # replace missing values of x with zeroes
    
    n <- length(y)                                            # number of rows
    
    if (nu>0) {                                               # if we regularize...:
      # We can't use one big regularization matrix for the whole problem, because it will lead to leakage between the end of one set of b coefficients and the start of another. 
      # Instead we need to create a segmented regularization matrix, with nx Tikhonov-Phillips matrices along the diagonal
      tt <- tikh(1+maxlag)                                    # start with one Tikhonov-Phillips matrix (this is our building block)
      H <- matrix(0, nrow = mm+h+1, ncol=mm+h+1)              # start with a matrix of all zeroes
      for (j in 1:nx)                                         # copy our building block into the right range(s) along the diagonal of H
        H[((j-1)*(1+maxlag)+1):(j*(1+maxlag)) , ((j-1)*(1+maxlag)+1):(j*(1+maxlag)) ] <- tt
    }    
    
    
    rw <- rep(1, n)                                           # initial values of robustness weights are all 1  NOTE: in contrast to earlier versions, here rw is of length n rather than n-maxlag.  That makes everything cleaner
    old_rw <- rw
    
    #////////////////////////////
    # robustness loop starts here
    iter <- 0                                                 # iteration counter
    xx <- NULL                                                # this will force evaluation of xx
    
    repeat{
      
      iter <- iter+1                                          # update counter
      
      
      wt <- rw*orig_wt                                        # multiply weights by robustness weight
      shortwt <- wt[(maxlag+1):n]                             # to avoid subsetting wt multiple times below
      
      sumwt <- sum(wt, na.rm=TRUE)
      
      
      # here is where we construct the xx matrix (with chunking if necessary) and calculate the crossproduct matrix C and the crossproduct vector xy
      # list2env will create local values of C, xy, n.nz, and xx.  xx will be NULL if chunking was required.
      list2env(    calc_crossprods(y=y, x=x, wt=wt, maxlag=maxlag, h=h, xx=xx, chunk.maxht=chunk.maxht, verbose=verbose),    envir=environment())
      
      C <- C/sumwt                                           # dividing by the sums of the weights does not change the regression parameters, 
      xy <- xy/sumwt                                         # but is necessary for consistency with the way uncertainties are calculated next
      
      if (nu==0) b.hat <- solve(C, xy)                       # if we *don't* regularize, solve for the un-regularized regression coefficients
      else {                                                 # if we *do* regularize...:
        # convert weighting parameter nu (values 0 - 1) to lambda as used in conventional Tikhonov-Phillips regularization formulas
        lambda <- nu / (1.0 - nu) * sum(diag(C)[1:(nx*(1+maxlag))]) / sum(diag(H))              # equation 50 of K2019 (omitting the diagonal elements that correspond to the constant term and AR terms)
        b.hat <- solve((C + (lambda * H)), xy) # now solve (regularized) linear system          # equation 46 of K2019, and equation A12 of KK2020
      }  
      
      # calculate predicted values -- note that rows with NA's will not count anyhow because they will have weight of zero in calculating error variance
      # note that our predictions don't include the first maxlag rows (because our xx matrix doesn't either...)
      
      if (!is.null(xx)) {                                     # if we have the xx matrix, we can calculate the predicted y values by matrix multiplication methods (fast with good BLAS)
        pred <- xx %*% b.hat
      } else {
        pred <- calc_pred(y=y, x=x, beta=b.hat, nx=nx, n=n, nlag=maxlag+1, h=h)
        pred[shortwt==0] <- NA
      }
      
      
      e <- orig_y[(maxlag+1):n] - pred                        # these are the residuals  (equation A13 of KK2020)
      e[shortwt==0] <- NA                                     # exclude all residuals for excluded rows (rows with zero weight)
      
      sumwt <- sum(shortwt)
      n.eff <- sumwt^2/sum(shortwt*shortwt)                   # calculate the effective sample size (accounting for uneven weighting)
      
      e2 <- e*e                                               # squared errors
      
      s2e <- (n.eff-1)/(n.eff-(mm+h+1)-1) * sum( e2 * shortwt, na.rm=TRUE)/sumwt   # weighted residual variance, corrected for degrees of freedom
      
      if (robust==FALSE) {break}                              # if we're not doing robust estimation, we can jump out of the loop here
      
      # now do iteratively reweighted least squares
      med.e <- sqrt(median(e2[shortwt!=0], na.rm=TRUE))                                                     # median absolute residual
      if (med.e==0) w <- ifelse(x==0, 1, 0)                                                                 # an unlikely exotic case: if MAR is zero, only 0 or 1 weights are possiblese
      else rw <- ifelse(is.na(e), 1, 1/(1+(e/(3.536*med.e))^2))                                             # Cauchy weight function
      
      rw <- c(rep(0, times=maxlag), rw)                       # add leading zeroes to bring rw to length n
      
      
      # exit criterion: 99 percent of the robustness weights have changed by less than 0.1 (absolute value) since previous iteration
      if (quantile(abs(rw-old_rw)[shortwt!=0], 0.99)<0.1) {  # don't count rw's and residuals for rows that have no weight!!
        if (verbose) cat("robustness iteration ", iter, " complete at time ", round(difftime(Sys.time(), TICK, units="secs"), 3),
                         "  s2e=", round(s2e, 8), " max e=", round(sqrt(max(e2, na.rm=TRUE)), 3),
                         "rw chg= ", round(quantile(abs(rw-old_rw)[shortwt!=0], 0.99), 3), "\n")
        {break}                            
      }
      
      if (verbose) cat("robustness iteration ", iter, " complete at time ", round(difftime(Sys.time(), TICK, units="secs"), 3),
                       "  s2e=", round(s2e, 8), " max e=", round(sqrt(max(e2, na.rm=TRUE)), 3),
                       "rw chg= ", round(quantile(abs(rw-old_rw)[shortwt!=0], 0.99), 3), "\n")
      
      old_rw <- rw            # save a copy of robustness weights to compare the next iteration with
      
      if (iter>iter_limit) stop("Error: iteration limit exceeded in robust estimation loop of RobustSolvAR")
      
    } # end robustness loop
    
    
    pacf <- wtd.pacf(e, shortwt)                     # we need to calculate *weighted* partial autocorrelations... using homebaked function above                             
    resid.pacf <- pacf$pacf
    resid.acf <- pacf$acf
    
    
    # Now calculate the parameter covariance matrix. Note that scaling needs to be correct here.  
    # If C is a covariance matrix (which it is), then s2e needs to be the variance of residuals (which it is).
    
    # Because we only want the uncertainties in the b's, treating the phi's as constant (which benchmark tests 
    # show is the right thing to do -- otherwise the standard errors are inflated), we need to truncate the covariance matrix C to exclude 
    # the rows and columns corresponding to the phi coefficients.
    
    if (h==0) clipC <- C
    else clipC <- C[-(mm+1):-(mm+h) , -(mm+1):-(mm+h)]                       # remove rows and columns corresponding to phi terms
    
    if (nu==0.0) Kbb = (s2e / n.eff) * solve(clipC)                          # covariance matrix of coefficients without regularization
    else {
      if (h==0) clipH <- H
      else clipH <- H[-(mm+1):-(mm+h) , -(mm+1):-(mm+h)]                     # remove rows and columns corresponding to phi terms
      regC.inverse = solve((clipC + (lambda * clipH)))                       # inverse of regularized C matrix
      Kbb = (s2e / n.eff) * regC.inverse %*% clipC %*% regC.inverse          # covariance matrix of regularized coefficients (equation A15 of KK2020)
    }
    
    
    
    # now shape output
    
    if (h==0) {    
      phi <- NA                                                              # if there is no AR correction, save phi as NA
    } else {                                                                 # if there is an AR correction...:
      phi <- b.hat[(length(b.hat)-h):(length(b.hat)-1)]                      # extract phi from fitted coefficients
    }
    
    
    b.hat <- matrix(as.vector(b.hat[1:mm]), ncol=nx)                         # wrap the b vector into a matrix (with one column for each x variable)
    b.hat <- b.hat[1:(m+1) , , drop=FALSE]                                   # truncate at maximum lag of m (since we don't need the other values any more)
    
    se <- matrix(sqrt(diag(Kbb))[1:mm], ncol=nx)                             # extract standard errors from diagonal of covariance matrix and wrap into a matrix with one column for each x
    se <- se[1:(m+1) , , drop=FALSE]                                         # truncate at maximum lag of m (since we don't need the other values any more)
    
    n.nz <- matrix(n.nz[1:mm], ncol=nx)                                      # wrap n.nz (from our original problem, not after AR correction) into a matrix
    n.nz <- n.nz[1:(m+1) , , drop=FALSE]                                     # truncate at maximum lag of m (since we don't need the other values any more)
    
    xKbb <- array(NA, dim=c(nx, 1+m, 1+m))                                   # create an array to export the covariance matrix
    for (j in 1:nx) {                                                        # to populate this array...
      ff <- (j-1)*(1+maxlag)+1                                               # we need to clip out the subset of the covariance matrix that pertains to each x
      ll <- ff+m
      xKbb[j, , ] <- Kbb[ff:ll , ff:ll] 
    }
    
    
    return(
      list(
        b.hat = b.hat ,    # regression coefficients (matrix of nx columns and m+1 rows)
        se = se ,           # standard errors (matrix of nx columns and m+1 rows)
        Kbb = xKbb ,         # stack of covariance matrices for each x (array of m+1 by m+1 matrices)
        n = n ,               # total number of time steps
        phi = phi ,            # fitted AR coefficients
        n.eff = n.eff ,         # effective number of time steps, accounting for uneven weighting
        n.nz = n.nz ,            # number of nonzero x values, with nonzero weight, in each column of design matrix (matrix of nx columns and m+1 rows)
        e = e ,                   # residuals 
        s2e = s2e ,                # weighted residual variance
        resid.acf = resid.acf ,     # autocorrelation function of residuals
        resid.pacf = resid.pacf ,    # partial autocorrelation function of residuals
        rw = rw                       # robustness weights
      )
    )
    
    
  } # end RobustSolvAR


#/////////////////////////////////////////////
# END OF robust solution to matrix equations #
#/////////////////////////////////////////////















#/////////////////////////////////////////////////////////////////
#### robust solution to matrix equations for broken stick IRF ####
#/////////////////////////////////////////////////////////////////
BrokenStickRobustSolvAR <-
  function(y, 
           x, 
           wt ,
           knots = c(0, 1, 2, 4, 8, 16, 32, 64) ,
           h = 0 ,
           nu = 0 ,
           robust = FALSE ,
           verbose = FALSE )
    
  {
    
    
    # takes as input:
    # y              vector of y values
    # x              matrix of x values
    # wt             vector of weight values
    # knots          vector of lags at knots for which piecewise linear IRF will be evaluated (in integer number of time steps)
    #                          must be positive integers in ascending order, with first value of zero (these are checked)
    # h              order of optional AR correction
    # nu             fractional weight, 0-1, to be given to Tikhonov-Phillips regularization (default is 0 = no regularization)
    # robust         flag controlling whether robust estimation by Iteratively Reweighted Least Squares (IRLS) will be used
    # verbose        controls whether progress reports are printed (TRUE) or suppressed (FALSE)
    
    # returns a list as output, with the following objects:
    # b.hat          regression coefficients (matrix of nx columns and nk rows)
    # se             standard errors (matrix of nx columns and nk rows)
    # Kbb            stack of covariance matrices of coefficients for each x (array of nk by nk matrices)
    # n              total number of time steps
    # n.eff          effective number of time steps, accounting for uneven weighting
    # n.nz           number of nonzero elements in the *weighted* x's (so values with zero weight don't count), stored as matrix of nx columns and m+1 rows
    # e              residuals
    # s2e            weighted residual variance
    # resid.acf      autocorrelation function of residuals
    # resid.pacf     partial autocorrelation function of residuals
    # rw             vector of robustness weights
    
    
    # Note that in order to keep the code relatively simple, the broken-stick solvAR does not compile the cross-products in chunks, but instead creates the
    # entire xx matrix all at once.  Thus it requires enough memory (either built-in or virtual) for this matrix (and potentially several copies of it).
    # On Windows, use memory.limit(size=###) to change the memory available to R, where ### is in MB. 
    
    
    iter_limit <- 200               # iteration limit for robust estimation
    
    
    if (min(knots) < 0) stop("Fatal error in BrokenStickRobustSolvAR: knots cannot be negative")
    if (any(knots != trunc(knots))) stop("Fatal error in BrokenStickRobustSolvAR: knots must be integer-valued")
    if (min(diff(knots)) <= 0) stop("Fatal error in BrokenStickRobustSolvAR: knots must be in ascending order, with no duplicates")
    if (knots[1] != 0) stop("Fatal error in BrokenStickRobustSolvAR: first knot must be zero")
    
    nk <- length(knots)             # number of knots
    
    maxlag <- knots[nk]+h          # maximum lag (not counting one for lag zero)
    
    
    #///////////////////////////////////////////////////
    # build list of filter vectors for each of the knots  (Equation 51 of K2022)
    
    # we will do this by brute force to make the underlying procedure explicit
    tauwt <- list()                                        # declare a list
    tauwt[[1]] <- 1                                        # start with the identity filter for lag zero
    for (k in 2:(nk-1)) {                                  # step through the other knots
      ww <- rep(0, knots[k+1])                             # note that because lags begin at zero (index=1 means lag=0) this vector will end at one lag *before* the next knot (which is OK)
      # ramp up from knots[k-1] to (and including) knots[k]
      for (i in knots[k-1]:knots[k]) ww[i+1] <- (i - knots[k-1])/(knots[k] - knots[k-1]) 
      
      # ramp down from knots[k] to knots[k+1] (not including either end, because knots[k] has already been done, and knots[k+1] would have a value of zero anyhow)
      if (knots[k+1]-knots[k] > 1) for (i in (knots[k]+1):(knots[k+1]-1)) ww[i+1] <- (knots[k+1] - i)/(knots[k+1] - knots[k]) 
      tauwt[[k]] <- ww                                     # add to list
    }
    
    ww <- rep(0, knots[nk]+1)                              # last filter vector
    for (i in knots[nk-1]:knots[nk]) ww[i+1] <- (i - knots[nk-1])/(knots[nk] - knots[nk-1]) 
    tauwt[[nk]] <- ww                                      # and add it to the list
    
    # end of filter vector construction
    
    
    nx <- NCOL(x)                   # this is the number of different x variables
    x <- as.matrix(x)               # convert x to a matrix if it isn't one
    
    n.terms <- nk+h                 # n.terms is the number of coefficients we need for each of the nx.
    
    mm <- n.terms*nx                # number of columns we will need in xx (plus h lags of y, plus 1 for the constant term)
    
    
    
    #///////////////////////////////////
    # zero out missing values of y and x
    
    orig_y <- y                                               # first make copies of y and x (for calculating residuals later)
    orig_x <- x
    orig_wt <- wt
    
    # note that all of these cases, and their lags, correspond to rows where wt has already been zeroed (in IRF), 
    # so replacing them with zeroes will not affect the values of the cross products (but will keep crossprod from choking on NA's)
    y[is.na(y)] <- 0                                          # replace missing values of y with zeroes
    x[is.na(x)] <- 0                                          # replace missing values of x with zeroes
    
    
    n <- length(y)                                            # number of rows
    
    xx <- matrix(0, nrow=n, ncol=mm+h+1)                      # design matrix with filtered time series for each knot and each x, plus AR terms and the constant term
    
    
    for (k in 1:nk) {
      xf <- stats::filter(x, filter=tauwt[[k]], method="c", sides=1)  # now apply convolutional filter for each knot to x's (xf is filtered x) -- note that doing this by FFT is MUCH slower!!
      for (i in 1:nx) xx[ , (i-1)*n.terms+k] <- xf[ , i]              # and distribute columns to the correct columns of xx
    }
    
    
    if (h>0) for (i in 1:nx) for (j in 1:h) xx[ , (i-1)*n.terms+nk+j] <- dplyr::lag(x[,i], knots[nk]+j)  # add one additional lag of each x for each AR term
    
    if (h>0) for (j in 1:h) xx[ , mm+j] <- dplyr::lag(y,j)     # copy AR terms, if any, into design matrix
    
    xx[ , ncol(xx)] <- 1                                       # don't forget the constant term...
    
    xx[1:maxlag, ] <- 0                                    # zero out any rows that will have NAs in filtered columns (these rows will all have zero weight anyhow)
    
    if (nu>0) {                                            # if we regularize...:
      # We can't use one big regularization matrix for the whole problem, because it will lead to leakage between the end of one set of b coefficients and the start of another. 
      # Instead we need to create a segmented regularization matrix, with nx Tikhonov-Phillips matrices along the diagonal
      tt <- tikh(n.terms)                                  # start with one Tikhonov-Phillips matrix (this is our building block)
      H <- matrix(0, nrow = mm+h+1, ncol=mm+h+1)           # start with a matrix of all zeroes
      for (j in 1:nx)                                      # copy our building block into the right range(s) along the diagonal of H
        H[((j-1)*n.terms+1):(j*n.terms) , ((j-1)*n.terms+1):(j*n.terms) ] <- tt
    }
    
    
    
    rw <- rep(1, n)                                          # robustness weights are initialized to 1
    old_rw <- rw                                             # save old robustness weights for comparison later
    
    #////////////////////////////
    # robustness loop starts here
    iter <- 0                                                # iteration counter
    repeat{
      
      iter <- iter+1                                         # update counter
      
      wt <- rw*orig_wt                                       # multiply weights by robustness weight
      sumwt <- sum(wt, na.rm=TRUE)
      sw <- as.vector(sqrt(wt))                              # sw is square root of weight (including robustness weights)
      
      # str(rw)
      
      C <- crossprod(xx*sw)                                  # use R's very fast cross-product routine (fast with a good BLAS, that is)
      xy <- crossprod(xx*sw, y*sw)                           # note that because x and yy are multiplied by sqrt(wt), these cross-products are weighted
      
      C <- C/sumwt                                           # dividing by the sums of the weights does not change the regression parameters, 
      xy <- xy/sumwt                                         # but is necessary for consistency with the way uncertainties are calculated next
      
      if (nu==0) b.hat <- solve(C, xy)                       # if we *don't* regularize, solve for the un-regularized regression coefficients
      else {                                                 # if we *do* regularize...:
        # convert weighting parameter nu (values 0 - 1) to lambda as used in conventional Tikhonov-Phillips regularization formulas
        lambda <- nu / (1.0 - nu) * sum(diag(C)[1:(nx*(1+n_lags))]) / sum(diag(H))                  # equation 50 of K2019 (omitting the diagonal elements that correspond to the constant term and AR terms)
        b.hat <- solve((C + (lambda * H)), xy) # now solve (regularized) linear system              # equation 46 of K2019, and equation A12 of KK2020
      }  
      
      # calculate predicted values -- note that rows with NA's will not count anyhow because they will have weight of zero in calculating error variance
      # need to use original x and y (which include NAs)
      pred <- rep(tail(b.hat,1), n)                          # initialize predicted values with the constant (the last element of b.hat)
      
      # now we need to create a matrix of beta values for each lag, linearly interpolated among the nk knots, for each of the x variables
      allb <- matrix(NA, nrow=knots[nk]+1, ncol=nx)
      
      for (j in 1:nx) {   # now step through each of the x's
        
        bk <- b.hat[((j-1)*n.terms+1):((j-1)*n.terms+nk), 1]     # subset the knot values of b for this x variable
        
        for (k in 2:nk) { 
          for (i in knots[k-1]:knots[k]) allb[i+1, j] <- bk[k-1]*(knots[k]-i)/(knots[k]-knots[k-1]) + bk[k]*(i-knots[k-1])/(knots[k]-knots[k-1])  # linear interpolation
        }
        
        pred <- pred + stats::filter(x[, j], filter=allb[ ,j], method="c", sides=1)  # use these interpolated b's to predict y (and keep them for later, for AR correction of beta...)
      } # next j
      
      if (h>0) for (j in 1:h) pred <- pred + dplyr::lag(y,j)*b.hat[mm+j] # include AR terms in predicted values
      
      
      e <- orig_y - pred                                                       # these are the residuals  (equation A13 of KK2020)
      
      n.nzwt <- sum((wt != 0), na.rm=TRUE)                                     # effective length of the x-matrix
      sumwt <- sum(wt)
      
      e2 <- e*e                                                                # squared errors
      
      s2e <- (n.nzwt-1)/(n.nzwt-(mm+h+1)-1) * sum( e2 * wt, na.rm=TRUE)/sumwt  # weighted residual variance, corrected for degrees of freedom
      
      if (robust==FALSE) {break}                                               # if we're not doing robust estimation, we can jump out of the loop here
      
      # now do iteratively reweighted least squares
      med.e <- sqrt(median(e2[wt!=0], na.rm=TRUE))                                                          # median absolute residual
      if (med.e==0) w <- ifelse(x==0, 1, 0)                                                                 # an unlikely exotic case: if MAR is zero, only 0 or 1 weights are possible
      else rw <- ifelse(is.na(e), 1, 1/(1+(e/(3.536*med.e))^2))                                             # Cauchy weight function
      

      # exit criterion: 99 percent of the robustness weights have changed by less than 0.1 (absolute value) since previous iteration
      if (quantile(abs(rw-old_rw)[wt!=0], 0.99)<0.1) {  # don't count rw's and residuals for rows that have no weight!!
        if (verbose) cat("robustness iteration ", iter, " complete at time ", round(difftime(Sys.time(), TICK, units="secs"), 3),
                         "  s2e=", round(s2e, 8), " max e=", round(sqrt(max(e2, na.rm=TRUE)), 3),
                         "rw chg= ", round(quantile(abs(rw-old_rw)[wt!=0], 0.99), 3), "\n")
        {break}                            
      }
      
      if (verbose) cat("robustness iteration ", iter, " complete at time ", round(difftime(Sys.time(), TICK, units="secs"), 3),
                       "  s2e=", round(s2e, 8), " max e=", round(sqrt(max(e2, na.rm=TRUE)), 3),
                       "rw chg= ", round(quantile(abs(rw-old_rw)[wt[(maxlag+1):n]!=0], 0.99), 3), "\n")
      
      old_rw <- rw
      
      if (iter>iter_limit) stop("Error: iteration limit exceeded in robust estimation loop of BrokenStickRobustSolvAR")
      
    } # end robustness loop
    
    
    
    resid.pacf <- as.vector(pacf(e, plot=FALSE, na.action=na.pass)$acf)      # residual ACF and PACF
    resid.acf <- as.vector(acf(e, plot=FALSE, na.action=na.pass)$acf)
    
    n.nzwt <- sum((wt[!is.na(e)] != 0), na.rm=TRUE)                          # effective length of the x-matrix (even though we've avoided creating one here...)
    sumwt <- sum(wt[!is.na(e)], na.rm=TRUE)
    
    s2e <- (n.nzwt-1)/(n.nzwt-(mm+h+1)-1) * sum(e * e * wt, na.rm=TRUE)/sumwt   # weighted residual variance, corrected for degrees of freedom
    
    n.eff <- sumwt^2/sum(wt*wt)                                              # calculate the effective sample size (accounting for uneven weighting)
    
    
    # Now calculate the parameter covariance matrix. Note that scaling needs to be correct here.  
    # If C is a covariance matrix (which it is), then s2e needs to be the variance of residuals (which it is).
    
    # Because we only want the uncertainties in the b's, treating the phi's as constant (which benchmark tests 
    # show is the right thing to do -- otherwise the standard errors are inflated), we need to truncate the covariance matrix C to exclude 
    # the rows and columns corresponding to the phi coefficients.
    
    if (h==0) clipC <- C
    else clipC <- C[-(mm+1):-(mm+h) , -(mm+1):-(mm+h)]                       # remove rows and columns corresponding to phi terms
    
    if (nu==0.0) Kbb = (s2e / n.eff) * solve(clipC)                          # covariance matrix of coefficients without regularization
    else {
      if (h==0) clipH <- H
      else clipH <- H[-(mm+1):-(mm+h) , -(mm+1):-(mm+h)]                     # remove rows and columns corresponding to phi terms
      regC.inverse = solve((clipC + (lambda * clipH)))                       # inverse of regularized C matrix
      Kbb = (s2e / n.eff) * regC.inverse %*% clipC %*% regC.inverse          # covariance matrix of regularized coefficients (equation A15 of KK2020)
    }
    
    
    
    # now shape output
    
    if (h==0) {
      phi <- NA                                                              # if there is no AR correction, save phi as NA
    } else {                                                                 # if there is an AR correction...:
      phi <- b.hat[(length(b.hat)-h):(length(b.hat)-1)]                      # extract phi from fitted coefficients
    }
    
    n.nz <- colSums(((xx*wt) != 0))                                          # count non-zero elements (with non-zero weight) in each x column
    
    
    b.hat <- matrix(as.vector(b.hat[1:mm]), ncol=nx)                         # wrap the b vector into a matrix (with one column for each x variable)
    b.hat <- b.hat[1:nk , , drop=FALSE]                                      # truncate at maximum knot index (since we don't need the other values any more)
    
    se <- matrix(sqrt(diag(Kbb))[1:mm], ncol=nx)                             # extract standard errors from diagonal of covariance matrix and wrap into a matrix with one column for each x
    se <- se[1:nk , , drop=FALSE]                                            # truncate at maximum knot index (since we don't need the other values any more)
    
    n.nz <- matrix(n.nz[1:mm], ncol=nx)                                      # wrap n.nz (from our original problem, not after AR correction) into a matrix
    n.nz <- n.nz[1:nk , , drop=FALSE]                                        # truncate at maximum knot index (since we don't need the other values any more)
    
    xKbb <- array(NA, dim=c(nx, nk, nk))                                     # create an array to export the covariance matrix
    for (j in 1:nx) {                                                        # to populate this array...
      ff <- (j-1)*n.terms+1                                                  # we need to clip out the subset of the covariance matrix that pertains to each x
      ll <- ff+nk-1
      xKbb[j, , ] <- Kbb[ff:ll , ff:ll] 
    }
    
    
    return(
      list(
        b.hat = b.hat ,   # regression coefficients (matrix of nx columns and m+1 rows)
        se = se ,          # standard errors (matrix of nx columns and nk rows)
        Kbb = xKbb ,        # stack of covariance matrices for each x (array of nk by nk matrices)
        allb = allb ,        # matrix of interpolated regression coefficients for each lag between 0 and m
        n = n ,               # total number of time steps
        phi = phi ,            # fitted AR coefficients
        n.eff = n.eff ,         # effective number of time steps, accounting for uneven weighting
        n.nz = n.nz ,            # number of nonzero x values, with nonzero weight, in each column of design matrix (matrix of nx columns and nk rows)
        e = e ,                   # residuals 
        s2e = s2e ,                # weighted residual variance
        resid.acf = resid.acf ,     # autocorrelation function of residuals
        resid.pacf = resid.pacf ,    # partial autocorrelation function of residuals
        rw = rw                       # vector of robustness weights
      )
    )
    
    
  } # end BrokenStickRobustSolvAR

#//////////////////////////////////////////////////////////
# END OF robust solution to broken stick matrix equations #
#//////////////////////////////////////////////////////////
















#////////////////////////////////////////////////////////////////////////////////////////////////
#### IRF - Impulse Response Function estimates by least squares with correction for AR noise ####
#////////////////////////////////////////////////////////////////////////////////////////////////


IRF <-
  function(y ,
           x ,
           xnames = NULL ,
           wt = rep(1, NROW(x)) ,
           m = 60 ,
           nk = 0 ,
           nu = 0 ,
           FD = FALSE ,
           h = NULL ,
           ARprob = 0.05 ,
           ARlim = 0.2 ,
           max.AR = 12 ,
           complete = FALSE ,
           verbose = FALSE ,
           robust = FALSE ,
           max.chunk = 2e8 )

  {
    
    # Estimates the impulse response of y to one or more x's over lags from 0 to m, with optional case weights "wt", while correcting for AR noise of arbitrary order.
    # Tikhonov-Phillips regularization can be optionally applied by setting nu>0.
    
    # takes as input:
    # y               a vector or time series of a single response variable.  May contain NA's.
    #
    # x               a vector, matrix, or time series containing one or more inputs of same length as y.  May contain NA's.
    #
    # xnames          optional vector of strings that name each input.  Length must equal the number of columns of x
    #
    # wt              optional vector of case weights of same length as y.  May contain NA's.
    #
    # m               maximum lag, in number of time steps.  Number of lags will be m+1 (the extra 1 is for lag zero).
    #
    # nk              number of knots in piecewise linear broken stick representation of beta as a function of lag.
    #                           Must be integer greater than 2 and less than or equal to m+1, or must be zero (the default).  
    #                           If nk>2 and nk<=m+1, nk knots will be created at lags of 0, 1, and a geometric progression 
    #                           (or as close to geometric as it can be, given that the lags are integers) between lags 1 and m.
    #                           If nk<=2, nk>m+1, or nk==0 (the default), the broken-stick approach is not used, and instead 
    #                           the IRF is evaluated for m+1 lags between 0 and m.
    #
    # nu              fractional weight, 0<=nu<1, to be given to Tikhonov-Phillips regularization (0 = no regularization)
    #
    # FD              flag for implementation of first differencing.  If FD==FALSE, no first differencing is applied.  
    #                          If FD==TRUE, y is first-differenced but x is not, and the coefficients and their covariance 
    #                          matrices are adjusted accordingly (as explained in K2022).
    #
    # h               integer order of autoregressive correction (non-negative integer).  AR corrections of sufficiently high order can, 
    #                          by the duality principle, be used to account for moving average (MA) noise as well.  
    #                          The value of h should be much less than m; otherwise strange results may occur.  If h==0, no correction is applied.
    #
    #                          If h==NULL, the order of autoregressive correction will be determined automatically, based on both statistical significance
    #                          and practical significance.  The practical significance criterion tests whether the absolute values of all of the acf
    #                          and pacf coefficients are less than the user-specified threshold ARlim, below which they are assumed to be practically
    #                          insignificant (that is, to have no substantial effect on the estimates of the IRFs).  The statistical
    #                          significance criterion tests whether the acf and pacf coefficients are statistically distinguishable from white noise
    #                          at the specified significance level ARprob.  First abs(acf) and abs(pacf) are compared at each lag against the 
    #                          two-tailed p<0.05 critical value for correlation coefficients, 1.96/sqrt(n.eff), where n.eff is the effective sample size
    #                          that accounts for uneven weighting.  Then the number of cumulative exceedances at each lag L (e.g., the number of cases 
    #                          where abs(acf)>1.96/sqrt(n.eff) at lags from 1 to L) is be compared to the critical number of exceedances predicted by 
    #                          binomial statistics (qbinom(p=ARprob, size=L, prob=0.05)).  The residuals are considered white if the cumulative number of 
    #                          exceedances is less than or equal to this critical number, over all lags, in both the acf and the pacf.  h is increased
    #                          sequentially until either the practical significance criterion or the statistical significance criterion is met, 
    #                          or until h==max.AR (which triggers a warning).  The practical significance test is needed because in large samples, 
    #                          even trivially small acf and pacf coefficients may still be statistically significant, triggering a pointless effort to
    #                          make them smaller than they need to be.  Conversely, the statistical significance test is needed because in small samples,
    #                          even true white noise processes may yield acf and pacf coefficients that do not meet the practical significance threshold 
    #                          (that is, ARlim may be less than 1.96/sqrt(n.eff)), triggering a pointless effort to further whiten residuals that are 
    #                          already white.
    #
    # ARprob          significance threshold for testing whether residuals are sufficiently white in automatic AR selection
    #
    # ARlim           threshold value of acf and pacf coefficients of residuals, used in practical significance test
    #
    # max.AR          maximum order of autoregressive correction that will be accepted in automatic AR order selection
    #
    # complete        flag for whether the number of lags will be *assumed* to be sufficient to hold the complete IRF (TRUE), meaning that any IRF
    #                          coefficients at longer lags are trivially small, or whether this cannot be assumed (FALSE).  Complete==TRUE will yield
    #                          smaller, and more accurately estimated, standard errors if the real-world IRF actually does converge to nearly zero before 
    #                          the maximum lag is reached.  But if this is not the case, complete==TRUE will force the IRF to artifactually
    #                          converge toward zero at the longest lag (with artificially small standard errors).  Complete==TRUE should thus
    #                          be invoked with caution.
    #
    # verbose         controls whether progress reports are printed (TRUE) or suppressed (FALSE)
    #
    # robust         flag controlling whether robust estimation by Iteratively Reweighted Least Squares (IRLS) will be used
    #
    # max.chunk       maximum size, in bytes, of the largest piece of the regression matrix (xx, the matrix of x and its lags) that will be created
    #                          at one time.  Such xx matrices that would be larger than max.chunk will be created and processed
    #                          in separate "chunks" to avoid triggering memory paging, which could substantially increase runtime,
    #                          or to avoid exceeding the available virtual memory, which will lead to a crash.
    #                          Keeping chunk.max relatively small (order 1e8 or 1e9) incurs only a small performance penalty, *if* 
    #                          one does not perform robust estimation.  But in robust estimation, the xx matrix can be re-used if it is
    #                          built in one piece, whereas if it is chunked it will need to re-built several times.  Therefore users should
    #                          not set max.chunk much smaller than it really needs to be.
    #                          Setting max.chunk=NULL suppresses chunking entirely.
    
    # y, x, and wt can contain missing values.  Missing values of y create 1+h missing rows in the regression matrix.
    #                          Missing values of x create m+1 missing rows (one for each lag) in the regression matrix.
    
    # returns a list as output, with the following objects
    #
    # lags        vector of lags (in number of time steps)
    #
    # IRF         impulse response function for each explanatory variable x (matrix of ncol(x) columns and m+1 rows, corresponding to lags 0 through m)
    #
    # se          standard errors of IRF coefficients (matrix of ncol(x) columns and m+1 rows, corresponding to lags 0 through m)
    #
    # Kbb         stack of covariance matrices for each x (array of m+1 by m+1 matrices, one for each column of x)
    #
    # n           length of original x and y series
    #
    # n.eff       effective sample size, accounting for uneven weighting
    #
    # n.nz        number of nonzero values (that also have nonzero weight) in each explanatory variable x at each lag  (matrix of ncol(x) columns and m+1 rows, corresponding to lags 0 through m)
    #
    # h           order of AR correction that was applied
    #
    # phi         fitted AR coefficients (vector of length h)
    #
    # resid       residuals (vector of length n)
    #
    # resid.acf   autocorrelation function of residuals
    #
    # resid.pacf  partial autocorrelation function of residuals
    
    
    
    
    #//////////////////////////////////////
    # do preliminary checks and conversions
    
    if (verbose) TICK <- Sys.time()
    
    if (is.ts(x)) {                                  # if x is a time series, coerce it to a vector or matrix, as appropriate
      if (is.vector(x)) x <- as.numeric(x)           # a one-dimensional time series is easy
      else x <- matrix(as.numeric(x), ncol=NCOL(x))  # otherwise we need to coerce to a vector and then reassemble the matrix
    }
    
    if (is.ts(y)) y <- as.numeric(y)                 # if y is a time series, coerce it to a vector
    
    if ((!is.vector(x)) & (!is.matrix(x))) stop("Fatal error in IRF: x must either be a matrix or a vector")
    
    nx <- NCOL(x)                                    # nx is the number of x variables
    
    n <- NROW(x)                                     # n is number of time steps (NROW works for both vectors and matrices)
    if (length(y)!=n |
        length(wt)!=n ) stop("Fatal error in IRF: input series x, y, and wt are not all equal in length")
    
    if( nu<0 | nu>=1 ) stop("Fatal error in IRF: nu must be >=0 and <1")
    
    if( (round(m)!=m) | m<0 ) stop("Fatal error in IRF: m must be a non-negative integer")
    
    if ((complete!=TRUE) & (complete!=FALSE)) stop("Fatal error in IRF: complete must be TRUE or FALSE")
    
    if (!is.null(h))
      if ((round(h)!=h) | (h<0)) stop("Fatal error in IRF: h must be a non-negative integer, or NULL")
    
    if (ARprob<=0 | ARprob>=1 ) stop("Fatal error in IRF: ARprob must be between 0 and 1")
    
    if (max.AR>20) warning("Large max.AR in IRF")
    
    if ((FD!=TRUE) & (FD!=FALSE)) stop("Fatal error in IRF: FD must be logical variable")
    
    if(round(nk)!=nk) stop("Fatal error in IRF: nk must be an integer")
    if (nk!=0) {
      if (nk>=(m+1)) {
        nk <- 0
        warning("nk too big -- reverting to evenly spaced lags in IRF")
      } else if (nk<3) {
        nk <- 0
        warning("nk too small -- reverting to evenly spaced lags in IRF")
      } else {
        if (FD==TRUE) {
          warning("first differencing turned off for broken-stick IRFs")
          FD <- FALSE
        }
      }
    }
    
    x <- as.matrix(x)                                # this converts x into a (potentially one-column) matrix, so that what follows will mostly not need to handle nx==1 and nx>1 as separate cases
    
    orig_y <- y                                      # save a copy of the original values of y, in case we do first-differencing
    
    if (FD==TRUE) {                                  # here we do our first-differencing if needed
      y <- y - dplyr::lag(y)                         # equation 28 of K2022
    }
    
    
    #//////////////////////////////////
    # here we build the lag knot series (if nk>0)
    
    if (nk>0) {
      # knot series should start as 0, 1, and then follow geometric progression to reach m at the nk'th knot, but with a minimum step size of 1 (no duplicate knots!)
      knots <- rep(NA, nk)                             # declare knot vector
      knots[1] <- 0
      knots[2] <- 1
      for (k in 3:nk) {
        ratio <- (m/knots[k-1])^(1/(nk-(k-1)))             # ratio for a geometric progression that will reach value of m by knot number nk
        if ((ratio-1)*knots[k-1]>1) knots[k] <- trunc(knots[k-1]*ratio)    # follow geometric progression if this would be a step bigger than 1
        else knots[k] <- knots[k-1] + 1                                # otherwise take a step of 1
      }
    } else {
      knots <- NULL
    }
    
    
    # determine whether we will need to do chunking or not
    byte.size <- 8       # nominal size (in bytes) of a numeric variable in R (not counting the one-time overhead, for example 48 bytes for any vector)
    # No need to change byte.size unless this fundamental property of R changes.  This is only used for calculating chunk sizes in solvAR.
    if (is.null(max.chunk)) chunk.maxht <- NULL else {
      if (nk>0) chunk.maxht <- floor(max.chunk / (nk*nx* byte.size))  # maximum height of a chunk of xx matrix, to be passed to solvAR
      else chunk.maxht <- floor(max.chunk / ((m+1)*nx* byte.size))  # maximum height of a chunk of xx matrix, to be passed to solvAR
    } 
    
    
    #////////////////////////////////////
    # here we start the AR selection loop
    # if autoAR==TRUE (i.e., h=NULL), this loop incrementally increases h until pacf.limit is met, or max.AR is reached
    # if autoAR==FALSE (i.e., h=integer), this loop runs just once at the user-supplied AR order
    
    
    if (is.null(h)) {
      h <- 0
      autoAR <- TRUE
    } else autoAR <- FALSE
    
    firstpass <- TRUE                  # flag for whether we are taking our first pass through the AR selection loop
    
    repeat {                           # start the AR selection loop
      
      n_lags <- m+h                    # total number of lags we will need (not counting lag zero)
      if (FD) n_lags <- n_lags+1
      
      # the strategy for handling NA's is as follows: we assign a weight of 0 to any row that will contain NA's (also for any lagged variables)
      # and in solvAR we set NA's to zero (so the crossprod function doesn't choke), but these will correspond to rows with zero weight so the
      # numeric values will be preserved
      
      wt[1:n_lags] <- 0                                               # assign weight of zero to first n_lags rows (because these will contain NAs)
      wt[is.na(wt)] <- 0                                              # assign weight of zero to rows where weight is missing
      
      # if a row is missing, assign weight of zero to that row and the next n_lags rows (total of 1+n_lags rows)
      # this is the non-vector way we did it previously: for (i in 1:n) if (missing.row[i]) wt[i:min(n, i+n_lags)] <- 0 # need to do this because those rows would contain a lagged NA value of x
      missing.row <- as.integer(is.na(rowSums(x, na.rm=FALSE)))       # boolean vector of rows that have NA in x
      if (sum(missing.row)>0) {
        missing.cancels <- rep(1, 1+n_lags)                           # this is a mask that propagates the effects of missing x's forward n_lags steps
        missing.row <- as.vector(stats::filter(missing.row, missing.cancels, sides=1, circular=TRUE))  # propagate the mask
        wt[(missing.row>0.5)] <- 0                                    # and set weight for masked rows to zero (using >0.5 here in case we have any near-zero values.  We should have only integers, but let's not take chances.)
      }
      
      # if a value of y is missing, assign weight of zero to that row and the next h rows (total of 1+h rows) because AR terms (lagged y's) will be NA at those rows
      # this is the non-vector way we did it previously:  for (i in 1:n) if (is.na(y[i])) wt[i:min(n, i+h)] <- 0        # need to do this because those rows would contain a lagged NA value of y
      missing.y <- as.integer(is.na(y) )                              # boolean vector of rows that have NA in y
      if (sum(missing.y)>0) {
        if (h==0) wt[(missing.y>0)] <- 0
        else {
          missing.cancels <- rep(1, 1+h)                              # this is a mask that propagates the effects of missing y's forward h steps
          missing.row <- as.vector(stats::filter(missing.y, missing.cancels, sides=1, circular=TRUE))  # propagate the mask
          wt[(missing.row>0.5)] <- 0                                    # and set weight for masked rows to zero
        }
      }
      
      
      #//////////////////////////////////////////////
      # now do cross-products and solve linear system
      
      if (robust & firstpass) {  
        
        # If we are doing robust estimation, we first need to determine robustness weights, setting h=0 because robust estimation doesn't play well with AR correction.
        # Note that we need h=0 here even if we subsequently are using a different level of AR correction.
        if (nk>0) s1 <- BrokenStickRobustSolvAR(y=y, x=x, wt=wt, knots=knots, h=0, nu=nu, robust=robust, verbose=verbose)
        else s1 <- RobustSolvAR(y=y, x=x, wt=wt, m=m, FD=FD, h=0, nu=nu, robust=robust, verbose=verbose, chunk.maxht=chunk.maxht)
        
        # now multiply weights by robustness weights, and use these in all further trips around this loop.
        # Don't change wt vector after this, except for zeroing out additional values as needed if h increases!
        wt <- wt*s1$rw  
        
        firstpass <- FALSE                                                 # now turn firstpass flag off, so we won't keep coming back here
        
        if ((autoAR==FALSE) & (h==0)) {
          list2env(s1, envir=environment())                                # unpack list from solvAR
          break                                                            # and jump out of loop          
        }
        
      }  # end if (robust & firstpass)
      
      
      # Here is where we do most of the work.  Note that we set robust=FALSE even if we are doing robust estimation, because we have already determined robustness weights above.
      if (nk>0) s1 <- BrokenStickRobustSolvAR(y=y, x=x, wt=wt, knots=knots, h=h, nu=nu, robust=FALSE, verbose=verbose)
      else s1 <- RobustSolvAR(y=y, x=x, wt=wt, m=m, FD=FD, h=h, nu=nu, robust=FALSE, verbose=verbose, chunk.maxht=chunk.maxht)
      
      list2env(s1, envir=environment())                                # unpack list from solvAR
      
      
      # s1 includes
      # b.hat          regression coefficients (matrix of nx columns and m+1 rows, or nk rows if knots!=NULL)
      # se             standard errors of regression coefficients
      # Kbb            stack of covariance matrices of regression coefficients for each x variable
      #                    (array of m+1 by m+1 matrices, or nk-by-nk matrices if knots!=NULL)
      # n              total number of time steps
      # phi            fitted AR coefficients
      # n.eff          effective number of time steps, accounting for uneven weighting
      # n.nz           a vector of mm+1 tallies of number of nonzero elements in the *weighted* x's (so values with zero weight don't count)
      # e              residuals
      # s2e            weighted residual variance
      # resid.acf      autocorrelation function of residuals
      # resid.pacf     partial autocorrelation function of residuals
      
      
      #///////////////////////////////////////////////
      # here we test for exiting the AR selection loop
      
      
      if (autoAR==FALSE) {break}                                                                 # if we weren't doing automatic AR order selection, exit here
      
      
      if (verbose) cat("at h=", h, " and time=", round(difftime(Sys.time(), TICK, units="secs"), 3), " residual PACF (first 5 terms) = ", round(resid.pacf[1:5], 4), "\n")
      
      if (verbose & (h>0)) {
        cat("AR coefficients (phi) = ", s1$phi, "\n")
      }
      
      # first we do the practical significance test
      
      if ( max(abs(resid.pacf), abs(resid.acf[-1])) < ARlim ) {                                  # if we pass the practical significance test...
        if (verbose) cat(paste0("practical (in)significance test passed at h=", h, "\n"))        # write a note
        {break}                                                                                  # and exit
      }
      
      
      # The statistical significance test inevitably involves multiple comparisons (multiple acf and pacf values to be compared to a threshold).
      # It does not use a Bonferroni-style alpha value to account for these multiple comparisons, because such small alphas might be vulnerable to distributional anomalies in the residuals.
      # Instead we tally the number of times that we exceed the p<0.05 threshold for individual acf and pacf estimates, and then test whether this number of exceedances
      # is improbable at a significance level of ARprob.  These comparisons are made sequentially: first we check whether the number of exceedances at the first lag is greater than
      # the expected number (which, if ARprob=0.05, is zero), then we check whether the number of exceedances in the first two lags is greater than the expected number (which, 
      # if ARprob=0.05, is one), then we check over the first three lags, and so on out to the length of the acf or pacf.  
      # We apply this criterion separately to the acf and pacf, because AR noise will cause many exceedances in the acf (but not the pacf), 
      # whereas MA noise will cause many exceedances in the pacf (but not the acf).
      
      threshold <- 1.96/sqrt(n.eff)                                                              # this is the (absolute) value of an *individual* acf or pacf element that would be *individually* statistically significant at p=0.05
      acf.ncrit <- qbinom(p=ARprob, size=seq(length(resid.acf)-1), prob=0.05, lower.tail=FALSE)  # running tally of the critical number of abs(acf)>acf.limit that we would expect to occur <ARprob fraction of the time 
      pacf.ncrit <- qbinom(p=ARprob, size=seq(length(resid.pacf)), prob=0.05, lower.tail=FALSE)  # running tally of the critical number of abs(pacf)>acf.limit that we would expect to occur <ARprob fraction of the time
      acf.exceedances <- cumsum(abs(resid.acf[-1])>threshold)                                    # tally the running total number of exceedances in the acf
      pacf.exceedances <- cumsum(abs(resid.pacf>threshold))                                      # tally the running total number of exceedances in the pacf
      
      if ((sum(acf.exceedances>acf.ncrit)==0) & (sum(pacf.exceedances>pacf.ncrit)==0)) {         # if we meet both of the statistical criteria...
        if (verbose) cat(paste0("statistical (in)significance test passed at h=", h, "\n"))      # write a note
        {break}                                                                                  # and exit
      }
      
      
      if (h >= max.AR) {                                                                         # if we have reached the maximum AR order allowed
        warning(paste0("Maximum AR order of ", max.AR, " reached in iterative loop. Consider aggregating time steps or first differencing."))            # print a warning
        {break}                                                                                  # and exit
      }
      
      h <- h + 1                                                                                 # if we haven't exited, increment h
      
    } # inner (AR selection) loop ends here                                                      # and continue around the AR selection loop again
    
    
    
    
    
    
    #///////////////////////////////////////////////////////////////////////////////////////////////
    # Now we correct the coefficients and their standard errors for AR noise
    # We use one of two different ways, depending on whether the IRF is assumed to be complete or not
    
    if (h>0) {                                                        # if h==0, then the first-stage results (s1) are already complete and we can go straight to the end
      
      #/////////////////////////////////////////////
      if (complete==FALSE) {                                          # if the IRF is not assumed complete, then...:
        
        # now we need to adjust the IRF coefficients and their standard errors to take account of the autoregressive noise
        # to do this we need to construct the psi matrix (the inverse of phi)
        
        v <- c(1, (-1*s1$phi), rep(0, m-h))                           # start with a vector whose first element is 1, followed by the phi values, and padded out with zeros to length m+1  (equation 20 of K2022)
        phi.mat <- toeplitz(v)                                        # create the corresponding Toeplitz matrix (symmetric, for now)
        phi.mat[upper.tri(phi.mat)] <- 0                              # and zero out the upper triangle... voila!
        psi <- solve(phi.mat)                                         # and take the inverse
        
        if (nk==0) {
          
          # if we have evenly spaced lags (conventional IRF, not broken-stick IRF)
          for (j in 1:nx) {                                             # for each x variable...
            b.hat[ , j] <- psi %*% b.hat[ , j]                          # convert b's using the inverse of the phi matrix   (equation 24 of K2022)
            Kbb[j, , ] <- psi %*% Kbb[j, , ] %*% t(psi)                 # use psi to propagate each covariance matrix (this is the standard way to do it, since psi is the Jacobian for the translation of the b's)  (equation 26 of K2022)
            se[ , j] <- sqrt(diag(Kbb[j, , ]))                          # and the square root of the diagonal gives the standard errors
          } 
          
        } else {
          
          # if we have unevenly spaced lags (broken-stick IRF)
          for (j in 1:nx) {                                             # for each x variable...
            allb[ , j] <- psi %*% allb[ , j]                            # convert b's using the inverse of the phi matrix   (equation 24 of K2022)
            b.hat[ ,j] <- allb[knots+1 , j]                             # sample the b's at the knot points to retrieve the AR-corrected beta values (need the +1 because lag zero is index=1...)
          }
          
          # Unfortunately the procedure above does not allow us to correct the standard error estimates for the autocorrelation in the noise.
          # For that, we will do something analogous to Cochrane-Orcutt, 
          # but without iterating to find the phi values -- we already have them
          new_y <- y
          for (i in 1:h) new_y <- new_y - s1$phi[i]*dplyr::lag(y, i)    # transform y to remove autoregressive errors (based on phi from first-stage results)
          
          new_x <- x                                                    
          for (i in 1:h) new_x <- new_x - s1$phi[i]*dplyr::lag(x, i)    # transform x to match the transformation of y                
          
          # now call BrokenStickSolvAR again, *without* AR terms, using AR-corrected x *AND* y, to get the standard error estimates (and new IRF coefficients, *if* we assume that the IRF is "complete"
          s2 <- BrokenStickRobustSolvAR(y=new_y, x=new_x, wt=wt, knots=knots, h=0, nu=nu, verbose=verbose, robust=FALSE)   # call solvAR again, now *without* AR terms, to get new IRF coefficients
          
          se <- s2$se                                                   # copy the standard errors from this run
          Kbb <- s2$Kbb                                                 # copy the parameter covariance matrices too
          
        } # end if (nk==0)
        
        #////////////////////////////////////////////        
      } else {                                                        
        
        # if we assume the IRF is complete, we will do something analogous to Cochrane-Orcutt, 
        # but without iterating to find the phi values -- we already have them
        new_y <- y
        for (i in 1:h) new_y <- new_y - s1$phi[i]*dplyr::lag(y, i)    # transform y to remove autoregressive errors (based on phi from first-stage results)
        
        new_x <- x                                                    
        for (i in 1:h) new_x <- new_x - s1$phi[i]*dplyr::lag(x, i)    # transform x to match the transformation of y                
        
        # now call solvAR again, *without* AR terms, using AR-corrected x *AND* y, to get new IRF coefficients
        if (nk>0) s2 <- BrokenStickRobustSolvAR(y=new_y, x=new_x, wt=wt, knots=knots, h=0, nu=nu, verbose=verbose, robust=FALSE)          # call solvAR again, now *without* AR terms, to get new IRF coefficients
        else s2 <- RobustSolvAR(y=new_y, x=new_x, wt=wt, m=m, FD=FD, h=0, nu=nu, verbose=verbose, robust=FALSE, chunk.maxht=chunk.maxht)  # call solvAR again, now *without* AR terms, to get new IRF coefficients
        
        list2env(s2, envir=environment())                             # unpack list from solvAR
        
      } # end else                                                   # and that's it!  we don't need to transform beta's or s.e.'s in this case
      
    } # end if (complete==FALSE)
    
    
    if (FD) {                                                         # if we are using first differences
      v <- rep(1, m+1)                                                # we need to construct the matrix that un-does the effects of first differencing on the b's  (equation 30 of K2022)
      FD.mat <- toeplitz(v)
      FD.mat[upper.tri(FD.mat)] <- 0
      
      for (j in 1:nx) {                                               # for each x variable...
        b.hat[ , j] <- FD.mat %*% b.hat[ , j]                         # convert b's using the inverse of the phi matrix  (equation 31 of K2022, with AR terms already accounted for in b.hat)
        Kbb[j, , ] <- FD.mat %*% Kbb[j, , ] %*% t(FD.mat)             # use psi to propagate each covariance matrix (this is the standard way to do it, since psi is the Jacobian for the translation of the b's) (equation 32 of K2022, with AR terms already accounted for in Kbb)
        se[ , j] <- sqrt(diag(Kbb[j, , ]))                            # and the square root of the diagonal gives the standard errors
      }
    }
    
    
    
    if (verbose) cat("IRF finished...:", round(difftime(Sys.time(), TICK, units="secs"), 3), "\n")
    
    if (verbose) cat(paste0("minimum column count of nonzero x values = ", min(n.nz), " (", round(100*min(n.nz)/max(n.nz)), "% of maximum count, and ", round(100*min(n.nz)/sum(wt>0), 1), "% of rows with nonzero weight)\n"))
    
    if (verbose & (h>0)) {
      cat("AR coefficients (phi) = ", s1$phi, "\n")
    }
    
    if (verbose) {
      cat("residual PACF (first 5 terms) =", round(resid.pacf[1:5], 4), "\n")
    }
    
    
    #//////////////////////////////////////////////
    # Now calculate predicted values of y from AR-adjusted coefficients, to see how well the statistical model fits the time-series behavior of y.
    # That may not be realistically reflected by the residuals in solvAR, because if y is strongly autocorrelated (and h>0), values of y will inherently
    # be close to their previous values.  This artifact is avoided if we use AR-adjusted coefficients and *only* x, and not lagged values of y,
    # to predict y.
    
    # first we make a naive estimate of y, without an intercept (because this has been lost in the AR correction)
    if (nk>0) ypred <- calc_ypred(x=as.vector(x), beta=as.vector(allb), nx=NCOL(x), n=NROW(x), m=m, h=(h+as.integer(FD)))
    else ypred <- calc_ypred(x=as.vector(x), beta=as.vector(b.hat), nx=NCOL(x), n=NROW(x), m=m, h=(h+as.integer(FD)))
    # then to add the intercept, we determine the weighted average difference between ypred and y, with robustness weights if robust==TRUE
    
    ww <- wt[(n_lags+1):n]
    ypred <- ypred + (weightedMean(orig_y[(n_lags+1):n], ww, na.rm=TRUE) - weightedMean(ypred, ww, na.rm=TRUE))
    ycomp <- data.table(timestep=(n_lags+1):n, y=orig_y[(n_lags+1):n], ypred=ifelse(ww>0, ypred, NA), yresid=ifelse(ww>0, orig_y[(n_lags+1):n]-ypred, NA), x=x[(n_lags+1):n,])
    
    # now create column names
    # if vector of xnames does not exist, create it
    if (is.null(xnames) | length(xnames)!=NCOL(x)) xnames <- paste0("x", 1:NCOL(x))
    
    colnames(b.hat) <- paste0("IRF_", xnames)
    colnames(se) <- paste0("se_", xnames)
    
    if (nk>0) lags <- knots 
    else lags <- 0:m
    
    return(
      list(
        lags = lags ,      # lags (in number of time steps)
        IRF = b.hat ,       # impulse response function for each column of x (matrix of ncol(x) columns and m+1 rows, corresponding to lags 0 through m)
        se = se ,            # standard errors (matrix of ncol(x) columns and m+1 rows, corresponding to lags 0 through m)
        Kbb = Kbb ,           # stack of covariance matrices for each x (array of ncol(x) m+1 by m+1 matrices)
        n = n ,                # length of original x and y series
        n.eff = s1$n.eff ,      # effective sample size, accounting for uneven weighting
        n.nz = s1$n.nz ,         # number of nonzero values with nonzero weight at each lag of each x variable
        h = h ,                   # order of AR correction that was applied
        phi = s1$phi ,             # fitted AR coefficients (vector of length h)
        resid = e ,                 # residuals
        resid.acf = resid.acf ,      # autocorrelation function of residuals
        resid.pacf = resid.pacf ,     # partial autocorrelation function of residuals
        ycomp = ycomp                  # data table comparing measured and fitted y time series (starting after initial lags)
        
      )
    ) # end return
    
    
    
  }  #end IRF

#///////////////////
# END of IRF
#///////////////////










#///////////////////////////////////////////////////////////////////////////////////////////////
#### nonlinIRF - nonlinear Impulse Response Function estimates with correction for AR noise ####
#///////////////////////////////////////////////////////////////////////////////////////////////


nonlinIRF <- function(y ,
                      x ,
                      xnames = NULL ,
                      wt = rep(1, NROW(x)) ,
                      m = 60 ,
                      nk = 0 ,
                      xknots = c(20, 50, 80, 90, 95) ,
                      pct_xknots = TRUE ,
                      nu = 0 ,
                      FD = FALSE ,
                      h = NULL ,
                      ARprob = 0.05 ,
                      ARlim = 0.2 ,
                      max.AR = 12 ,
                      complete = FALSE ,
                      verbose = FALSE ,
                      robust = FALSE ,
                      max.chunk = 2e8 )
{
  # Shell that calls IRF to estimate nonlinear impulse response functions, i.e., IRFs that depend on the input (x)
  
  
  # takes as input:
  # y               a vector or time series of a single response variable.  May contain NA's.
  #
  # x               a vector, matrix, or time series containing one or more inputs of same length as y.  Must be non-negative.  May contain NA's.
  #
  # xnames          optional vector of strings that name each input.  Length must equal the number of columns of x
  #
  # xknots          a vector or matrix of knots for the piecewise linear approximation of the RRD's nonlinear dependence on x.
  #                       Knots can be specified as fixed values or as percentiles of the x distribution (depending on the pct_knots flag -- see below).
  #                                  Values of p=0 are ignored when these percentiles are later converted to fixed values
  #                       If xknots is a matrix, it must have ncol(x), and each column of knots will be applied to the corresponding column of x.
  #                                  Each column of knots must have the same number of knots, although the knot values themselves may differ.
  #                       If xknots is a vector or single-column matrix and x is a multi-column matrix, the same knots will be applied to each column of x.
  #                       Knots must be between (and not include) the minimum and maximum values of the corresponding column of x.  
  #                       If xknots == NULL, potential nonlinear responses to input intensity are ignored and a single IRF is estimated for each 
  #                                 column of x.
  #                       If xknots != NULL, separate IRF's are estimated for each specified knot point (and also the maximum value) of each column of x.
  #
  # pct_xknots      a flag indicating whether nonlinearity knots are expressed values of x (FALSE) or percentiles (TRUE, 1, 2, or 3).
  #                       If pct_knots==TRUE or pct_xknots==1, knots are calculated as percentiles of x.
  #                       If pct_knots==2, knots are calculated as percentiles of the cumulative sum of x (so if xknots=80, for example, the knot will be the value of x 
  #                                 for which all smaller x's sum to 80 percent of the sum of all x's).  Thus these knots will delimit fractions of the total input x.
  #                       If pct_knots==3, knots are calculated as percentiles of the cumulative sum of squares of x (so if xknots=20, for example, the knot will be the value of x for which all smaller x's, 
  #                                 squared and summed, add up to 20 percent of the sum of squared x's).  Thus these knots will delimit fractions of the total sum of squared inputs x^2.
  #                                 These will roughly approximate the corresponding fractions of the total leverage in the data set, if the distribution of x's is strongly skewed with a peak near zero and a long right tail.
  #                       Any values other than TRUE, FALSE, 0, 1, 2, or 3 will trigger an error.
  #
  #                       knots should be chosen so that (these are not checked):
  #                       (a) there are enough points in each interval between knots, and sufficient variability in their values, to define the dependence of y on x,
  #                       (b) intervals between individual pairs of knots do not span major changes in slope in the nonlinear relationship between y and x
  #  
  #                       y, x, and wt can contain missing values.  Missing values of y create 1+h missing rows in the design matrix.
  #                       Missing values of x create m+1 missing rows (one for each lag) in the regression matrix.
  #
  #
  #
  # The remaining inputs are identical to those in IRF:
  #
  # wt              optional vector of case weights of same length as y.  May contain NA's.
  #
  # m               maximum lag, in number of time steps.  Number of lags will be m+1 (one extra for lag zero).
  #
  # nk              number of knots in piecewise linear broken stick representation of beta as a function of lag.
  #                           Must be integer greater than 2 and less than or equal to m+1, or must be zero (the default).  
  #                           If nk>2 and nk<=m+1, nk knots will be created at lags of 0, 1, and a geometric progression 
  #                           (or as close to geometric as it can be, given that the lags are integers) between lags 1 and m.
  #                           If nk<=2, nk>m+1, or nk==0 (the default), the broken-stick approach is not used, and instead 
  #                           the IRF is evaluated for m+1 lags between 0 and m.
  #
  # nu              fractional weight, 0 <= nu < 1, to be given to Tikhonov-Phillips regularization (0 = no regularization)
  #
  # FD              flag for implementation of first differencing.  If FD==FALSE, no first differencing is applied.  
  #                          If FD==TRUE, y is first-differenced but x is not, and the coefficients and their covariance 
  #                          matrices are adjusted accordingly (as explained in K2022).
  #
  # h               integer order of autoregressive correction (non-negative integer).  AR corrections of sufficiently high order can, 
  #                          by the duality principle, be used to account for moving average (MA) noise as well.  
  #                          The value of h should be much less than m; otherwise strange results may occur.  If h==0, no correction is applied.
  #
  #                          If h==NULL, the order of autoregressive correction will be determined automatically, based on both statistical significance
  #                          and practical significance.  The practical significance criterion tests whether the absolute values of all of the acf
  #                          and pacf coefficients are less than the user-specified threshold ARlim, below which they are assumed to be practically
  #                          insignificant (that is, to have no substantial effect on the estimates of the IRFs).  The statistical
  #                          significance criterion tests whether the acf and pacf coefficients are statistically distinguishable from white noise
  #                          at the specified significance level ARprob.  First abs(acf) and abs(pacf) are compared at each lag against the 
  #                          two-tailed p<0.05 critical value for correlation coefficients, 1.96/sqrt(n.eff), where n.eff is the effective sample size
  #                          that accounts for uneven weighting.  Then the number of cumulative exceedances at each lag L (e.g., the number of cases 
  #                          where abs(acf)>1.96/sqrt(n.eff) at lags from 1 to L) is compared to the critical number of exceedances predicted by 
  #                          binomial statistics (qbinom(p=ARprob, size=L, prob=0.05)).  The residuals are considered white if the cumulative number of 
  #                          exceedances is less than or equal to this critical number, over all lags, in both the acf and the pacf.  h is increased
  #                          sequentially until either the practical significance criterion or the statistical significance criterion is met, 
  #                          or until h==max.AR (which triggers a warning).  The practical significance test is needed because in large samples, 
  #                          even trivially small acf and pacf coefficients may still be statistically significant, triggering a pointless effort to
  #                          make them smaller than they need to be.  Conversely, the statistical significance test is needed because in small samples,
  #                          even true white noise processes may yield acf and pacf coefficients that do not meet the practical significance threshold 
  #                          (that is, ARlim may be less than 1.96/sqrt(n.eff)), triggering a pointless effort to further whiten residuals that are 
  #                          already white.
  #
  # ARprob          significance threshold for testing whether residuals are sufficiently white in automatic AR selection
  #
  # ARlim           threshold value of acf and pacf coefficients of residuals, used in practical significance test
  #
  # max.AR          maximum order of autoregressive correction that will be accepted in automatic AR order selection
  #
  # complete        flag for whether the number of lags will be *assumed* to be sufficient to hold the complete IRF (TRUE), meaning that any IRF
  #                          coefficients at longer lags are trivially small, or whether this cannot be assumed (FALSE).  Complete==TRUE will yield
  #                          smaller, and more accurately estimated, standard errors if the real-world IRF actually does converge to nearly zero before 
  #                          the maximum lag is reached.  But if this is not the case, complete==TRUE will force the IRF to artifactually
  #                          converge toward zero at the longest lag (with artificially small standard errors).  Complete==TRUE should thus
  #                          be invoked with caution.
  #
  # robust          flag controlling whether robust estimation by Iteratively Reweighted Least Squares (IRLS) will be used
  #
  # verbose         controls whether progress reports are printed (TRUE) or suppressed (FALSE)
  #
  # max.chunk       maximum size, in bytes, of the largest piece of the design matrix (the matrix of x and its lags) that will be created
  #                          at one time.  Design matrices that would be larger than max.chunk will be created and processed
  #                          in separate "chunks" to avoid triggering memory paging, which could substantially increase runtime.  
  #                          Keeping chunk.max relatively small (order 1e8 or 1e9) incurs only a small performance penalty, except in the case of robust estimation,
  #                          where the need to iterate means that the solution will be faster (by about a factor of 2 or so) *if* one can keep the whole design matrix
  #                          in one chunk without triggering memory paging.
  
  
  # returns a list as output, with the following objects
  #
  # knots       knots, in matrix of nk rows and nx columns (not including each x's first knot, which is zero)
  #
  # IRF         nonlinear impulse response function (beta) evaluated for each x variable at each knot point (except zero): matrix of nx*nk columns and m+1 rows, corresponding to lags 0 through m)
  #
  # se          standard errors of IRF coefficients (matrix of nx*nk columns and m+1 rows, corresponding to lags 0 through m)
  #
  # ykx         nonlinear impulse response, expressed as contribution to y from x, evaluated at each knot point (except zero): matrix of nx*nk columns and m+1 rows, corresponding to lags 0 through m)
  #
  # ykx_se      standard errors of ykx (matrix of nx*nk columns and m+1 rows, corresponding to lags 0 through m)
  #
  #             Note that ykx will show the shape of the nonlinearity more intuitively than IRF.
  #             For a linear system, ykx will be a straight line, whereas IRF will be a constant (within uncertainty)
  #             For a quadratic dependence, ykx will be parabolic, whereas IRF will be a straight line
  #
  # avg_IRF     arithmetic average of IRFs over all time steps with nonzero x's (use this when IRFs converge toward zero as x approaches zero)
  #
  # avg_se      standard error of avg_IRF
  #
  # alt_avg_IRF arithmetic average of IRFs by alternate method (for IRFs that do not converge to zero as x approaches zero)
  #
  # alt_avg_se  standard error of average IRF by alternate method
  #
  #             Note that arithmetic averages of IRFs will often depend strongly on the behavior of the IRF near x=0, since the distribution of x
  #             is strongly skewed with most values near zero and a long upper tail.  In such cases, the arithmetic average IRF will be close to
  #             the IRF for x near zero (which may be difficult to estimate accurately, because small values of x will yield small system responses).
  #
  #             Also in such cases, the average IRF will depend strongly on whether the IRF converges to 0, or to a finite nonzero value, as x approaches 0.
  #             Different methods are needed to yield reasonable estimates average IRFs in the two cases.  avg_IRF yields reasonable estimates for most
  #             nonlinear cases in which the IRF converges to zero as x approaches zero (note that this refers to the IRF=ykx/x, not ykx, which should 
  #             *always* converge toward zero!).  alt_avg_IRF yields reasonable estimates for most cases in which the IRF converges to a finite nonzero
  #             value as x approaches zero (one example would be a sublinear system, in which the IRF is finite near zero and declines with increasing x).
  #             
  #
  # avg_ykx     time-averaged ykx (including time steps with x=0)
  #
  # avg_ykx_se  standard error of avg_ykx
  #
  # wtd_avg_IRF weighted average of IRFs (weighted by input x)
  #
  # wtd_avg_se  standard error of wtd_avg_IRF
  # 
  # wtd_avg_ykx weighted average of ykx (weighted by the input x)
  #
  # wtd_avg_ykx_se  standard error of wtd_avg_ykx
  #
  #
  #
  # the rest of these outputs are passed from IRF
  #
  # Kbb         stack of covariance matrices for each xprime (array of nx*nk matrices, each m+1 by m+1)
  #
  # n           length of original x and y series
  #
  # n.eff       effective sample size, accounting for uneven weighting
  #
  # n.nz        number of nonzero values (that also have nonzero weight) in each explanatory variable x at each lag  (matrix of ncol(x)*nk columns and m+1 rows, corresponding to lags 0 through m)
  #
  # h           order of AR correction that was applied
  #
  # phi         fitted AR coefficients (vector of length h)
  #
  # resid       residuals (vector of length n)
  #
  # resid.acf   autocorrelation function of residuals
  #
  # resid.pacf  partial autocorrelation function of residuals
  
  
  #############################
  # potential confusion alert!!
  #
  # Here we use broken-stick piecewise linear approximations for two different purposes.  We use one broken-stick model
  # to approximate the nonlinear dependence of y on x, and another to approximate the dependence of the IRF on lag time
  # (thus allowing users to capture both the rapid changes in the IRF at short lag times, and the slow changes in the
  # IRF at long lag times, without needing to estimate IRF coefficients for huge numbers of individual lags.
  #
  # This is potentially confusing because both of these broken-stick approximations have knots.  To disambiguate these 
  # two types of knots, we use "xknots" to refer to the knots in the broken-stick representation of y's dependence on x,
  # and "knots" to refer to the knots in the broken-stick representation of the relationship between the IRF and
  # lag time.  The corresponding numbers of these knots are "nxk" and "nk".  
  
  
  
  
  # check range of x
  if (min(x, na.rm=TRUE)<0) stop("Fatal error in nonlinIRF: x's may not be negative")
  
  
  # if x is not a matrix, make it one
  x <- as.matrix(x)
  nx <- ncol(x)            # number of x variables
  
  
  
  # prepare xknots vector
  # if nonlinearity knots are specified as a vector, make it a matrix with the same number of columns as x
  xknots <- as.matrix(xknots)
  if (ncol(xknots)!=nx) {
    if (ncol(xknots)!=1) stop(paste0("Fatal error in nonlinIRF: ", ncol(xknots), " knot columns but ", nx, "x columns"))
    else xknots <- matrix(data=xknots, nrow=length(xknots), ncol=nx)   # if only a single column of xknots is supplied, then copy it for every column of x's
  }
  
  nxk <- NROW(xknots)+1                                           # number of knots above zero, for which nonlinear response will need to be estimated
  
  xprime <- matrix(data=0, nrow=NROW(x), ncol=nx*nxk)             # this is the x_prime matrix (matrix of increments of x between each pair of knots)
  
  for (i in 1:nx) sort(xknots[,i])                                # sort into ascending order
  
  kpts <- matrix(data=0, nrow=nrow(xknots)+2, ncol=ncol(xknots))  # temporary matrix that we will use to prepare knots
  for (i in 1:nx) {
    
    xx <- x[,i]
    kp <- xknots[,i]
    
    if (pct_xknots!=FALSE) {
      
      xx <- xx[!is.na(xx)]                                 # discard nonzero x's and get rid of na's
      xx <- sort(xx[xx>0])                                 # discard zeroes and sort remaining values in ascending order 
      
      if (pct_xknots==1) kp <- quantile(xx, probs=kp/100, na.rm=TRUE)   # convert percentile knots to values of x (note that these are quantiles of *nonzero* x's...)
      
      else if ((pct_xknots==2) | (pct_xknots==3)) { 
        # if pct_xknots==2, we set the knots according the percentile of the *running sum* of the x distribution.  That is, xknots=50 should find the value of X for which sum(x<X) equals sum(x>X).
        # if pct_xknots==3, we set the knots according the percentile of the *running sum* of the *squares* of the x distribution.  That is, xknots=50 should find the value of X for which sum((x<X)^2) equals sum((x>X)^2).
        
        xs <- cumsum(xx^(pct_xknots-1))                    # xs is the cumulative sum of the sorted values of xx (after they have been raised to the appropriate power)
        for (j in 1:length(kp)) kp[j] <- xx[which.min( abs(xs - max(xs)*kp[j]/100) )]    # this picks the values of xx that most closely delimit the corresponding fractions of the cumulative sum (or sum of squares if pct_xknots=3)
        
      } else stop("Fatal error in nonlinIRF: pct_xknots must be TRUE, FALSE, or an integer >=0 and <=3")
      
    }
    
     
    minx <- min(xx, na.rm=TRUE)
    maxx <- max(xx, na.rm=TRUE)
    
    if (sum(kp<=minx)>0) stop(paste0("Fatal error in nonlinIRF: one or more knots is <= minimum value of x for column ", i))
    if (sum(kp>=maxx)>0) stop(paste0("Fatal error in nonlinIRF: one or more knots is >= maximum value of x for column ", i))
    
    kpts[,i] <- c(0, kp, maxx)                # add knots of zero and max(x)
  } 
  xknots <- kpts # copy back to knots matrix
  
  if (sum(is.na(xknots))>0) stop("Fatal error in nonlinIRF: one or more knots is NA")
  
  # here we expand each column of x's into columns of xprimes (to detect nonlinear dependence on x)
  for (i in 1:ncol(x)) {
    
    xx <- x[,i] # take each column of x separately
    kp <- xknots[,i] # take each column of knots separately
    
    # now write the appropriate columns of the matrix of xprimes (increments of x between knots)
    
    for (el in 1:nxk) {                                                  # calling this "el" instead of l so that it is not confused with 1
      interval <- kp[el+1]-kp[el]                                        # interval between knots
      xprime[ ,(el+(i-1)*nxk)] <- pmax(0, pmin(xx-kp[el], interval))     # this should give the correct increment values within each interval between knots (equation 42 of K2022)
    }
  }
  
  
  
  ##############################################################
  # now call IRF with xprime, passing along all other parameters
  # since we are supplying xprimes here, the IRF routine will return beta-primes instead of betas
  
  zz <- IRF(y = y ,                                            
            x = xprime ,
            wt = wt ,
            m = m ,
            nk = nk ,
            nu = nu ,
            FD = FD ,
            h = h ,
            ARprob = ARprob ,
            ARlim = ARlim ,
            max.AR = max.AR ,
            complete = complete ,
            verbose = verbose ,
            robust = robust ,
            max.chunk = max.chunk )
  
  
  # now we need to convert betaprime to beta and system response y_k(x)
  betaprime <- zz$IRF                                              # zz$IRF is beta-prime, because we supplied xprime rather than x to IRF
  seprime <- zz$se                                                 # seprime is the se of beta-prime
  
  if (nk>0) nbeta <- nk
  else nbeta <- m+1
  
  beta <- matrix(data=0, nrow=nbeta, ncol=nx*nxk)                   # these will be the effective betas: accumulated averages over each interval
  beta_se <- matrix(data=0, nrow=nbeta, ncol=nx*nxk)                # and their standard errors
  ykx <- matrix(data=0, nrow=nbeta, ncol=nx*nxk)                    # these will be the total effect of x on y (defining the nonlinear response)
  ykx_se <- matrix(data=0, nrow=nbeta, ncol=nx*nxk)                 # and their standard errors
  
  for (i in 1:nx) {   # here's where we do the conversion itself...                                        
    
    ykx[ ,(1+(i-1)*nxk)] <- betaprime[ ,(1+(i-1)*nxk)]*(xknots[2,i]-xknots[1,i])
    ykx_se[ ,(1+(i-1)*nxk)] <- seprime[ ,(1+(i-1)*nxk)]*(xknots[2,i]-xknots[1,i])
    beta[ ,(1+(i-1)*nxk)] <- betaprime[ ,(1+(i-1)*nxk)]
    beta_se[ ,(1+(i-1)*nxk)] <- seprime[ ,(1+(i-1)*nxk)]
    
    for (el in 2:nxk) {        # equations 43 and 46 of K2022
      ykx[ ,(el+(i-1)*nxk)] <- ykx[ ,((el-1)+(i-1)*nxk)] + betaprime[ ,(el+(i-1)*nxk)]*(xknots[el+1,i]-xknots[el,i])                   # integrate over the piecewise linear approximation to the nonlinear dependence of y on x
      ykx_se[ ,(el+(i-1)*nxk)] <- sqrt(ykx_se[ ,(el-1+(i-1)*nxk)]^2 + (seprime[ ,(el+(i-1)*nxk)]*(xknots[el+1, i]-xknots[el,i]))^2)    # Gaussian error propagation (we may want to replace this with the full covariance matrix) 
      beta[ ,(el+(i-1)*nxk)] <- ykx[ ,(el+(i-1)*nxk)] / xknots[el+1,i]
      beta_se[ ,(el+(i-1)*nxk)] <- ykx_se[ ,(el+(i-1)*nxk)] / xknots[el+1,i]
    }
    
  } #next i
  
  
  ######################################################################
  # now calculate averages, and weighted averages, of IRF's and y_k(x)'s
  # first declare the necessary arrays...
  avg_IRF <- matrix(data=0, nrow=nbeta, ncol=nx)
  avg_se <- matrix(data=0, nrow=nbeta, ncol=nx)
  alt_avg_IRF <- matrix(data=0, nrow=nbeta, ncol=nx)
  alt_avg_se <- matrix(data=0, nrow=nbeta, ncol=nx)
  wtd_avg_IRF <- matrix(data=0, nrow=nbeta, ncol=nx)
  wtd_avg_se <- matrix(data=0, nrow=nbeta, ncol=nx)
  avg_ykx <- matrix(data=0, nrow=nbeta, ncol=nx)
  avg_ykx_se <- matrix(data=0, nrow=nbeta, ncol=nx)
  wtd_avg_ykx <- matrix(data=0, nrow=nbeta, ncol=nx)
  wtd_avg_ykx_se <- matrix(data=0, nrow=nbeta, ncol=nx)
  sumwt <- rep(0, nx)
  countx <- rep(0, nx)
  countx_nz <- rep(0, nx)
  avgxi <- rep(0, nxk)
  countxi <- rep(0, nxk)
  
  
  # beta is indexed by [lag, nx*nxk]
  # xprime is indexed by [time, nx*nxk]
  # betaprime is indexed by [lag, nx*nxk]
  # IRF is indexed by [lag, nx*nxk]
  # avg_IRF is indexed by [lag, nx]
  # 
  # ykx is the crossproduct (across rows) of xprime[time, nxk] and betaprime[lag, nxk]
  
  
  
  for (i in 1:nx) { # step through columns of x
    xi <- x[,i]                                        # select column of x
    xp <- xprime[!is.na(xi) , (1+(i-1)*nxk):(i*nxk)]   # excerpt from xprime (the relevant columns for this x), removing rows where xi is missing
    xi <- xi[!is.na(xi)]                               # remove missing values of xi
    b <- beta[ , (1+(i-1)*nxk):(i*nxk)]                # excerpt from beta (the relevant columns for this x)
    bs <- beta_se[ , (1+(i-1)*nxk):(i*nxk)]            # excerpt from beta_se (the relevant columns for this x)
    bs2 <- bs*bs                                       # squared standard errors of beta coefficients
    kp <- betaprime[ , (1+(i-1)*nxk):(i*nxk)]          # excerpt from betaprime (the relevant columns for this x)
    sp <- seprime[ , (1+(i-1)*nxk):(i*nxk)]            # excerpt from seprime (the relevant columns for this x)
    sp2 <- sp*sp                                       # square seprime (elementwise) because this will be necessary for error propagation
    xp2 <- xp*xp                                       # square xprime (for error propagation)
    
    # the error propagation here will look a bit odd because within each time step, the uncertainties in beta-prime are independent of one another (and thus these are added in quadrature)
    # but the uncertainties are not independent from one time step to the next (because errors in beta-prime affect all time steps equally) and thus are not averaged in quadrature
    
    xpkp <- xp %*% t(kp)                                    # this is ykx (rows are time steps, columns are lags)
    xp2sp2 <- xp2 %*% t(sp2)                                # error propagation for ykx
    xpsp <- sqrt(xp2sp2)                                    # standard error of ykx for each lag and time step
    avg_ykx[,i] <- colMeans(xpkp)                           # average the ykx's across all time steps (with non-missing x's)
    avg_ykx_se[,i] <- colMeans(xpsp)                        # error propagation for the mean of ykx (not averaged in quadrature -- see note above)
    wtd_avg_ykx[,i] <- colWeightedMeans(xpkp, w=xi)         # weighted sum of ykx's
    wtd_avg_ykx_se[,i] <- colWeightedMeans(xpsp, w=xi)      # error propagation for the weighted mean of ykx (not averaged in quadrature -- see note above)
    
    wtd_avg_IRF[,i] <- avg_ykx[,i]/mean(xi)                 # since ykx is IRF*x, the mean of ykx divided by the mean of x is the weighted mean of the IRF
    wtd_avg_se[,i] <- avg_ykx_se[,i]/mean(xi)                # error propagation for the weighted average of IRF (not averaged in quadrature -- see note above)
    
    
    # here we calculate average IRFs from linear interpolation of beta between pairs of knots, assuming beta=0 at x=0
    # this will be more reliable than the alternative method below if IRF approaches zero as x approaches zero
    ukp <- xknots[-1, i]       # vector of lower knots
    lkp <- xknots[-(nxk+1), i]  # vector of upper knots
    nxi <- length(xi)            # count of xi
    wj <- rep(0, nxk)             # weights for each of the betas in estimates of average beta
    
    for (j in 1:nxk) { # step through knots and calculate counts and averages of xi values between knots
      avgxi[j] <- mean(xi[(xi>lkp[j])&(xi<=ukp[j])])  # average xi between knots
      countxi[j] <- sum(((xi>lkp[j])&(xi<=ukp[j])))   # sum of booleans is count of xi between knots
    } 
    
    for (j in 1:nxk) {                  # now assign weights to each of the beta values, according to linear interpolation between knots
      if (j==nxk) wj[j] <- (avgxi[j]-lkp[j])/(ukp[j]-lkp[j])*countxi[j]/nxi    # weight for the last beta value
      else wj[j] <- (avgxi[j]-lkp[j])/(ukp[j]-lkp[j])*countxi[j]/nxi + (ukp[j+1]-avgxi[j+1])/(ukp[j+1]-lkp[j+1])*countxi[j+1]/nxi   # weights for all other beta values
    }
    
    avg_IRF[,i] <- rowSums(sweep(b, MARGIN=2, wj, FUN="*"))  # avg_IRF is weighted average
    avg_se[,i] <- sqrt(rowSums(sweep(bs2, MARGIN=2, wj*wj, FUN="*"))/sum(wj*wj))  # error propagation by weighted average of (squared) standard errors
    
    
    
    # here we calculate average IRFs by alternative method (will be more reliable if IRF does not approach zero as x approaches zero)
    xpkp_nz <- xpkp[(xi>0), ]                 # exclude any rows with x=0, since IRF at x=0 is undefined                           
    xpsp_nz <- xpsp[(xi>0), ]                 # exclude any rows with x=0, since IRF at x=0 is undefined
    xi_nz <- xi[(xi>0)]                       # exclude any rows with x=0, since IRF at x=0 is undefined
    alt_avg_IRF[,i] <- colMeans(sweep(xpkp_nz, MARGIN=1, as.array(xi_nz), FUN='/'))  #divide each row of ykx by x (nonzero values only) to get IRFs
    alt_avg_se[,i] <- colMeans(sweep(xpsp_nz, MARGIN=1, as.array(xi_nz), FUN='/'))   #error propagation (not averaged in quadrature -- see note above)
    
  }  # next column of x
  
  
  xknots = xknots[-1, , drop=FALSE] # now discard first row (zeroes) from knots
  
  
  #########################
  # now create column names
  
  # if vector of xnames does not exist, create it
  if (is.null(xnames) | length(xnames)!=ncol(x)) xnames <- paste0("x", 1:nx)
  
  ####################################
  # create vector of xnames with knots
  xknames <- rep("", nx*nxk)  # create vector for column names with x's and knots
  for (i in 1:nx) for (k in 1:nxk) {
    xstr <- xnames[i]
    nch <- nchar(xstr)
    ii <- 0  # need to search for location of first "|" by brute force, because | is a special character so we can't use regexpr!
    repeat{
      ii <- ii+1
      if (substr(xstr, ii, ii)=="|") {break}
      if (ii==nch) {break}
    }
    if ((ii==nch) || (ii==1)) xknames[k+(i-1)*nxk] <- paste0(xstr, "=", sprintf("%g", xknots[k,i]))   # if there is no "|" divider or it comes at the end or beginning
    else xknames[k+(i-1)*nxk] <- paste0(substr(xstr, 1, ii-1), "=", sprintf("%g", xknots[k,i]), substr(xstr, ii, nch))
  }
  
  colnames(beta) <- paste0("IRF_", xknames)
  colnames(beta_se) <- paste0("se_", xknames)
  colnames(ykx) <- paste0("ykx_", xknames)
  colnames(ykx_se) <- paste0("ykx_se_", xknames)
  colnames(avg_IRF) <- paste0("avg_IRF_", xnames)
  colnames(avg_se) <- paste0("avg_se_", xnames)
  colnames(avg_ykx) <- paste0("avg_ykx_", xnames)
  colnames(avg_ykx_se) <- paste0("avg_ykx_se_", xnames)
  colnames(wtd_avg_IRF) <- paste0("wtd_avg_IRF_", xnames)
  colnames(wtd_avg_se) <- paste0("wtd_avg_se_", xnames)
  colnames(wtd_avg_ykx) <- paste0("wtd_avg_ykx_", xnames)
  colnames(wtd_avg_ykx_se) <- paste0("wtd_avg_ykx_se_", xnames)
  colnames(alt_avg_IRF) <- paste0("alt_avg_IRF", xnames)
  colnames(alt_avg_se) <- paste0("alt_avg_se", xnames)
  
  
  
  
  
  
  # return results
  return(
    list(
      lags = zz$lags ,    # lags (in number of time steps)
      nxk = rep(nxk, nx) , # number of nonlinearity knots (vector of values for each x)
      xknots = xknots ,     # nonlinearity knots, with first row (zeroes) removed
      IRF = beta ,           # nonlinear impulse response function (beta) evaluated at each knot (except zero): matrix of nx*nxk columns and m+1 rows, corresponding to lags 0 through m)
      se = beta_se ,          # standard errors of IRF coefficients (matrix of nx*nxk columns and m+1 rows, corresponding to lags 0 through m)
      ykx = ykx ,              # nonlinear impulse response, expressed as contribution to y from x, evaluated at each knot (except zero): matrix of nx*nxk columns and m+1 rows, corresponding to lags 0 through m)
      ykx_se = ykx_se ,         # standard errors of ykx (matrix of nx*nxk columns and m+1 rows, corresponding to lags 0 through m)
      avg_IRF = avg_IRF ,        # arithmetic average of IRFs for nonzero x's (when x=0, IRF is undefined)
      avg_se = avg_se ,           # standard error of average of IRFs
      avg_ykx = avg_ykx ,          # time-averaged ykx
      avg_ykx_se = avg_ykx_se ,     # standard error of time-averaged ykx
      wtd_avg_IRF = wtd_avg_IRF ,    # weighted average of IRF
      wtd_avg_se = wtd_avg_se ,       # standard error of weighted average of IRF
      wtd_avg_ykx = wtd_avg_ykx ,      # averaged ykx weighted by input
      wtd_avg_ykx_se = wtd_avg_ykx_se , # standard error of weighted average ykx
      alt_avg_IRF = alt_avg_IRF ,        # average IRF by alternate method (better for IRFs that do not converge to zero at x=0)
      alt_avg_se = alt_avg_se ,           # standard error of average IRF by alternate method
      Kbb = zz$Kbb ,            # stack of covariance matrices for each x and each knot (except zero) (array of nx*nxk m+1 by m+1 matrices)
      n = zz$n ,                 # length of original x and y series
      n.eff = zz$n.eff ,          # effective sample size, accounting for uneven weighting
      n.nz = zz$n.nz ,             # number of nonzero values with nonzero weight at each lag of each x variable
      h = zz$h ,                    # order of AR correction that was applied
      phi = zz$phi ,                 # fitted AR coefficients (vector of length h)
      resid = zz$e ,                  # residuals
      resid.acf = zz$resid.acf ,       # autocorrelation function of residuals
      resid.pacf = zz$resid.pacf ,      # partial autocorrelation function of residuals
      ycomp = zz$ycomp                   # data table comparing measured and fitted y time series
      
      
    )
  )
  
  
}  # end nonlinIRF

#///////////////////
# END of nonlinIRF
#///////////////////








