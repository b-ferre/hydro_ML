import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import scipy as scipy

## helper for scipy.optimize.curve_fit
def gamma_pdf(x, shape, loc, rate):
    scale = 1/rate
    return scipy.stats.gamma.pdf(x, shape, loc, scale)

def weight_fun(n):
    ## this function returns a weight vector of length n
    ## TODO: make these weights depend on len(x) for a more consistent weighting scheme for different length kernels
    a = 20          # controls where the "center" of the flipped sigmoid is
    b = 10          # controls difference between heighest and lowest weighted example
    k = 0.15        # controls "steepness" of transition from high weight to low weight
    x = np.arange(n)
    return -(b * (1 / (1 + np.exp(-k * (x - a))))) + b + 0.5

## front-focused, highly outlier-sensitive loss function for optimization
def custom_loss(params, x, y):
    assert(len(x) == len(y))
    # params are parameters for a gamma_pdf in order shape, loc, rate
    fit = gamma_pdf(x, params[0], params[1], params[2])
    # apply custom weighting scheme (see above)
    residuals = np.abs(y - fit)
    weighted_residuals = weight_fun(len(x)) * (residuals**2)

    ## part of loss comes from weighted residual calculation, other part comes from a term designed to emphasize hitting the highest point on the gamma dist.
    peak_y_i = np.argmax(y)
    peak_fit_i = np.argmax(fit)
    peak_loss = np.abs(fit[peak_y_i] - y[peak_y_i]) #+ np.abs(fit[peak_fit_i] - y[peak_fit_i])
    pre_peak_loss = np.sum(np.abs(fit[:peak_y_i] - y[:peak_y_i]))
    post_peak_range = 5
    post_peak_loss = np.sum(np.abs(fit[peak_y_i:peak_y_i + post_peak_range] - y[peak_y_i:peak_y_i + post_peak_range]))  ## encourages fit to immediately return to "flat" if the kernel does
    loss = np.sum(weighted_residuals) + (100 * peak_loss) + (70 * pre_peak_loss) #+ (50 * post_peak_loss / post_peak_range)
    print(loss)
    return loss

## load catchment nums
to_test = [s.replace(".csv", "") for s in os.listdir("./data/raw_data")]    # remove .csv from catchment data
to_test = to_test[:-3]                                                     # remove non-catchment data
n = len(to_test)

## get list of bad catchments (do not have kernel)
bad_catchments = list(pd.read_csv("./FD_kernels/bad_catchments.csv").iloc[0])

## for each catchment
for i in to_test :
    if not int(i) in bad_catchments:
        ## for debug
        print("working on catchment " + i + "...")

        ## get kernel as a vector
        kernel = pd.read_csv("./FD_kernels/" + i + ".csv").iloc[: , 1].values
        lags = pd.read_csv("./FD_kernels/" + i + ".csv").iloc[: , 0].values

        ## fit a gamma function to data using custom weighting scheme
        initial_params = [1.05, 0, 0.5]
        bounds = [(1, 5), (-20, 50), (0.1, 10)]
        best_params = scipy.optimize.minimize(custom_loss, initial_params, args = (lags, kernel), bounds = bounds)["x"]
        print(best_params)
        fit_1 = gamma_pdf(lags, best_params[0], best_params[1], best_params[2])

        ## remove for main run, just doing some visualization  right now
        plt.scatter(lags, kernel, s = 2, c = "blue")
        plt.scatter(lags, fit_1, s = 1, c = "red")
        plt.show()