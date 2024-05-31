library(SlidingWindowReg)
# package to handle date formats
library(readr, include.only = "read_csv")

# package to handle date formats
library(lubridate, include.only = "ymd")

# useful packages for catching bugs
library(testthat, include.only = "expect_equal")
library(here)

library(dplyr)

# package for handling file/path names as strings
library(stringr, include.only = "str_replace_all")

# packages for parallelization
library(foreach)
library(doParallel)

## ----load catchment data in easy-to-iterate form------------------------------
catchment_data <- list.files(path = here("..", "data", "WaterInputOutputData"),                                                                 # nolint
                            full.names = TRUE)

## SET MINI-BATCH TESTING PARAMS (n_catchmnets, cluster_size, training_iters, n_criteria)                                                       # nolint
criteria_to_test <- c("best_bic", "best_aic")
n_criteria_to_test <- length(criteria_to_test)
cluster_size <- 20
training_iters <- 3
cv_folds <- 10
print(paste("num catchments to process:", length(catchment_data)))
print(paste("num cores in cluster:", cluster_size))
print(paste("no. of training criteria being tested:", n_criteria_to_test))
print(paste("training_iters:", training_iters))
print(paste("cv_folds:", cv_folds))

## ----set up parallelization---------------------------------------------------
cl <- makeCluster(cluster_size)
registerDoParallel(cl, cores = cluster_size)
print(paste("number of workers registered:", getDoParWorkers()))

start_time <- Sys.time()
##----run (& time) for-loop parallelized over each (catchment, criteria) pair---
full_results <- foreach(i = 1:(length(catchment_data) * n_criteria_to_test),
.packages = c("SlidingWindowReg", "readr", "lubridate", "dplyr", "stringr"), .combine = rbind, .errorhandling = "stop") %dopar% {               # nolint

## ----map unique (catch., crit.) index to catchment and criteria indexes-------
cr_i <- as.integer((i %% n_criteria_to_test) + 1)
c_i <- as.integer(ceiling(i / length(criteria_to_test)))

## ----get catchment no, lat, and lon-------------------------------------------
    catchment_no <- as.integer(str_replace_all(catchment_data[c_i], "[a-zA-Z/_:.]", ""))                                                        # nolint

    ## ----load and format/clean raw data---------------------------------------
    data <- read_csv(catchment_data[c_i])
    composite_dates <- with(data, ymd(paste(data$Year, data$Month, data$Day)))
    data$date <- composite_dates

    ## ----run temporally-exclusive cross-validation------------------------
    cv_results <- data.frame(n_windows = rep(0.0, cv_folds),
                            nrmse = rep(NA, cv_folds),
                            r2 = rep(NA, cv_folds),
                            nse = rep(NA, cv_folds),
                            kge = rep(NA, cv_folds))

    for (fold in 1:cv_folds) {
        ## ----get training and validation sets-----------------------------
        # get random index in first ~30% of sample indexes (edge exclusive
        # to prevent slice "funkiness")
        start_ind <- sample(2:((0.3 * nrow(data)) - 1), 1)
        end_ind <- start_ind + (0.7 * nrow(data))
        # use 60% of data, starting from random index, as training set
        training_set <- data[start_ind:end_ind, ]
        # other 20%, regardless of discontinuity, is validation set
        # mess is to handle weird cases
        validation_set <- rbind(data[1:(start_ind - 1),], data[(end_ind + 1):nrow(data),])                                                      # nolint                                                               # nolint

        ## ----train model--------------------------------------------------
        mod <- trainSWR(training_set$"Forcing_Precipitation_mm/d",
                                    training_set$"Observed_Streamflow_mm/d",
                                    iter = training_iters,
                                    param_selection = criteria_to_test[cr_i])

        ## ----get validation errors----------------------------------------
        pred <- predict(mod, validation_set$"Forcing_Precipitation_mm/d")
        metrics <- eval_all(pred, validation_set$"Observed_Streamflow_mm/d")

        ## ----record (all but rmse) errors for this iteration of CV--------
        cv_results[fold, c("nrmse",
                            "r2",
                            "nse",
                            "kge",
                            "n_windows")] <- c(metrics$nrmse,
                                                metrics$r2,
                                                metrics$nse,
                                                metrics$kge,
                                                nrow(mod$param))
    }

    ## return vectorized results for this (catchment, criteria) pairing-----
    # FORMAT: c_num, lat, lon, criteria, n_windows, window_sd, nrmse, r2, nse, kge                                                              # nolint
   c(catchment_no, criteria_to_test[cr_i],
        median(cv_results$n_windows), sd(cv_results$n_windows),
        mean(cv_results$nrmse), mean(cv_results$r2),
        mean(cv_results$nse), mean(cv_results$kge))
}

stopCluster(cl)

colnames(full_results) <- c("catchment_no", "criteria",
"n_windows", "window_sd", "nrmse", "r2", "nse", "kge")

print(paste("elapsed time:", Sys.time() - start_time))
save(full_results, file = here("..", "data", "results", "hail_mary_results.Rda"))