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
catchment_data <- list.files(path = here("..", "data", "raw_data"),                                                                             # nolint
                            full.names = FALSE)

## ----SET TESTING PARAMS (cluster_size, training_iters, n_criteria)------------                                                                # nolint
criteria_to_test <- c("best_bic", "best_aic", "best_rmse")
n_criteria_to_test <- length(criteria_to_test)
cluster_size <- 10
training_iters <- 3
cv_folds <- 10
paste("TESTING PARAMS")
print(paste("num cores in cluster:", cluster_size))
print(paste("training_iters:", training_iters))
print(paste("num criteria to test (tested sequentially):", n_criteria_to_test))
print(paste("cv_folds (run in parallel):", cv_folds))

## ----set up parallelization---------------------------------------------------
cl <- makeCluster(cluster_size)
registerDoParallel(cl, cores = cluster_size)
print("setting up parallelization...")
print(paste("number of workers registered:", getDoParWorkers()))

## ---- get catchment index-----------------------------------------------------
c_i  <- as.integer(commandArgs(trailingOnly = TRUE)[1])

## ----get catchment no, lat, and lon-------------------------------------------
catchment_no <- as.integer(str_replace_all(catchment_data[c_i], "[a-zA-Z/_:.]", ""))                                                            # nolint
print(paste("catchment number:", catchment_no))

## ----load and format/clean raw data-------------------------------------------
data <- read_csv(paste(here("..", "data", "raw_data"), "/", catchment_data[c_i], sep = ""))                                                          # nolint
composite_dates <- with(data, ymd(paste(data$Year, data$Month, data$Day)))
data$date <- composite_dates

for (criteria in criteria_to_test) {
print(paste("testing", criteria, "as training criteria..."))

## ----run (& time) parralelized, temporally-exclusive cross-validation---------
cv_start_time <- Sys.time()
cv_results <- foreach(fold = 1:cv_folds, .packages = "SlidingWindowReg", .combine = rbind) %dopar% {                                            # nolint
    ## ----get training and validation sets-------------------------------------
    # get random index in first ~30% of sample indexes (edge exclusive
    # to prevent slice "funkiness")
    start_ind <- sample(2:((0.3 * nrow(data)) - 1), 1)
    end_ind <- start_ind + (0.7 * nrow(data))
    # use 60% of data, starting from random index, as training set
    training_set <- data[start_ind:end_ind, ]
    # other 20%, regardless of discontinuity, is validation set
    # mess is to handle weird cases
    validation_set <- rbind(data[1:(start_ind - 1),], data[(end_ind + 1):nrow(data),])                                                          # nolint                                                               # nolint

    ## ----train model----------------------------------------------------------
    mod <- trainSWR(training_set$"Forcing_Precipitation_mm/d",
                                training_set$"Observed_Streamflow_mm/d",
                                iter = training_iters,
                                param_selection = criteria)

    ## ----get validation errors------------------------------------------------
    pred <- predict(mod, validation_set$"Forcing_Precipitation_mm/d")
    metrics <- eval_all(pred, validation_set$"Observed_Streamflow_mm/d")

    ## ----return (all but rmse) errors for this iteration of CV----------------
    c(nrow(mod$param), metrics$nrmse, metrics$r2, metrics$nse,
    metrics$kge)
}

print(paste(cv_folds, "-fold cross-validation completed. elapsed time: ", Sys.time() - cv_start_time, sep = ""))

colnames(cv_results) <- c("n_windows", "nrmse", "r2", "nse", "kge")
cv_results <- as.data.frame(cv_results)

catchment_results <- data.frame(
    catchment_no = catchment_no,
    criteria = criteria,
    n_windows = median(cv_results$n_windows),
    nrmse =  mean(cv_results$nrmse),
    r2 = mean(cv_results$r2),
    nse = mean(cv_results$nse),
    kge = mean(cv_results$kge)
)

save(catchment_results, file = here("..", "data", "results",
paste(catchment_no, "_", criteria, ".Rda", sep = "")))
}

stopCluster(cl)
print(paste("batch finished. total elapsed time:", Sys.time() - t_0))