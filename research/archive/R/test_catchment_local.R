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
catchment_data <- list.files(path = here(".", "data", "raw_data"),                                                                 # nolint
                            full.names = TRUE)

## ----SET TESTING PARAMS (cluster_size, training_iters, n_criteria)------------                                                                # nolint
criteria_to_test <- c("best_bic", "best_aic")
n_criteria_to_test <- length(criteria_to_test)
cluster_size <- 5
training_iters <- 3
cv_folds <- 5
paste("TESTING PARAMS")
print(paste("num cores in cluster:", cluster_size))
print(paste("training_iters:", training_iters))
print(paste("cv_folds:", cv_folds))

## ----set up parallelization---------------------------------------------------
cl <- makeCluster(cluster_size)
registerDoParallel(cl, cores = cluster_size)
print(paste("setting up parallelization...\n",
"number of workers registered:", getDoParWorkers()))

## ---- get i & map to unique (catch., crit.) index-----------------------------
i  <- 1
cr_i <- as.integer((i %% n_criteria_to_test) + 1)
c_i <- as.integer(ceiling(i / length(criteria_to_test)))
print(paste("test index:", i))

## ----get catchment no, lat, and lon-------------------------------------------
catchment_no <- as.integer(str_replace_all(catchment_data[c_i], "[a-zA-Z/_:.]", ""))                                                            # nolint
print(paste("catchment number:", catchment_no))
print(paste("training criteria to test:", criteria_to_test[cr_i]))

## ----load and format/clean raw data---------------------------------------
data <- read_csv(catchment_data[c_i])
composite_dates <- with(data, ymd(paste(data$Year, data$Month, data$Day)))
data$date <- composite_dates

## ----run (&time) parralelized, temporally-exclusive cross-validation------
start_time <- Sys.time()
cv_results <- foreach(fold = 1:cv_folds, .packages = "SlidingWindowReg", .combine = rbind) %dopar% {                                                           # nolint
    ## ----get training and validation sets-----------------------------
    # get random index in first ~30% of sample indexes (edge exclusive
    # to prevent slice "funkiness")
    start_ind <- sample(2:((0.3 * nrow(data)) - 1), 1)
    end_ind <- start_ind + (0.7 * nrow(data))
    # use 60% of data, starting from random index, as training set
    training_set <- data[start_ind:end_ind, ]
    # other 20%, regardless of discontinuity, is validation set
    # mess is to handle weird cases
    validation_set <- rbind(data[1:(start_ind - 1),], data[(end_ind + 1):nrow(data),])                                                          # nolint                                                               # nolint

    ## ----train model--------------------------------------------------
    mod <- trainSWR(training_set$"Forcing_Precipitation_mm/d",
                                training_set$"Observed_Streamflow_mm/d",
                                iter = training_iters,
                                param_selection = criteria_to_test[cr_i])

    ## ----get validation errors----------------------------------------
    pred <- predict(mod, validation_set$"Forcing_Precipitation_mm/d")
    metrics <- eval_all(pred, validation_set$"Observed_Streamflow_mm/d")

    ## ----return (all but rmse) errors for this iteration of CV--------
    c(nrow(mod$param), metrics$nrmse, metrics$r2, metrics$nse,
    metrics$kge)
}

stopCluster(cl)
print(paste("elapsed time:", Sys.time() - start_time))

colnames(cv_results) <- c("n_windows", "nrmse", "r2", "nse", "kge")
cv_results <- as.data.frame(cv_results)

catchment_results <- data.frame(
    catchment_no = catchment_no,
    criteria = criteria_to_test[cr_i],
    n_windows = median(cv_results$n_windows),
    nrmse =  mean(cv_results$nrmse),
    r2 = mean(cv_results$r2),
    nse = mean(cv_results$nse),
    kge = mean(cv_results$kge)
)

save(catchment_results, file = here(".", "data", "results",
paste(catchment_no, "_", criteria_to_test[cr_i], ".Rda", sep = "")))