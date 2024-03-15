library(SlidingWindowReg)                                                                                                                       # nolint
library(readr)
library(lubridate) # package to handle date formats
library(testthat)
library(stringr)
library(dplyr)
library(here)

verbose <- FALSE

## ----load catchment data in easy-to-iterate form------------------------------
catchment_data <- list.files(path = "./my_work/data/WaterInputOutputData",
                            full.names = TRUE)

## ----load "catchment atlas"---------------------------------------------------
catchment_atlas <- read.csv("./my_work/data/WaterInputOutputData/all_gauges_metadata.csv")                                                        # nolint

## FOR TEST RUNS ONLY: select subset of catchments to run on
catchment_data <- catchment_data[300:301]
expect_equal(length(catchment_data), 2)

## ----init results table for this batch of catchments--------------------------
n_catchments <- length(catchment_data)
results <- data.frame(matrix(nrow = n_catchments, ncol = 10))
colnames(results) <- c("catchment_no",                                                                                                          # nolint
                    "catchment_latitude", "catchment_longitude",
                    "n_windows", "window_param_cv_sd",
                    "best_criteria",
                    "nrmse", "r2", "nse", "kge")
results$kge <- rep(-Inf, n_catchments)

## for loop iterating through each catchment------------------------------------
system.time(
for (i in 1:n_catchments) {
    ## ----get catchment number, latitude, & longitude and save them------------
    catchment_no <- as.integer(str_replace_all(catchment_data[i], c("./my_work/data/WaterInputOutputData/" = "", ".csv" = "")))                 # nolint
    results[i, "catchment_no"] <- catchment_no
    results[i, "catchment_latitude"] <- catchment_atlas$latitude[catchment_no]
    results[i, "catchment_longitude"] <- catchment_atlas$longitude[catchment_no]

    ## ----inform user (if verbose)---------------------------------------------
    if (verbose) { print(paste("processing catchment number ", catchment_no, " (", i, "/", n_catchments, ")", sep = "")) }                      # nolint

    ## ----load and format/clean raw data---------------------------------------
    data <- read_csv(catchment_data[i])
    composite_dates <- with(data, ymd(paste(data$Year, data$Month, data$Day)))
    data$date <- composite_dates

    ## ----for each training criteria-------------------------------------------
    for (criteria in as.list(c("best_rmse", "best_bic", "best_aic"))) { 
        ## ----inform user (if verbose)-----------------------------------------
        if (verbose) { print(paste("testing model using", criteria, "as training criteria.")) }                                                 # nolint
        ## ----run temporally-exclusive cross-validation------------------------
        n_folds <- 25
        cv_results <- data.frame(n_windows = rep(0.0, n_folds),
                                nrmse = rep(NA, n_folds),
                                r2 = rep(NA, n_folds),
                                nse = rep(NA, n_folds),
                                kge = rep(NA, n_folds))

        for (fold in 1:n_folds) {
            ## ----inform user (if verbose)-------------------------------------
            if (verbose) { print(paste("running cross-validation round ", fold, "...", sep = "")) }                                             # nolint

            ## ----get training and validation sets-----------------------------
            # get random index in first ~20% of sample indexes (edge exclusive
            # to prevent slice "funkiness")
            start_ind <- sample(2:((0.4 * nrow(data)) - 1), 1)
            end_ind <- start_ind + (0.6 * nrow(data))
            # use 60% of data, starting from random index, as training set
            training_set <- data[start_ind:end_ind, ]
            # other 20%, regardless of discontinuity, is validation set
            # mess is to handle weird cases
            validation_set <- rbind(data[1:(start_ind - 1),], data[(end_ind + 1):nrow(data),])                                                  # nolint
            # double check I did my slices right
            expect_equal((nrow(validation_set) + nrow(training_set)), nrow(data))                                                               # nolint

            ## ----train model--------------------------------------------------
            mod <- trainSWR(training_set$"Forcing_Precipitation_mm/d",
                                        training_set$"Observed_Streamflow_mm/d",
                                        iter = 3,
                                        param_selection = criteria)

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
                                                    nrow(mod$param))                                                                            # nolint
        }

        ## ----summarize cv results and update results depending on KGE---------
        final_results <- c(median(cv_results$n_windows),
                            sd(cv_results$n_windows),
                            cv_results[["nrmse"]],
                            cv_results[["r2"]],
                            cv_results[["nse"]],
                            cv_results[["kge"]])

        # guard against cases where results contain NaN/NA, which cannot be
        # compared; assuming this indicates poor fit, even just for one metric
        usable_results <- !anyNA(cv_results) & !sum(is.nan(as.matrix(cv_results)))                                                              # nolint
        if (usable_results) {
            if ((is.na(results[i, ]$kge)) ||
            (final_results[6] > results[i, ]$kge)) {
                results[i, 4:10] <- final_results
                results["best_criteria"] <- criteria
            }
        }
    }
}
)

save(results, file = "./my_work/data/results/results.Rda")
print(paste("Done.", n_catchments, "catchments tested and documented."))
