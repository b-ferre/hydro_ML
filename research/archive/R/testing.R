library(combinat)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(knitr)
library(methods)
library(nloptr)
library(parallel)
library(pbapply)
library(rdist)
library(Rdpack)
library(stats)
library(lifecycle)
library(here)
library(readr, include.only = "read_csv")
library(lubridate, include.only = "ymd")
library(testthat, include.only = "expect_equal")
library(here)
library(dplyr)
library(stringr, include.only = "str_replace_all")
library(foreach)
library(doParallel)

#install.packages(here("model_code", "SlidingWindowReg_0.1.1.tar.gz"))
#library(SlidingWindowReg)
install.packages(here(".", "model_code", "SWRGamma_0.1.0.tar.gz"))
library(SWRGamma)

c_i <- 3
criteria <- "best_bic"
training_iters <- 2

## ----load catchment data in easy-to-iterate form------------------------------
catchment_data <- list.files(path = here(".", "data", "raw_data"),                                                                             # nolint
                            full.names = FALSE)

## ----get catchment no, lat, and lon-------------------------------------------
catchment_no <- as.integer(str_replace_all(catchment_data[c_i], "[a-zA-Z/_:.]", ""))                                                           # nolint
print(paste("testing catchment number", catchment_no, "..."))

## ----load and format/clean raw data-------------------------------------------
data <- read_csv(paste(here(".", "data", "raw_data"), "/", catchment_data[c_i], sep = ""))                                                     # nolint
composite_dates <- with(data, ymd(paste(data$Year, data$Month, data$Day)))
data$date <- composite_dates

## ----get training and validation sets-------------------------------------
    # get random index in first ~30% of sample indexes (edge exclusive
    # to prevent slice "funkiness")
    start_ind <- sample(2:(floor(0.3 * nrow(data)) - 1), 1)
    end_ind <- start_ind + ceiling(0.7 * nrow(data))
    # use 60% of data, starting from random index, as training set
    training_set <- data[start_ind:end_ind, ]
    # other 20%, regardless of discontinuity, is validation set
    # mess is to handle weird cases
    validation_set <- rbind(data[1:(start_ind - 1),], data[(end_ind + 1):nrow(data),])                                                          # nolint

    ## ----train model----------------------------------------------------------
    mod <- SWRGamma::trainSWR(training_set$"Forcing_Precipitation_mm/d",
                                training_set$"Observed_Streamflow_mm/d",
                                iter = training_iters,
                                param_selection = criteria)
    ## ----get validation errors------------------------------------------------
    pred <- predict.SWR(mod, validation_set$"Forcing_Precipitation_mm/d")
    metrics <- eval_all(pred, validation_set$"Observed_Streamflow_mm/d")
    first_window <- which.max(mod$param[, "delta"])
    second_window <- which.min(mod$param[, "delta"])
    results <- c(nrow(mod$param), mod$param[first_window, 1], mod$param[first_window, 2], mod$param[first_window, 3],                           # nolint
      mod$param[second_window, 1], mod$param[second_window, 2], mod$param[second_window, 3],                                                    # nolint
      metrics$nrmse, metrics$r2, metrics$nse, metrics$kge)
    #colnames(results) <- c("n_windows", "delta1", "shape1", "rate1", "delta2", "shape2", "rate2", "nrmse", "r2", "nse", "kge")
    print(results)
