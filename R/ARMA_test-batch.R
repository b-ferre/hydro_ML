library(here)
library(dplyr)
library(data.table)
library(nloptr)
library(stringr, include.only = "str_replace_all")
library(hydroGOF)
library(readr)

source(here("kirschner_ar_stuff", "irfnnhs.r"))

## ----load catchment data in easy-to-iterate form------------------------------
to_test <- str_replace_all(list.files(path = here(".", "data", "raw_data"),                                                         # nolint
                            full.names = FALSE), "\\D", "")
to_test <- to_test[to_test != ""]
n <- length(to_test)

## ----set testing parameters---------------------------------------------------
max_lag <- 100                # maximum lag to consider for regression
train_size <- 0.7             # determines the training/validation split
random_subsample_size <- 0.8  # ...
folds <- 3                    # ...

results <- data.frame(catchment_no = to_test,
                      n_bad_runs = rep(0, n),
                      mean_kge = rep(0, n),
                      median_kge = rep(0, n),
                      AR_order = rep(0, n),
                      my_kge = rep(NA, n))

## ----do validation testing AND graph full set kernel for each catchment-------
source(here("R", "ARMA_test_func.R"))
for (i in to_test) {
    ## plotting run
    plot_IRF(i, max_lag = max_lag)
    ## validation testing
    res <- robust_validate_IRF(i, max_lag = max_lag,
                        train_size = train_size,
                        random_subsample_size = random_subsample_size,
                        folds = folds)
    results[results$catchment_no == i, "mean_kge"] <- res$mean_kge
    results[results$catchment_no == i, "median_kge"] <- res$median_kge
    results[results$catchment_no == i, "AR_order"] <- res$median_AR_order
    results[results$catchment_no == i, "n_bad_runs"] <- res$problems
}

## ----append *my* model's median kge for each of the tested catchments
load("data/results/window_model/clean/2w_gamma25.Rda")
load("data/raw_data/catchment_atlas.Rda")

for (i in to_test) {
    kge_ <- clean_results[(clean_results$catchment_no == i) &
                            (clean_results$criteria == "best_bic"), "median_mod_kge"]       # nolint
    if (!(identical(kge_, numeric(0)))) {
        results[results$catchment_no == i, "my_kge"] <- kge_
    }
    country <-  catchment_atlas[catchment_atlas$gridcode == i, "country"]
    results[results$catchment_no == i, "country"] <- country

}

View(results, "m=100_ALL_results")