library(here)
source(here("research", "kirschner_ar_model", "irfnnhs.r"))
source(here("research", "R", "helpers.R"))
library(readr)
library(hydroGOF)
library(data.table)
library(stats)
library(stringr)

############################################################################################################################################# nolint
##                                                        HELPER TEST FUNCTION                                                             ## nolint
############################################################################################################################################# nolint


## This function  performs k-fold cross validation to test Jim Kirschner's model
test_fdirf <- function(catchment_no, max_lag = 100, dir = "fwd", k = 25) {
    ## load catchment data
    all_data <- get_data(catchment_no)
    n <- nrow(all_data)

    cv_scores <- data.frame(fold = seq_len(k), kge = rep(0, k))
    for (fold in seq_len(k)) {
        ## get a random 60% continuous subsample of time series
        ## TODO: make subsample % scale with k folds
        subsample_size <- 0.4
        max <- floor(n * (1 - subsample_size))
        i_0 <- floor(runif(n = 1, min = 1, max = max))
        i_f <- i_0 + floor(n * subsample_size)
        data <- all_data[i_0:i_f, ]

        ## get training/testing split according to dir flag
        split <- floor(0.75 * nrow(data))
        if (dir == "fwd") {
            train <- data[1:split, ]
            test <- data[split:nrow(data), ]
        }
        if (dir == "bkwd") {
            n <- nrow(data)
            train <- data[(n - split):n, ]
            test <- data[1:(n -split), ]
        }

        tryCatch({
            ## train model, get kernels
            mod <- IRF(y = train$streamflow, x = train$precip,
                    FD = TRUE, m = max_lag)
            irf_ker <- rev(sapply(mod$IRF, as.numeric))
            ar_ker <- rev(sapply(mod$phi, as.numeric))},
        error = function(e) {
            print(paste("error (see below) while training model on catchment", catchment_no))                                                   ## nolint
            print(e)
            cv_scores[cv_scores$fold == fold, "kge"] <- NA
        })


        ## IMPORTANT NOTE:
        ## KGE'S obtained by this method are NOT representative of the ability of this model to predict streamflow                               ## nolint
        ## from rainfall alone - they are ONLY an approximate measure of the correctness of the kernel (see below).                             ## nolint
        #########################################################################################################                               ## nolint
        ## I use h days of measured streamflow to substitute for missing prediction data during the initial                                     ## nolint
        ## stages of prediction. Not doing so would cause a cascading effect of error in prediction - error                                     ## nolint
        ## that would NOT be the result of incorrectness in the kernel. As the purpose of this testing is to                                    ## nolint
        ## test the correctness of the kernel and the kernel only, this "peaking" at test data is justified, but                                ## nolint
        ## users/viewers of this project should be very aware of this fact. I do not want to misconstrue my results.                            ## nolint


        tryCatch({
        ## iteratively convolve coeffs with rainfall to predict streamflow
        ## TODO: make this not iterative... lol
        rain <- sapply(test$precip, as.numeric)
        t <- length(rain)
        pred <- rep(NA, t)
        m <- length(mod$IRF)
        for (i in (seq_len(t - m) + m)) {
            if (i == m + 1) {
                reg <- sum((irf_ker %*% rain[(i - m) : (i - 1)]), na.rm = TRUE)
                ar <- sum((ar_ker %*% test$streamflow[(i - mod$h) : (i - 1)]), na.rm = TRUE)                                                    ## nolint
                pred[i] <- reg + ar
            }
            else if (i - m <= mod$h) {
                h_acc <- i - m - 1
                h_fill <- mod$h - h_acc
                reg <- sum((irf_ker %*% rain[(i - m) : (i - 1)]), na.rm = TRUE)
                padded_streamflow <- append(test$streamflow[(i - mod$h) : (i - h_acc - 1)],                                                     ## nolint
                                            pred[(i - h_acc): (i - 1)])
                ar <- sum(ar_ker %*% padded_streamflow)
                pred[i] <- reg + ar
            } else {
                reg <- sum((irf_ker %*% rain[(i - m) : (i - 1)]), na.rm = TRUE)
                ar <- sum((ar_ker %*% pred[(i - mod$h) : (i - 1)]), na.rm = TRUE)                                                               ## nolint
                pred[i] <- reg + ar
            }
        }
        },
        error = function(e) {
            print(paste("error (see below) while testing model on catchment", catchment_no))                                                    ## nolint
            print(e)
            cv_scores[cv_scores$fold == fold, "kge"] <- NA
        })

        tryCatch({
            pred <- pred[!(is.na(test$streamflow))]
            streamflow <- test$streamflow[!(is.na(test$streamflow))]
            kge <- KGE(pred, streamflow)
            # print(paste("fold :", fold, "; kge :", kge))
            cv_scores[cv_scores$fold == fold, "kge"] <- kge
        },
        error = function(e) {
            print(paste("error (see below) when calculating kge in catchment", catchment_no))                                                   ## nolint
            cv_scores[cv_scores$fold == fold, "kge"] <- NA
        })
    }

    median_kge <- median(cv_scores$kge, na.rm = TRUE)
    mean_kge <- mean(cv_scores$kge, na.rm = TRUE)
    return(data.frame(median_kge = median_kge, mean_kge = mean_kge))
}


############################################################################################################################################    ## nolint
##                                                  RUN TEST ON ALL CATCHMENTS                                                                  ## nolint
############################################################################################################################################    ## nolint


## ----load catchment data in easy-to-iterate form------------------------------
to_test <- str_replace_all(list.files(path = here(".", "research", "data", "new_data",                                                          ## nolint
                            "raw_data", "by_catchment"), full.names = FALSE), "\\D", "")                                                        ## nolint
to_test <- to_test[to_test != ""]
n <- length(to_test)

test_scores <- data.frame(catchment_no = to_test,
                        median_kge = rep(NA, n),
                        mean_kge = rep(NA, n))

## TODO: remove code block above
load("./research/data/new_data/results/arma_model/max_FD_lags.Rda")
load("./research/data/new_data/results/arma_model/validation_testing/bad.Rda")

to_test <- to_test[!(to_test %in% bad)]
n <- length(to_test)
t0 <- as.numeric(Sys.time())
for (j in to_test[1:n]) {
    test <- test_fdirf(j, max_lag = max_practical_lags$max_practical_lag[max_practical_lags$catchment_no == j] + 20)                            ## nolint
    test_scores[test_scores$catchment_no == j, "median_kge"] <- test$median_kge
    test_scores[test_scores$catchment_no == j, "mean_kge"] <- test$mean_kge
    if ((as.numeric(match(j, to_test)) %% 10.0) == 0) {
        print(paste(round(match(j, to_test) / n, digits = 4) * 100, "% done;   ",                                                               ## nolint
        "estimated time remaining: ",
        round((((as.numeric(Sys.time()) - t0) / (match(j, to_test[1:n]))) * (n - match(j, to_test[1:n]))  / 60)),                               ## nolint
        " minutes...", sep = ""))
    }
}

View(test_scores)