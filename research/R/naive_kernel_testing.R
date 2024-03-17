library(here)
source(here("kirschner_ar_model", "irfnnhs.r"))
source(here("R", "helpers.R"))
library(readr)
library(hydroGOF)
library(data.table)
library(stats)
library(stringr)


############################################################################################################################################# nolint
##                                                        HELPER TEST FUNCTION                                                             ## nolint
############################################################################################################################################# nolint


## This function is a naive way to test Jim Kirschner's FD model
test_fdirf <- function(catchment_no, dir = "fwd") {
    ## load catchment data
    data <- get_data(catchment_no)
    max_lag <- get_max_practical_lag(catchment_no)

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

    ## train model, get kernels
    mod <- IRF(y = train$streamflow, x = train$precip,
            FD = TRUE, m = max_lag)
    irf_ker <- rev(sapply(mod$IRF, as.numeric))
    ar_ker <- rev(sapply(mod$phi, as.numeric))

    ## iteratively convolve coeffs with rainfall to predict streamflow
    ## TODO: make this not iterative... lol
    rain <- sapply(test$precip, as.numeric)
    t <- length(rain)
    pred <- rep(NA, t)
    m <- length(mod$IRF)
    for (i in (seq_len(t - m) + m)) {
        if (i - m == 1) {
            pred[i] <- sum(irf_ker %*% rain[(i - m) : (i - 1)])
        } else if (i - m <= mod$h) {
            h_acc <- i - m - 1
            reg <- sum(irf_ker %*% rain[(i - m) : (i - 1)])
            ar <- sum(ar_ker[1:h_acc] %*% pred[(i - h_acc): (i - 1)])
            pred[i] <- reg + ar
        } else {
            reg <- sum(irf_ker %*% rain[(i - m) : (i - 1)])
            ar <- sum(ar_ker %*% pred[(i - mod$h) : (i - 1)])
            pred[i] <- reg + ar
        }
    }

    tryCatch({
        kge <- KGE(pred, test$streamflow)
        return(kge)}, error = function(e) {
            print(paste("fatal error when calculating kge in catchment", catchment_no))                                                     # nolint
            return(-Inf)
        }
    )
}


############################################################################################################################################# nolint
##                                                  RUN TEST ON ALL CATCHMENTS                                                             ## nolint
############################################################################################################################################# nolint


## ----load catchment data in easy-to-iterate form------------------------------
to_test <- str_replace_all(list.files(path = here(".", "data", "raw_data"),                                                                 # nolint
                            full.names = FALSE), "\\D", "")
to_test <- to_test[to_test != ""]
n <- length(to_test)

## load bad catchments
load("./zipped_map/data/bad.Rda")

test_scores <- data.frame(catchment_no = to_test,
                        kge = rep(-Inf, n))

for (i in to_test) {
    if (!(i %in% bad)) {
        test_scores[test_scores$catchment_no == i, "kge"] <- test_fdirf(i)                                                   # nolint

    }
    if ((as.numeric(match(i, to_test)) %% 10.0) == 0) {
        print(paste(round(match(i, to_test) / n, digits = 4) * 100, "% done...", sep = ""))                                          # nolint
    }
}

View(test_scores)