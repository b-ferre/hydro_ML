library(here)
source(here("kirschner_ar_model", "irfnnhs.r"))
source(here("R", "helpers.R"))
library(readr)
library(hydroGOF)
library(data.table)
library(stats)

## This function is a naive way to test Jim Kirschner's FD model
naive_FD_IRF_test <- function(catchment_no) {
    ## load catchment data
    data <- get_data(catchment_no)
    max_lag <- get_max_practical_lag(catchment_no)

    ## get training/testing split (preferences forward prediction)
    split <- floor(0.75 * nrow(data))
    train <- data[1:split, ]
    test <- data[split:nrow(data), ]

    mod <- IRF(y = train$streamflow, x = train$precip,
            FD = TRUE, m = max_lag)
    irf_ker <- rev(sapply(mod$IRF, as.numeric))
    ar_ker <- rev(sapply(mod$phi, as.numeric))

    ## iteratively convolve coeffs with rainfall to predict streamflow
    ## TODO: make this not iterative... lol
    rain <- sapply(test$precip, as.numeric)
    pred <- rep(0, length(rain) - length(irf_ker))
    for (i in seq_len(length(pred))) {
        if (i == 1) {
            print(irf_ker)
            print(rain[i : (i + length(irf_ker) - 1)])
            print(length(irf_ker))
            print(length(rain[i : (i + length(irf_ker) - 1)]))
            pred[i] <- irf_ker %*% rain[i :(i + length(irf_ker) - 1)]
            print(pred[i])
        }
        else if (i <= mod$h + 1) {
            reg <- (irf_ker %*% rain[i :(i + length(irf_ker) - 1)])
            print(reg)
            ar <- (ar_ker[1:(mod$h - i)] %*% pred[1:(mod$h - 1)])
            print(ar)
            continue
        }
    }

    kge <- KGE(pred, test$streamflow[length(irf_ker):length(test$streamflow)])
    return(kge)
}

test_scores <- data.frame(catchment_no = to_test,
                        kge = rep(-Inf, n))

for (i in to_test) {
    if (!(i %in% bad)) {
        test_scores[test_scores$catchment_no == i, "kge"] <- naive_FD_IRF_test(i)                                                   # nolint

    }
    if ((as.numeric(match(i, to_test)) %% 10.0) == 0) {
        print(paste(round(match(i, to_test) / n, digits = 4) * 100, "% done...", sep = ""))                                          # nolint
    }
}


View(test_scores)