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


## This function is a naive way to test Jim Kirschner's FD model
test_fdirf <- function(catchment_no, dir = "fwd") {
    ## load catchment data
    data <- get_data(catchment_no)
    ## find max practical lag
    max_lag <- find_max_practical_lag(catchment_no)
    if (max_lag == 5) {
        return(c(NA, 5))
    }

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
        return(c(NA, max_lag))
    })

    tryCatch({
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
    },
    error = function(e) {
        print(paste("error (see below) while testing model on catchment", catchment_no))                                                   ## nolint
        print(e)
        return(c(NA, max_lag))
    })

    tryCatch({ 
        kge <- KGE(pred, test$streamflow)
        return(c(kge, max_lag))},
    error = function(e) {
        print(paste("error (see below) when calculating kge in catchment", catchment_no))                                                   ## nolint
        print(e)
        return(c(NA, max_lag))
    }
    )
}


############################################################################################################################################## nolint
##                                                  RUN TEST ON ALL CATCHMENTS                                                              ## nolint
############################################################################################################################################## nolint


## ----load catchment data in easy-to-iterate form------------------------------
to_test <- str_replace_all(list.files(path = here(".", "research", "data", "new_data",                                                      ## nolint
                            "raw_data", "by_catchment"), full.names = FALSE), "\\D", "")
to_test <- to_test[to_test != ""]
n <- length(to_test)

# max_practical_lags <- data.frame(catchment_no = to_test,
#                                 max_practical_lag = rep(5, n))

# test_scores <- data.frame(catchment_no = to_test,
#                         kge = rep(NA, n))

t0 <- as.numeric(Sys.time())
for (i in to_test[1753:n]) {
    test <- test_fdirf(i)[1]
    kge <- test[1]
    max_lag <- test[2]
    test_scores[test_scores$catchment_no == i, "kge"] <- kge
    max_practical_lags[max_practical_lags$catchment_no == i, "max_practical_lag"] <- max_lag                                                ## nolint
    if ((as.numeric(match(i, to_test)) %% 10.0) == 0) {
        print(paste(round(match(i, to_test) / n, digits = 4) * 100, "% done;   ",                                                           ## nolint
        "estimated time remaining: ",
        round((((as.numeric(Sys.time()) - t0) / (match(i, to_test) - 1753)) * (n - match(i, to_test))  / 60)),                              ## nolint
        " minutes...", sep = ""))
    }
}

View(test_scores)

print("DEAL WITH SAVING MAX LAGS")