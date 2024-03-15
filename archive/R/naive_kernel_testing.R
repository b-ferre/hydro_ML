source(here("kirschner_ar_stuff", "irfnnhs.r"))
source(here("R", "helpers.R"))
library(readr)
library(hydroGOF)
library(data.table)
library(stats)

## This function is a naive way to test Jim Kirschner's FD model
## FIXME: currently only testing the actual kernel, ignoring autocorellation
naive_FD_IRF_test <- function(catchment_no) {
## load catchment data
catchment_data <- read_csv(
    paste(here(".", "data", "raw_data"), "/", catchment_no, ".csv", sep = ""),
    show_col_types = FALSE)

n <- nrow(catchment_data)

load("max_FD_lags.Rda")
max_lag = max_practical_lags[max_practical_lags$catchment_no == i, "max_practical_lag"]

## get training/testing split (preferences forward prediction)
split <- floor(0.75 * n)
train <- catchment_data[1:split, ]
test <- catchment_data[split:n, ]

## FIXME: change both instances of HBV below to forcing_precip if needed
mod <- IRF(y = train$"Observed_Streamflow_mm/d", x = train$"Forcing_Precipitation_mm/d",
        FD = TRUE, m = max_lag)
ker <- sapply(rev(mod$IRF), as.numeric)

## convolve rainfall with kernel to rudimentarily predict streamflow
rain <- sapply(test$"Forcing_Precipitation_mm/d", as.numeric)
pred <- stats::filter(rain, ker)

kge <- KGE(pred, test$"Observed_Streamflow_mm/d")
return(kge)
}


## get (bad) catchments
to_test <- str_replace_all(list.files(path = here(".", "data", "raw_data"),                                                         # nolint
                            full.names = FALSE), "\\D", "")
to_test <- to_test[to_test != ""]
n <- length(to_test)
bad_catchments <- read.csv(here(".", "FD_kernels", "bad_catchments.csv"))
bad <- list()
for (i in seq_len(ncol(bad_catchments))) {
    bad[[i]] <- as.numeric(bad_catchments[1, i])
}

test_scores <- data.frame(catchment_no = to_test,
                        kge = rep(-Inf, n))

for (i in to_test) {
    if (!(i %in% bad)) {
        test_scores[test_scores$catchment_no == i, "kge"] <- naive_FD_IRF_test(i)                                                   # nolint

    }
    if ((as.numeric(match(i, to_test)) %% 10.0) == 0) {
        print(paste(round(match(i, to_test) / n, digits = 4) * 100, "% done...", sep = ""))                                               # nolint
    }
}


View(test_scores)