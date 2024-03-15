source("./kirschner_ar_model/irfnnhs.r")
source("./R/helpers.R")
library(gridExtra)
library(ggplot2)
library(here)
library(dplyr)
library(data.table)

## only doing this for catchments that have decent performance;
## get a list of such catchments
load("./data/results/arma_model/validation_testing/ARMA_m=100_results.Rda")
to_do <- results$catchment_no[results$median_kge > 0]
instabilities <- data.frame(id = to_do,
                        translational_instability = rep(0, length(to_do),
                        sparsity_instability = rep(0, length(to_do))))

## FIXME: there are (seemingly random) catchments that break the code below

## annual stability
for (i in to_do) {
    print(i)
    data <- get_data(i)
    n <- nrow(data)
    splits <- seq(1, n, floor(n / 6))

    kernels <- matrix(ncol = get_max_practical_lag(i) + 1, nrow = length(splits) - 1)                                                   # nolint
    plots <- list()
    for (j in seq_len(length(splits) - 1)) {
        print(j)
        data_ <- data.frame(data[splits[j] : splits[j + 1], ])
        tryCatch({
            print("here")
            mod <- IRF(data_$streamflow, data_$precip, FD = TRUE, h = NULL,
                        m = get_max_practical_lag(i), complete = FALSE)
            print("here now")
            kernels[j, ] <- as.numeric(mod$IRF)},
        error = function(e) {
            print(e)
            kernels[j, ] <- rep(0, get_max_practical_lag(i) + 1) })
    }

    ## record kernel instability
    mean_kernel <- colMeans(kernels, na.rm = TRUE)
    instabilities$translational_instability[instabilities$id == i] <- sum(sweep(kernels, 2, mean_kernel, "-") ^ 2, na.rm = TRUE)        # nolint

    ## plotting
    pdf(paste("./data/results/arma_model/stabilitiy_testing/translational/graphs(n=6)", i, ".pdf", sep = ""))                           # nolint
    par(mfrow = c(2, ceiling((length(splits) - 1) / 4)))
    for (j in seq_len(length(splits) - 1)) {
        print(j)
        print(kernels[j, ])
        print("oh my god I'm gonna plot...")
        plot(lag = seq_len(length(kernels[j, ])), kernels[j, ], xlab = "lag", ylab = "IRF_coeff",                                       # nolint
            sub = paste("trained on", data_[1, "date"], "-", data_[nrow(data_), "date"]),                                               # nolint
            ylim = c(min(kernels[j, ]), min(100, max(kernels[j, ]))))
    }
    dev.off()
}
