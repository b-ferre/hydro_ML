## load packages
library(here)
library(data.table)             # data.table used inside model training
library(readr)                  # read_csv
library(hydroGOF)               # KGE
library(stringr)                # str_replace_all

## load Jim's model code
source(here("kirschner_ar_stuff", "irfnnhs.r"))

## load catchment data
to_test <- str_replace_all(list.files(path = here(".", "data", "raw_data"),                                                         # nolint
                            full.names = FALSE), "\\D", "")
to_test <- to_test[to_test != ""]
n <- nrow(to_test)

## load max lag dictionary
load("max_FD_lags.Rda")
res <- max_practical_lags

## turn off plotting if requested
plot <- FALSE

## set up list of "bad" catchments
bad_catchments <- list()

## for each catchment...
for (i in to_test) {
    ## get max practical lag from dictionary
    max_lag <- res[res$catchment_no == i, "max_practical_lag"]

    ## if max_lag = 5 this means that NO tested lag allowed for no negative
    ## kernel weights, so append this to the list of bad catchments and skip it
    if (max_lag <= 5) {
        bad_catchments <- c(bad_catchments, i)
        next
    }

    ## get data and make model using max lag (same params as during max lag run)
    mod <- IRF_i(i, FD = TRUE, lag = max_lag)

    ## if there are any negative kernel weights something is wrong, as model is
    ## deterministic and the goal of the discovery process was to avoid this, so
    ## throw a warning in the terminal for later use
    if (length(mod$IRF[mod$IRF < 0]) > 0) {
        print("ERROR: negative kernel weights... c_no:", i)
    }

    ## plot and save kernel with padded kernel weights so that x-axis of all
    ## kernels is comprable (100 max possible, given nature of discovery test)
    if (plot) {
        pdf(here("FD_kernel_graphs", paste(i, ".pdf", sep = "")))
        padded_coeffs <- mod$IRF
        if (length(mod$IRF) < 100) {
            padded_coeffs <- c(mod$IRF, rep(0, 100 - length(mod$IRF))) }
        plot(padded_coeffs)
        dev.off()
    }

    ## save IRF coeffs for plotting later (each in its own Rda)
    write.csv(mod$IRF, file = here("FD_kernels", paste(i, ".csv", sep = "")))
}
