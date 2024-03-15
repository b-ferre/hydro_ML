## load packages
library(here)
library(data.table)             # data.table used inside model training
library(readr)                  # read_csv
library(hydroGOF)               # KGE
library(stringr)                # str_replace_all

## load catchment data
to_test <- str_replace_all(list.files(path = here(".", "data", "raw_data"),                                                         # nolint
                            full.names = FALSE), "\\D", "")
to_test <- to_test[to_test != ""]
n <- nrow(to_test)

for (i in to_test) {
    data <- get_data(i)

    # generate regular ACF
    pdf(here("ACF_plots", "regular", paste(i, ".pdf", sep="")))
    acf(data$streamflow, type = "correlation", na.action = na.pass)
    dev.off()

    # first difference target time-series
    k <- length(data$streamflow)
    fd_series <- data$streamflow[2:k] - data$streamflow[1:(k-1)]

    # generate first-diffed ACF
    pdf(here("ACF_plots", "first-diffed", paste(i, ".pdf", sep="")))
    acf(fd_series, type = "correlation", na.action = na.pass)
    dev.off()
}
