library(here)
library(stringr)
library(data.table)

## load Jim's model code & helpers
source(here("research", "kirschner_ar_model", "irfnnhs.r"))
source("./research/R/helpers.R")

to_test <- str_replace_all(list.files(path = here(".", "research", "data", "new_data",                                                      ## nolint
                            "raw_data", "by_catchment"), full.names = FALSE), "\\D", "")                                                    ## nolint
n <- length(to_test)

load("./research/data/new_data/results/arma_model/max_FD_lags.Rda")

for (i in to_test[2578:n]) {
    data <- get_data(i)
    max_lag <- get_max_practical_lag(i)
    mod <- IRF(x = data$precip, y = data$streamflow, m = max_lag, FD = TRUE)
    write.csv(mod$IRF, file = here("research/FD_kernels", paste(i, ".csv", sep = "")))
}