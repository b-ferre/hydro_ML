library(here)
library(dplyr)
library(data.table)
library(nloptr)
library(stringr, include.only = "str_replace_all")
library(hydroGOF)
library(readr)

source(here("kirschner_ar_stuff", "irfnnhs.r"))

to_test <- str_replace_all(list.files(path = here(".", "data", "raw_data"),                                                         # nolint
                            full.names = FALSE), "\\D", "")
to_test <- to_test[to_test != ""]

to_test <- to_test[2000:length(to_test)]

n <- length(to_test)

for (i in to_test) {
    m_p_lag <- find_max_practical_lag(i, first_diff = TRUE)
    res[res$catchment_no == i, "max_practical_lag"] <- m_p_lag

    print(paste("done with catchment", i, "     max practical lag =", m_p_lag))
}

max_practical_lags <- res

View(max_practical_lags, "Max Practical First-Diffed Lags")


## RUN THIS DUMMY
save(max_practical_lags, file = "max_FD_lags.Rda")