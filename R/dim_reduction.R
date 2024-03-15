## fitting 3 gammma functions, sequentially, to kernels
## saving parameters and goodness of fit for analysis later

library(fitdistrplus)
library(here)
library(stringr)

## load catchment nums
to_test <- str_replace_all(list.files(path = here(".", "data", "raw_data"),                                                         # nolint
                            full.names = FALSE), "\\D", "")
to_test <- to_test[to_test != ""]
n <- nrow(to_test)

for (i in to_test) {
    # load kernel
    load(here("FD_kernels", paste(i, ".Rda", sep = "")))
    kernel <- as.numeric(kernel)

    # normalize kernel so that I can fit a gamma distr to it
    sum <- as.numeric(sum(kernel))
    n_kernel <- kernel * 100

    # fit gamma_1 (ideally should be furthest left)
    gamma_1 <- fitdist(n_kernel, distr = "gamma", method = "mle")

    plot(gamma_1)

    # get residuals and fit gamma_2 to them
    gamma_2 <- fitdist(kernel - gamma_1, distr = "gamma", method = "mle")

    plot(gamma_1 + gamma_2)
}
