
## TODO: implement dependency tags (here and read.csv)
get_data <- function(catchment_no, old_data = FALSE) {
    if (old_data) {
    print("warning: old data is being retrieved and presumably used - is this intentional?")                                            # nolint
    raw <- read.csv(here("data", "raw_data",
                    paste(catchment_no, ".csv", sep = "")))
    data <- data.frame(
        date = as.Date(raw$Month, raw$Day, raw$Year, format = "%m.%d.%Y"),
        precip = raw$"Forcing_Precipitation_mm.d",
        snowmelt = raw$"HBV_Snowmelt_mm.d",
        streamflow = raw$"Observed_Streamflow_mm.d"
    )
    return(data)
    } else {
        raw <- read.csv(here("research", "data", "new_data", "raw_data",
                        "by_catchment", paste(catchment_no, ".csv", sep = ""))) 
        data <- data.frame(
        date = raw$date,
        precip = raw$precipitation_mmd,
        streamflow = raw$streamflow_mmd)
    }
}

## the goal of this function is to (somewhat heuristically) test the maximum
## practical lag for a given catchment by iterating through possible lags
## I do this by running the standard IRF model with iteratively larger lags
## (and complete = TRUE)
find_max_practical_lag <- function(catchment_no,
                                   max_considered = 100,
                                   first_diff = TRUE,
                                   old_data = FALSE) {
    # minimum lag to test is set at 5 to avoid numerical errors later on
    max_practical_lag <- 5

    # get data
    data <- get_data(catchment_no, old_data = old_data)

    # iterate through lags to check, training a
    for (lag in seq(10, max_considered, step = 10)) {
        mod <- IRF(data$streamflow, data$precip, FD = first_diff, m = lag, complete = FALSE)                                             # nolint

        if (length(mod$IRF[mod$IRF < 0]) == 0) {
            max_practical_lag <- lag
        } else if (lag - max_practical_lag > 20) {
            # early stopping condition to fight optimization bias
            break()
        }
    }

    for (lag in seq(max_practical_lag, max_practical_lag + 9)){
        mod <- IRF(data$streamflow, data$precip, FD = first_diff, m = lag, complete = FALSE)                                             # nolint

        if (length(mod$IRF[mod$IRF < 0]) == 0) {
            max_practical_lag <- lag
        } else if (lag - max_practical_lag > 20) {
            # early stopping condition to fight optimization bias
            break()
        }
    }

    # print(paste("finished testing catchment ",
    #             catchment_no,
    #             "... max practical lag: ",
    #             max_practical_lag, sep = ""))
    return(max_practical_lag)
}

get_max_practical_lag <- function(catchment_no, old_data = FALSE) {
    if (old_data) {
        load("./research/data/old_data/results/arma_model/max_FD_lags.Rda")
        return(max_practical_lags[max_practical_lags$catchment_no == catchment_no, "max_practical_lag"])                                # nolint
    } else {
        load("./data/new_data/results/arma_model/max_FD_lags.Rda")
        return(max_practical_lags[max_practical_lags$catchment_no == catchment_no, "max_practical_lag"])                                # nolint
    }
}

## return an IRF model trained on the entire dataset for catchment number i
irf_i <- function(catchment_no, FD = FALSE, h = NULL, lag = 100, old_data = FALSE) {                                                    # nolint
    data <- get_data(catchment_no, old_data = old_data)
    mod <- IRF(data$streamflow, data$precip, FD = FD, h = NULL, m = lag, complete = FALSE)                                              # nolint

    return(mod)
}