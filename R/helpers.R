
## TODO: implement dependency tags (here and read.csv)
get_data <- function(catchment_no) {
    raw <- read.csv(here("data", "raw_data",
                    paste(catchment_no, ".csv", sep = "")))

    data <- data.frame(
        precip = raw$"Forcing_Precipitation_mm.d",
        snowmelt = raw$"HBV_Snowmelt_mm.d",
        streamflow = raw$"Observed_Streamflow_mm.d"
    )

    return(data)
}

## the goal of this function is to (somewhat heuristically) test the maximum
## practical lag for a given catchment by iterating through possible lags
## I do this by running the standard IRF model with iteratively larger lags
## (and complete = TRUE)
find_max_practical_lag <- function(catchment_no,
                                   max_considered = 100,
                                   first_diff = FALSE) {
    # minimum lag to test is set at 5 to avoid numerical errors later on
    max_practical_lag <- 5

    # get data
    data <- get_data(catchment_no)

    # iterate through lags to check, training a
    for (lag in 10:max_considered) {
        mod <- IRF(data$streamflow, data$precip, FD = first_diff, m = lag, complete = FALSE)                                             # nolint

        if (length(mod$IRF[mod$IRF < 0]) == 0) {
            max_practical_lag <- lag
        } else if (lag - max_practical_lag > 20) {
            # early stopping condition to fight optimization bias
            break()
        }
    }

    print(paste("finished testing catchment ",
                catchment_no,
                "... max practical lag: ",
                max_practical_lag, sep = ""))
    return(max_practical_lag)
}


get_learnability <- function(catchment_no) {
    ## perfect score (1.0) means that every single day in the time series has
    ## both target and input datapoints. otherwise I just calculate the percent
    ## of days that have sufficient data
    data <- get_data(catchment_no)

    input_density <- length(data$precip[!(is.na(data$precip) | is.nan(data$precip))]) / length(data$precip)                             # nolint
    target_density <- length(data$streamflow[!(is.na(data$streamflow) | is.nan(data$streamflow))]) / length(data$streamflow)            # nolint

    return((0.5 * input_density) + (0.5 * target_density))
}


## return an IRF model trained on the entire dataset for catchment number i
IRF_i <- function(catchment_no, FD = FALSE, h = NULL, lag = 100) {
    data <- get_data(catchment_no)
    mod <- IRF(data$streamflow, data$precip, FD = FD, h = NULL, m = lag, complete = FALSE)

    return(mod)
}