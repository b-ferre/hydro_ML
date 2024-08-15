library(here)
library(stringr)

## get a list of all catchment data in CSV form
to_process <- list.files("./Sync/Benjamin/NewData_2024/Daily Data")                                                 ## nolint

## get seasonality classifications for all catchments
season_classifications <- read.csv("./seasonality_testing/data/months_probability_aggregated.csv")                      ## nolint

n <- length(to_process)
pb <- txtProgressBar(min = 0, max = n, style = 3, width = 50, char = "=")

snow_dominated <- c()

for (i in 2926:n) { ## FIXME: replace with 1:n once done debug
    print(i)
    file <- to_process[i]
    raw_data <- as.data.frame(read.csv(paste0("./Sync/Benjamin/NewData_2024/Daily Data/", file)))                   ## nolint
    gridcode <- as.integer(str_replace_all(file, ".csv", ""))

    ## remove leap days if any are present
    raw_data <- raw_data[!str_detect(raw_data$date, "02-29"), ]

    ## determine whether catchment is rain dominated; if not, skip it
    temp_filter <- (raw_data$temperature_C > 0)
    snow_sum <- sum(raw_data$precipitation_mmd[!temp_filter])
    precip_sum <- sum(raw_data$precipitation_mmd[temp_filter])
    snow_frac <- snow_sum / precip_sum
    if (snow_frac >= 0.1) {
        setTxtProgressBar(pb, i)
        snow_dominated <- c(snow_dominated, gridcode)
        next
    }

    ## set precip to 0 for any days with negative temperature as to exclude snowfall
    raw_data[!temp_filter, "precipitation_mmd"] <- 0

    ## FIXME: make sure the below approach to pruning time series is not too naive
    ## (setting in/output to 0 for off season)
        ## note: seasons[i + 1] represents the season class for the ith month  
    seasons <- season_classifications[season_classifications$gridcode == gridcode, ]
    if (any(is.na(seasons))) {
        print(paste("seasonality data corrupted for catchment", gridcode))
        next
    }

    num_to_mon <- function(num) {
        months <- c("Jan", "Feb", "Mar", "Apr", "May", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
        return(months[num])
    }

    ## first isolate dormant season data
    mon_ids <- num_to_mon(as.integer(substr(raw_data$date, 6, 7)))

    mon_to_grow_bool <- function(month) {
        return(seasons[month] >= 0)
    }

    grow_index <- mon_to_grow_bool(mon_ids)
    dorm_index <- !grow_index
    grow_raw <- data.frame(raw_data)
    dorm_raw <- data.frame(raw_data)
    grow_raw[dorm_index, "precipitation_mmd"] <- 0
    grow_raw[dorm_index, "streamflow_mmd"] <- 0
    dorm_raw[grow_index, "precipitation_mmd"] <- 0
    dorm_raw[grow_index, "streamflow_mmd" <- 0]

    for (season in c("grow", "dorm")) {

        data <- grow_raw
        if (season == "dorm") {
            data <- dorm_raw
        }

        ## get rainfall and streamflow and format them so that there is one year per row                             ## nolint
        years <- unique(format(as.Date(data$date), "%Y"))
        rain <- matrix(nrow = length(years), ncol = 365)
        streamflow <- matrix(nrow = length(years), ncol = 365)
        for (j in seq_len(length(years))) {
            syear <- paste0(years[j])
            year_filter <- str_detect(data$date, syear)
            rain[j, ] <- data$precipitation_mmd[year_filter]
            streamflow[j, ] <- data$streamflow_mmd[year_filter]                                                     ## nolint
        }

        rownames(rain) <- years
        rownames(streamflow) <- years

        ## remove seasonal signals such that the mean for each day of the year is 0
        for (j in seq_len(365)) {
            rain[, j] <- rain[, j] - mean(rain[, j], na.rm = TRUE)
            streamflow[, j] <- streamflow[, j] - mean(streamflow[, j], na.rm = TRUE)
        }

        ## delete any year that is >30% NA for either x or y
        ## (fixes issues in UK where last 4 years are completely missing)
        learnable <- function(row) {
            return(!((sum(is.na(row)) / length(row)) > 0.3))
        }
        learnable_rain <- apply(X = rain, FUN = learnable, MARGIN = 1)
        learnable_sflow <- apply(X = streamflow, FUN = learnable, MARGIN = 1)
        learnable_years <- learnable_rain & learnable_sflow
        rain <- rain[learnable_years, ]
        streamflow <- streamflow[learnable_years, ]

        data_dir_path <- paste0("./seasonality_testing/data/pp", season, "_data/")

        saveRDS(rain, file = paste0(data_dir_path, "grid", gridcode, "_rain.rds"))                                      ## nolint
        saveRDS(streamflow, file = paste0(data_dir_path, "grid", gridcode, "_streamflow.rds"))                          ## nolint
    }

    setTxtProgressBar(pb, i)
}

print()
print()
print("dont forget to save snow_dominated")
