library(here)
library(stringr)

## get a list of all catchment data in CSV form
to_process <- list.files("./Sync/Benjamin/NewData_2024/Daily Data")                                                 ## nolint

n <- length(to_process)
pb <- txtProgressBar(min = 0, max = n, style = 3, width = 50, char = "=")

snow_dominated <- c()

for (i in 1:n) {
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

    ## get rainfall and streamflow and format them so that there is one year per column                             ## nolint
    years <- unique(format(as.Date(raw_data$date), "%Y"))
    rain <- matrix(nrow = length(years), ncol = 365)
    streamflow <- matrix(nrow = length(years), ncol = 365)
    for (j in seq_len(length(years))) {
        syear <- paste0(years[j])
        ## TODO: clean/double check lines below
        year_filter <- str_detect(raw_data$date, syear)
        rain[j, ] <- raw_data$precipitation_mmd[year_filter]
        streamflow[j, ] <- raw_data$streamflow_mmd[year_filter]                                                     ## nolint
    }

    ## remove seasonal signals such that the mean for each day of the year is 0
    for (j in seq_len(365)) {
        rain[, j] <- rain[, j] - mean(rain[, j], na.rm = TRUE)
        streamflow[, j] <- streamflow[, j] - mean(streamflow[, j], na.rm = TRUE)
    }

    data_dir_path <- "./research/hflm-ds_model/slapdash_test/data/"
    saveRDS(rain, file = paste0(data_dir_path, "grid", gridcode, "_rain.rds"))                                      ## nolint
    saveRDS(streamflow, file = paste0(data_dir_path, "grid", gridcode, "_streamflow.rds"))                          ## nolint

    setTxtProgressBar(pb, i)
}

save(snow_dominated, file = "./snow_dom.rda")