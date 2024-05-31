library(here)

catchment_atlas <- read.csv(here("data", "raw_data", "all_gauges_metadata.csv"))

result_files <- list.files(here("data", "results", "raw", "2w_gamma25_NEW_HEUR_INIT"), full.names = TRUE)
load(result_files[1])
clean_results <- catchment_results
for (i in 2:length(result_files) - 1) {
    load(result_files[i])
    clean_results <- rbind(clean_results, catchment_results)
}

## remove NaN's
clean_results <- na.omit(clean_results)

get_lon <- function(c_num) {
    matching_row <- catchment_atlas[catchment_atlas$gridcode == c_num, ]
    return(matching_row$longitude)
}

get_lat <- function(c_num) {
    matching_row <- catchment_atlas[catchment_atlas$gridcode == c_num, ]
    return(matching_row$latitude)
}

get_loc <- function(c_num) {
    matching_row <- catchment_atlas[catchment_atlas$gridcode == c_num, ]
    return(c(matching_row$longitude, matching_row$latitude))
}

clean_results$longitude <- as.numeric(lapply(clean_results$catchment_no, get_lon))              # nolint
clean_results$latitude <- as.numeric(lapply(clean_results$catchment_no, get_lat))               # nolint

get_region <- function(tuple) {
    lon <- tuple[1]
    lat <- tuple[2]
    if (lat > 10) {
        if (lon > -150 & lon < -20) {
            return("n.a.")
        } else {
            return("e.u.")
        }
    }
    else if (lon > -100 & lon < -10) {
        return("s.a.")
    }
    else {
        return("a.u.")
    }
}

locs <- lapply(clean_results$catchment_no, get_loc)
regions  <- lapply(locs, get_region)
clean_results$region <- regions

locs <- lapply(catchment_atlas$gridcode, get_loc)
regions  <- lapply(locs, get_region)
catchment_atlas$region <- regions
catchment_atlas <- as.data.frame(catchment_atlas)

save(clean_results, file = here("data", "results", "clean", "2w_gamma25_NEW_HEUR_INIT.Rda"))
#write.csv(as.matrix(clean_results), here("data", "results", "clean", "3w_gauss25.csv"))
