# loading the required packages
library(ggplot2)
library(ggmap)
citation("ggmap")
library(here)
library(foreach)
library(dplyr)

## ----global tuning params-----------------------------------------------------
## ----NOTE: you also need to HARDCODE region choice below (LL54 & 56)----------
region_codes <- c("n.a.", "SA", "EU", "AU")
selected_criteria <- "best_aic"
selected_metric <- "kge"
metric_threshold <- 0
reload_basemaps <- TRUE

## ----load catchment results/atlas---------------------------------------------
gamma_three <- read.csv(here("data", "results", "clean", "3w_gauss25.csv"))
gamma_three$longitude <- as.numeric(gamma_three$longitude)
gamma_three$latitude <- as.numeric(gamma_three$latitude)
gamma_three$kge <- as.numeric(gamma_three[, "kge"])
gamma_three[, "nse"] <- as.numeric(gamma_three[, "nse"])
gamma_three[, "nrmse"] <- as.numeric(gamma_three[, "nrmse"])
gamma_three[, "r2"] <- as.numeric(gamma_three[, "r2"])
gamma_three[, "X"] <- NULL
gamma_three[is.na(gamma_three$region), "region"] <- "n.a."

gauss_three <- read.csv(here("data", "results", "clean", "3w_gamma25.csv"))
gauss_three$longitude <- as.numeric(gauss_three$longitude)
gauss_three$latitude <- as.numeric(gauss_three$latitude)
gauss_three$kge <- as.numeric(gauss_three[, "kge"])
gauss_three[, "nse"] <- as.numeric(gauss_three[, "nse"])
gauss_three[, "nrmse"] <- as.numeric(gauss_three[, "nrmse"])
gauss_three[, "r2"] <- as.numeric(gauss_three[, "r2"])
gauss_three[, "X"] <- NULL
gauss_three[is.na(gauss_three$region), "region"] <- "n.a." as.da

gamma_two <- read.csv(here("data", "results", "clean", "2w_gamma25.csv"))
gamma_two$longitude <- as.numeric(gamma_two$longitude)
gamma_two$latitude <- as.numeric(gamma_two$latitude)
gamma_two[, "median_mod_kge"] <- as.numeric(gamma_two[, "median_mod_kge"])
gamma_two[, "median_mod_nse"] <- as.numeric(gamma_two[, "median_mod_nse"])
gamma_two[, "median_mod_nrmse"] <- as.numeric(gamma_two[, "median_mod_nrmse"])
gamma_two[, "median_mod_r2"] <- as.numeric(gamma_two[, "median_mod_r2"])
gamma_two[, "X"] <- NULL
gamma_two[is.na(gamma_two$region), "region"] <- "n.a."

gamma_one <- read.csv(here("data", "results", "clean", "1w_gamma25.csv"))
gamma_one$longitude <- as.numeric(gamma_one$longitude)
gamma_one$latitude <- as.numeric(gamma_one$latitude)
gamma_one[, "kge"] <- as.numeric(gamma_one[, "kge"])
gamma_one[, "nse"] <- as.numeric(gamma_one[, "nse"])
gamma_one[, "nrmse"] <- as.numeric(gamma_one[, "nrmse"])
gamma_one[, "r2"] <- as.numeric(gamma_one[, "r2"])
gamma_one[, "X"] <- NULL
gamma_one[is.na(gamma_one$region), "region"] <- "n.a."

load(here("data", "catchment_atlas.Rda"))

## ----log in to google API and get regional basemaps, IF NECESSARY-------------
if (reload_basemaps) {
register_google(key = "AIzaSyAjm095BT1OytCgE6DJCvNkuJrr_pxIRMQ")
NA_basemap <- get_map(location = c(lon = mean(catchment_atlas$longitude[catchment_atlas$region == "NA"], na.rm = TRUE),          # nolint
                                lat = mean(catchment_atlas$latitude[catchment_atlas$region == "NA"], na.rm = TRUE)),             # nolint
                   maptype = "satellite",
                   zoom = 3)
SA_basemap <- get_map(location = c(lon = mean(catchment_atlas$longitude[catchment_atlas$region == "SA"], na.rm = TRUE),          # nolint
                                lat = mean(catchment_atlas$latitude[catchment_atlas$region == "SA"], na.rm = TRUE)),             # nolint
                   maptype = "satellite",
                   zoom = 3)
AU_basemap <- get_map(location = c(lon = mean(catchment_atlas$longitude[catchment_atlas$region == "AU"], na.rm = TRUE),          # nolint
                                lat = mean(catchment_atlas$latitude[catchment_atlas$region == "AU"], na.rm = TRUE)),             # nolint
                   maptype = "satellite",
                   zoom = 4)
EU_basemap <- get_map(location = c(lon = mean(catchment_atlas$longitude[catchment_atlas$region == "EU"], na.rm = TRUE),          # nolint
                                lat = mean(catchment_atlas$latitude[catchment_atlas$region == "EU"], na.rm = TRUE)),             # nolint
                   maptype = "satellite",
                   zoom = 5)

basemaps <- data.frame("NA" = NA_basemap,
                       "SA" = SA_basemap,
                       "AU" = AU_basemap,
                       "EU" = EU_basemap)
}

## ----select test set to work on-----------------------------------------------
clean_results <- gauss_three # gamma_two, gamma_one, gauss_three

## ----select and plot catchments, color-coded according to kge-----------------
data_subset <- clean_results[((clean_results$criteria == "best_bic")
                    & (clean_results[, "median_mod_kge"] >= 0)), ]                                                   # nolint


## ----SHOW MAP (MAP PARAMS ARE HARD-CODED. MAKE SURE THEY MATCH DATA PARAMS!!!)
ggmap(basemaps$"NA") +
  geom_point(data = clean_results[clean_results$median_mod_kge > 0,],
             aes(x = longitude, y = latitude, colour = median_mod_kge),
             size = 2) +
  scale_color_gradient(low = "#b10d0d", high = "#21e825") +
  guides(fill = FALSE, alpha = FALSE, size = FALSE)

boxplot <- boxplot(data_subset$kge)
