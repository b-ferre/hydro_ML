require(leaflet)
require(leafpop)
require(ggplot2)
require(sf)
require(here)
require(stringr)
require(jpeg)
require(htmltools)
require(base64enc)
require(leaflet)
require(leafpop)
source("./testing_funs.R")

M <- 364
max_lag <- 150

map <- function() {
  load("./data/atlas.Rda")
  load("./data/mappable.rda")
  load("./data/scores.rda")
  load("./data/snow_dom.rda")

  ## pre-process data for mapping
  locs <- data.frame(index = atlas$gridcode,
                      lng = atlas$longitude,
                      lat = atlas$latitude)

  ## pull out snow dominated catchment locations
  snow_locs <- locs[locs$index %in% snow_dominated, ]

  ## trim locs such that only mappable 
  good_locs <- locs[locs$index %in% mappable, ]

  ## create list of kernels which is in the same order as locs
  ## so that binding works correctly
  ## TODO: replace with code that generates ggplot objects because this makes final file smaller
  print("plotting kernels")                                                                                       # nolint
  pb <- txtProgressBar(min = 1, max = nrow(good_locs), style = 3, width = 150, char = "~")                        # nolint
  kernels <- list()
  for (i in seq_len(nrow(good_locs))) {
    grid <- good_locs$index[i]
    kern_i <- as.numeric(read.csv(paste0("./data/raw_estimates/beta_est_", grid, ".csv"))$pred_beta)
    plot_i <- plot_bvec_better(kern_i, M, max_lag)
    kernels[[i]] <- plot_i
    setTxtProgressBar(pb, i)
  }
  print("done!")
  print("\n")

  map_to_colors <- function(x) {
    # NA values are black
    if (x < 0 || is.na(x) || x == -Inf) {
      return(rgb(1, 0, 0))
    }
    # values < 0 are red; values in [0, max_score] scaled linearly from red to green;                             # nolint
    green <- max(0, min(1, x / max(scores[, metric])))
    return(rgb((1 - green), green, 0.2))                                                                          # nolint
  }

  ## build map
  print("building map (this may take ~5 minutes, just wait for the > to appear in the console)...")               # nolint
  pb <- txtProgressBar(min = 1, max = nrow(good_locs), style = 3, width = 150, char = "~")                        # nolint
    map <- leaflet()    %>%
    addTiles()          %>%
    addCircleMarkers(data = snow_locs,
                    radius = 2,
                    color = "grey",
                    opacity = 0.85)

    for (i in 1:nrow(good_locs)) {
      loc <- good_locs[i, ]
      score <- round(scores[scores$gridcode == loc$index, metric], digits = 4)                             # nolint
      suppressMessages({
      map <- addCircleMarkers(map, data = loc,
                  radius = 3,
                  group = paste(loc$index),
                  color = map_to_colors(score),
                  opacity = 1,
                  popup = paste("gridcode,  ", metric, "  :", loc$index, ",  ", score))                                      # nolint
      map <- addPopupGraphs(map, graph = list(kernels[[i]]), group = paste(loc$index), width = 400, height = 300)})           # nolint
      setTxtProgressBar(pb, i)
    }

  print("done!")


  return(map)
}

library(htmlwidgets)
metric <- "kge"
kgemap <- map()
saveWidget(kgemap, "hflm-ds_map_kge.html")


metric <- "r2"
r2map <- map()
saveWidget(r2map, "hflm-ds_map_r2.html")

metric <- "nse"
nsemap <- map()
saveWidget(nsemap, "hflm-ds_map_nse.html")
