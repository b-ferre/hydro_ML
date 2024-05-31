require(leaflet)
require(leafpop)
require(ggplot2)
require(sf)
require(here)

map <- function() {
  ## load data
  load("./data/scores.Rda")
  load("./data/bad.Rda")
  load("./data/atlas.Rda")

  ## get list of all catchment numbers I have a (valid) kge data for
  good_ids <- scores$id

  ## get a list of all catchment numbers which I tested and failed on (in some capacity)                          # nolint 
  bad_ids <- bad

  ## pre-process data for mapping
  # good_locs given score and kernel
  locs <- data.frame(index = atlas$gridcode,
                      lng = atlas$longitude,
                      lat = atlas$latitude)

  ## remove all data where lng, lat are invalid 
  ## TODO: talk to dr ameli about this and try to recover these locations somehow                                 # nolint
  locs <- locs[((locs$lat < 87) & (locs$lat > -87)) &
               ((locs$lng > -181) & (locs$lng < 181)), ]

  # only get locations for catchments which I tested, seperated out depending on how testing went                 # nolint
  bad_locs <- locs[locs$index %in% bad_ids, ]
  good_locs <- locs[locs$index %in% good_ids, ]

  ## set up kge color-scaling
  ## input x's must be \leq 1
  map_to_colors <- function(x) {
    # NA values are black
    if (x < 0 || is.na(x) || x == -Inf) {
      return("000000")
    }
    # values < 0 are red; values in [0, 1] scaled linearly from red to green;
    # values > 1 are green;
    green <- max(0, min(1, x / 2))
    return(rgb((1 - green), green, 0))                                                                            # nolint
  }
  colors <- vector(mode = "character", length = nrow(good_locs))
  scaled_scores <- data.frame(id = scores$id,
                            score = (scores$median_kge + 1.2))                                                    # nolint


  ## build kernel graphs for good catchments
  print("plotting kernels...")
  kernels <- list()
  for (j in seq_len(length(good_locs$index))) {
      i <- good_locs$index[j]
      data <- data.frame(sapply(read.csv(paste("./data/kernels/", i, ".csv", sep = "")), as.numeric))             # nolint
      colnames(data) <- c("lag", "IRF_coeff")
      plt <- ggplot(data = data, aes(x = lag, y = IRF_coeff)) +
                  geom_point() +
                  geom_line(linetype = "dashed", color = "#2d14d3") +
                  xlim(0, max(30, nrow(data))) +
                  labs(title = paste("gridcode", i),
                      subtitle = paste("median_kge :", round(scores[scores$id == i, "median_kge"], digits = 4)))  # nolint
      kernels[[j]] <- plt
      colors[j] <- map_to_colors(scaled_scores[scaled_scores$id == i, "score"])                                   # nolint
  }

  ## build map
  t0 <- Sys.time()
  print("building map (this may take ~5 minutes, just wait for the > to appear in the console)...")               # nolint
  #suppressMessages(
    map <- leaflet()    %>%
    addTiles()      %>%
    addCircleMarkers(data = bad_locs,
                  radius = 1,
                  group = "bad",
                  color = "#2f2d2d",
                  opacity = 0.5,
                  popup = sapply(bad_locs$index, as.character)) %>%
    addCircleMarkers(data = good_locs,
                  radius = 3,
                  group = "good",
                  color = colors,
                  opacity = 1) %>%
    addPopupGraphs(kernels, group = "good", width = 300, height = 300)                                            # nolint
  #)
  tf <- Sys.time()
  print(paste("elapsed time:", tf - t0))


  return(map)
}