map <- function() {
## load libraries (quietly)
suppressMessages(require(leaflet))
suppressMessages(require(leafpop))
suppressMessages(require(here))
suppressMessages(require(stringr))
suppressMessages(require(ggplot2))
suppressMessages(require(sf))

## load data
load("./data/scores.Rda")
load("./data/bad.Rda")
load("./data/atlas.Rda")

## get list of all catchment numbers I have (good) data for
ids <- scores$id[scores$kge != -Inf]

## pre-process data for mapping
locs <- data.frame(index = atlas$gridcode,
                    lng = atlas$longitude,
                    lat = atlas$latitude)
locs <- locs[locs$index %in% ids, ]
bad_locs <- locs[locs$index %in% bad, ]
good_locs <- locs[!(locs$index %in% bad), ]

## set up kge color-scaling
## input x's must be \leq 1
map_to_colors <- function(x) {
  # bad values/values less than 0 are assigned "red"
  if (x < 0 || is.na(x) || x == -Inf) {
    return("ff0000")
  }
  # Values between 0 and 1 are scaled linearly from red to green
  green <- max(0, min(1, x / 2))
  return(rgb((0.5 - green) * 7 / 5, green * 11 / 10, 0))                         # nolint
}
colors <- vector(mode = "character", length = nrow(good_locs))
scaled_scores <- data.frame(id = scores$id,
                          score = (scores$kge + (1 - max(scores$kge, na.rm = TRUE))))                       # nolint
## build kernel graphs for good catchments
print("plotting kernels...")
kernels <- list()
for (j in seq_len(length(good_locs$index))) {
    i <- good_locs$index[j]
    data <- data.frame(sapply(read.csv(paste("./data/kernels/", i, ".csv", sep = "")), as.numeric))         # nolint
    colnames(data) <- c("lag", "IRF_coeff")
    plt <- ggplot(data = data, aes(x = lag, y = IRF_coeff)) +
                geom_point() +
                geom_line(linetype = "dashed", color = "#2200ff") +
                xlim(0, max(30, nrow(data))) +
                labs(title = paste("gridcode", i),
                    subtitle = paste("kge :", round(scores[scores$id == i, "kge"], digits = 4)))            # nolint
    kernels[[j]] <- plt
    colors[j] <- map_to_colors(scaled_scores[scaled_scores$id == i, "score"])                               # nolint
}

## build map
t0 <- Sys.time()
print("building map (this may take ~5 minutes, just wait for the > to appear in the console)...")           # nolint
suppressMessages(
  map <- leaflet()    %>%
  addTiles()      %>%
  addCircleMarkers(data = bad_locs,
                radius = 1,
                group = "bad",
                color = "#2f2d2d",
                opacity = 1,
                popup = sapply(bad_locs$index, as.character)) %>%
  addCircleMarkers(data = good_locs,
                radius = 1.2,
                group = "good",
                color = colors,
                opacity = 1) %>%
  addPopupGraphs(kernels, group = "good", width = 300, height = 300)                                       # nolint
)
tf <- Sys.time()
print(paste("elapsed time:", tf - t0))


return(map)
}
