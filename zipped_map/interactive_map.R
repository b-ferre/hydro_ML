## library checks for running on new machines
if (system.file(package = "leaflet") == "") {
  install.packages("leaflet")
}
if (system.file(package = "leafpop") == "") {
  install.packages("leafpop")
}
if (system.file(package = "here") == "") {
  install.packages("here")
}
if (system.file(package = "stringr") == "") {
  install.packages("stringr")
}
if (system.file(package = "ggplot2") == "") {
  install.packages("ggplot2")
}
if (system.file(package = "sf") == "") {
  install.packages("sf")
}

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

## get list of all catchment numbers I have data for
ids <- scores$id

## pre-process data for mapping
locs <- data.frame(index = atlas$gridcode,
                    lng = atlas$longitude,
                    lat = atlas$latitude)
locs <- locs[locs$index %in% ids, ]
bad_locs <- locs[locs$index %in% bad, ]
good_locs <- locs[!(locs$index %in% bad), ]

## set up kge color-scaling
map_to_colors <- function(x) {
  # Values less than 0 are assigned "red"
  if (x < 0 || is.na(x)) {
    return("ff0000")
  }
  # Values between 0 and 1 are scaled linearly from red to green
  green <- max(0, min(1, x / 2))
  return(rgb((0.5 - green) * 8 / 5, green * 8 / 5, 0))                         # nolint
}
colors <- vector(mode = "character", length = nrow(good_locs))

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
                    subtitle = paste("kge :", round(scores[i, "kge"], digits = 4)))                         # nolint
    kernels[[j]] <- plt
    colors[j] <- map_to_colors(scores[i, "kge"])
}


print("building map (this may take ~5 minutes, just wait for the > to appear in the console)...")           # nolint
## build map
suppressMessages(
  map <- leaflet()    %>%
  addTiles()      %>%
  addCircleMarkers(data = good_locs,
                radius = 1,
                group = "good",
                color = colors,
                opacity = 1) %>%
  addCircleMarkers(data = bad_locs,
                radius = 1,
                group = "bad",
                color = "#d50000dc",
              opacity = 0.7,
              popup = sapply(bad_locs$index, as.character)) %>%
  addPopupGraphs(kernels, group = "good", width = 300, height = 300)
)

print("done building map. type 'map' in the console to view.")
