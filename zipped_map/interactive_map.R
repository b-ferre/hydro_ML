require(leaflet)
require(leafpop)
require(here)
require(stringr)
require(ggplot2)
require(sf)

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
ii <- cut(scores[!(scores$id %in% bad), "kge"],
          breaks = seq(-1, 1, len = nrow(good_locs)),
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("Red", "Green"), bias = 5)(nrow(good_locs))[ii]

## build kernel graphs for good catchments
kernels <- list()
for (i in good_locs$index) {
    data <- data.frame(sapply(read.csv(paste("./data/kernels/", i, ".csv", sep = "")), as.numeric))      # nolint
    colnames(data) <- c("lag", "IRF_coeff")
    plt <- ggplot(data = data, aes(x = lag, y = IRF_coeff)) +
                geom_point() +
                geom_line(linetype = "dashed", color = "#2200ff") +
                xlim(0, max(30, nrow(data))) +
                labs(title = paste("gridcode", i),
                    subtitle = paste("kge :", round(scores[i, "kge"], digits = 4)))                         # nolint
    kernels[[length(kernels) + 1]] <- plt
}

## build map
map <- leaflet()    %>%
    addTiles()      %>%
    addCircleMarkers(data = good_locs[0:100, ],
                  radius = 1,
                  group = "good",
                  color = colors[0:100],
                  opacity = 1) %>%
    addCircleMarkers(data = bad_locs,
                 radius = 0.5,
                 group = "bad",
                 color = "#323232b7",
                opacity = 0.2,
                popup = sapply(bad_locs$index, as.character)) %>%
    addPopupGraphs(kernels[0:100], group = "good", width = 300, height = 300)

map