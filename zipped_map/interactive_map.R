library(leaflet)
library(leafpop)
library(here)
library(stringr)
library(ggplot2)
library(sf)

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
dryness_pallete <- colorNumeric(palette = "plasma", domain = scores$kge)


## build kernel graphs for good catchments
kernels <- list()
for (i in seq_len(nrow(good_locs))) {
    data <- data.frame(sapply(read.csv(paste("./FD_kernels/", good_locs$index[i], ".csv", sep = "")), as.numeric))      # nolint
    colnames(data) <- c("lag", "IRF_coeff")
    plt <- ggplot(data = data, aes(x = lag, y = IRF_coeff)) +
                geom_point() +
                geom_line(linetype = "dashed", color = "red") +
                xlim(0, max(50, nrow(data))) +
                labs(title = paste("gridcode", good_locs$index[i]),
                    subtitle = paste("'kge' :", scores[good_locs$index[i]]))                         # nolint
    kernels[[i]] <- plt
}

## get model performances

map <- leaflet()    %>%
    addTiles()      %>%
    #setView(lat = 54.074372, lng = -104.535173, zoom = 3) %>%
    addCircleMarkers(data = good_locs,
                  radius = 1,
                  group = "pt",
                  color = ~dryness_pallete(reactivities),
                  opacity = 1) %>%
    ## uncomment lines below to show "bad" catchment locations, colored red
    # addCircleMarkers(data = bad_locs,
    #              radius = 0.5,
    #              group = "bad",
    #              color = "red",
    #             opacity = 0.2) %>%
    addPopupGraphs(kernels, group = "pt", width = 300, height = 300)

map