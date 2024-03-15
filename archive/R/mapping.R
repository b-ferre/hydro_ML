library(maps)
library(ggplot2)
library(here)

# get map data
world_map <- map_data("world")
na_map <- subset(world_map, ((world_map$lat > 0) & (world_map$lat < 100)
                            & (world_map$long > -150) & (world_map$long < -50)))

# create a base plot with gpplot2
p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")

# add world map to base plot
base_world <- p + geom_polygon(data = world_map, aes(x=long, y=lat, group=group),       # nolint
                               colour = "light green", fill = "light green")
base_na <- p + geom_polygon(data = na_map, aes(x=long, y=lat, group=group),
                               colour = "light green", fill = "light green")

# get catchment locations
load(here("data", "raw_data", "catchment_atlas.Rda"))

# plot ALL catchment locations
catchment_map <-
  base_world +
  geom_point(data = catchment_atlas,
             aes(x = longitude, y = latitude), colour = "Deep Pink",
             fill = "Pink", pch = 21, size = 0.2, alpha = I(0.7))

# get NA catchments
na_catchments <- catchment_atlas[(catchment_atlas$country == "CA" |
                                  catchment_atlas$country == "US"), ]

# plot NA catchment locations
na_catchment_map <-
  base_na +
  geom_point(data = na_catchments,
             aes(x = longitude, y = latitude), colour = "Deep Pink",
             fill = "Pink", pch = 21, size = 0.5, alpha = I(0.7))
