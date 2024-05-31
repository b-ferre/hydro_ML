require(leaflet)
require(leafpop)
require(ggplot2)
require(sf)
require(here)
require(stringr)


## TODO: add code that plots snow dominated catchments


map <- function() {
  ## load catchment atlas
  load("./data/atlas.Rda")

  ## get list of all catchment numbers I have a kernel for
  extract_number <- function(filename) {
    return(as.integer(str_remove_all(filename, "[^\\d]")))
  }
  plotted_catchments <- lapply(list.files("./data/plots"), extract_number)

  ## pre-process data for mapping
  locs <- data.frame(index = atlas$gridcode,
                      lng = atlas$longitude,
                      lat = atlas$latitude)

  ## trim locs such that only plotted catchments are included
  locs <- locs[locs$index %in% plotted_catchments, ]

  ## create list of kernels which is in the same order as locs
  ## so that binding works correctly
  kernels <- c()
  for (i in length(locs$index)) {
    kern_i <- paste0("./data/plots/beta_est_", locs$index[i], ".jpeg")
    kernels <- append(kernels, kern_i)
  }

  ##TODO: add code which creates colors based on scores

  ## build map
  t0 <- Sys.time()
  print("building map (this may take ~5 minutes, just wait for the > to appear in the console)...")               # nolint
  #suppressMessages(
    map <- leaflet()    %>%
    addTiles()      %>%
    addCircleMarkers(data = good_locs,
                  radius = 3,
                  group = "good",
                  #color = colors,
                  opacity = 1) %>%
    addPopupImages(kernels, group = "good", width = 300, height = 300)                                            # nolint
  #)
  tf <- Sys.time()
  print(paste("elapsed time:", tf - t0))


  return(map)
}