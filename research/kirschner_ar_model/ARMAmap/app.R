library(shiny)
library(leaflet)

# Load necessary R files
source("R/interactive_map.R")

# Define UI
ui <- fluidPage(
  # Leaflet map output
  leafletOutput("map")
)

# Define server logic
server <- function(input, output, session) {
  # Render the leaflet map
  output$map <- renderLeaflet({
    # Call your map function and store the leaflet object
    map_object <- map()
    # Return the leaflet object
    map_object
  })
}

# Run the application
shinyApp(ui, server)
