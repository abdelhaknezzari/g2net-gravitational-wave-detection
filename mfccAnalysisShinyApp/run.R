library(shiny)

port <- Sys.getenv('PORT')

shiny::runApp(
  appDir = "mfccAnalysisShinyApp",
  host = '0.0.0.0',
  port = 8085
)
