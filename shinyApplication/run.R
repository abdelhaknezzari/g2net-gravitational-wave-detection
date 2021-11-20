library(shiny)

port <- Sys.getenv('PORT')

shiny::runApp(
  appDir = "shinyApplication",
  host = '0.0.0.0',
  port = 8085
)
