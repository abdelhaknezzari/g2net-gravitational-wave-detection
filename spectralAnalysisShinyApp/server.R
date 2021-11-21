server <- function(input, output) {
  output$specrogramPlot1 <- renderPlot({
    paste0( "data/",input$waveName) %>% readWaveFromNpy(1) %>% showPlots(input$fmax/100, 2^(input$wl), input$wn, input$ovlp )
  })
  output$specrogramPlot2 <- renderPlot({
    paste0( "data/",input$waveName) %>% readWaveFromNpy(2) %>% showPlots(input$fmax/100, 2^(input$wl), input$wn, input$ovlp )
  })

  output$specrogramPlot3 <- renderPlot({
    paste0( "data/",input$waveName) %>% readWaveFromNpy(3) %>% showPlots(input$fmax/100, 2^(input$wl), input$wn, input$ovlp )
  })

}