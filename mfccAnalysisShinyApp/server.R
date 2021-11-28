server <- function(input, output) {
  listOfWaves <- reactive({
    "../../machineLearningData/gravitationalWaves/train/file_labels.csv" %>% getListOfWaves()
  })

  oneRandomWaveFilePath <- reactive( {
    listOfWaves() %>%
      filter( target == input$classTypeRadioId) %>%
      sample_n( size = 1) %>%
      select(-X) %>%
      select(filePath) %>%
      slice(1) %>%
      unlist() %>%
      as.character() %>%
      paste0("../",.)
  } )


  output$mfccPlot1 <- renderPlot({
    oneRandomWaveFilePath( ) %>%
      plotOneWaveMFCC2(  1 , input$classTypeRadioId, ncep = input$ncep ,wl = 2^(input$wl), fbtype=input$fbtype, dcttype=input$dcttype )
  })

  output$mfccPlot2 <- renderPlot({
    oneRandomWaveFilePath( ) %>%
      plotOneWaveMFCC2(  2 , input$classTypeRadioId, ncep = input$ncep ,wl = 2^(input$wl), fbtype=input$fbtype, dcttype=input$dcttype )
  })

  output$mfccPlot3 <- renderPlot({
    oneRandomWaveFilePath( ) %>%
      plotOneWaveMFCC2(  3 , input$classTypeRadioId, ncep = input$ncep ,wl = 2^(input$wl), fbtype=input$fbtype, dcttype=input$dcttype )
  })

  output$text <- renderText({" Calculation of the MFCCs imlcudes the following steps:
      1. Preemphasis filtering/n
      2. Take the absolute value of the STFT (usage of Hamming window)
      3. Warp to auditory frequency scale (Mel/Bark)
      4. Take the DCT of the log-auditory-spectrum
      5. Return the first ‘ncep’ components\n" })


}