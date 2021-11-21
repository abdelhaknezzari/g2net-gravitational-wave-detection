ui <- fluidPage(

  # App title ----
  titlePanel("Spectral analysis"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      selectInput("waveName",
                  "The name of the wave:",
                  c("00000e74ad.npy" = "00000e74ad.npy",
                    "00001f4945.npy" = "00001f4945.npy",
                    "0860082986.npy" = "0860082986.npy",
                    "0c0019670b.npy" = "0c0019670b.npy",
                    "0c00123722.npy" = "0c00123722.npy" ,
                    "0c0008c89c.npy" = "0c0008c89c.npy")),


      selectInput("wn",
                  "Window name:",
                  c("hanning" = "hanning",
                    "bartlett" = "bartlett",
                    "blackman" = "blackman",
                    "flattop" = "flattop",
                    "hamming" = "hamming" ,
                    "rectangle" = "rectangle")),

      sliderInput(inputId = "fmax",
                  label = "Uper frequency:",
                  min = 1,
                  max = 100,
                  value = 30),


      sliderInput(inputId = "wl",
                  label = "Frequency window length = 2^n, choose n:",
                  min = 7,
                  max = 15,
                  value = 11),

      sliderInput(inputId = "ovlp",
                  label = "Choose overlap in %",
                  min = 0,
                  max = 100,
                  value = 11)

    ),

    # Main panel for displaying outputs ----
    mainPanel(
      h3("First detector:"),
      plotOutput(outputId = "specrogramPlot1"),
      h3("Second detector:"),
      plotOutput(outputId = "specrogramPlot2"),
      h3("Third detector:"),
      plotOutput(outputId = "specrogramPlot3"),

    )
  )
)