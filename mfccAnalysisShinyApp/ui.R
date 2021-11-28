ui <- fluidPage(

  # App title ----
  titlePanel("MFCC analysis"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      textOutput("text"),
      radioButtons("classTypeRadioId", "Choose the class type:",
                   c( "Black Hole not Detected" = 0, "Black Hole Detected" = 1)
                   ),
     selectInput("fbtype", "Auditory frequency scale to use:",
                  c("htkmel" = "htkmel",
                    "mel" = "mel",
                    "bark" = "bark",
                    "fcmel" = "fcmel")),

      selectInput("dcttype", "Discrete cosine transform type:",
                  c("t1" = "t1",
                    "t2" = "t2",
                    "t3" = "t3",
                    "t4" = "t4")),

      sliderInput("ncep", "Number of cepstra to return:",
                  min = 5, max = 35,
                  value = 13, step = 1),
      sliderInput("wl", "Window length 2^n specify n",
                  min = 1, max = 20,
                  value = 10, step = 1),
  ),
    # Main panel for displaying outputs ----
    mainPanel(
      # h3("MFCC coeficients"),
      plotOutput(outputId = "mfccPlot1",width = "100%"),
      plotOutput(outputId = "mfccPlot2",width = "100%"),
      plotOutput(outputId = "mfccPlot3",width = "100%")
    )
)
)