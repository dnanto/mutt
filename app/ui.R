shinyUI(
  navbarPage(
    "mutt",
    tabPanel(
      "taxa",
      shinyDirButton("import", "Import", "Please select a directory...", F),
      verbatimTextOutput("root"),
      withSpinner(DT::DTOutput("taxa"), type = 8)
    ),
    tabPanel(
      "gmm",
      selectizeInput("reference", "reference", NULL),
      selectizeInput("product", "product", NULL),
      sliderInput("range", "range", min = 1, max = 1, value = c(1, 1)),
      actionButton("run", "run"),
      withSpinner(plotOutput("gmm", height = 800), type = 8)
    ),
    tabPanel(
      "calls",
      withSpinner(DT::DTOutput("calls"), type = 8)
    ),
    tabPanel(
      "conf",
      fluidPage(
        fluidRow(
          column(
            4,
            spectrumInput("color_trv", "color_trv", selected = "magenta"),
            spectrumInput("color_trs", "color_trs", selected = "cyan"),
            spectrumInput("color_ins", "color_ins", selected = "black"),
            spectrumInput("color_del", "color_del", selected = "black")
          ),
          column(
            4,
            numericInput("size_trv", "size_trv", val = 10, min = 1),
            numericInput("size_trs", "size_trs", val = 10, min = 1),
            numericInput("size_ins", "size_ins", val = 2, min = 1),
            numericInput("size_del", "size_del", val = 2, min = 1)
          ),
          column(
            4,
            numericInput("shape_trv", "shape_trv", val = 124, max = 127),
            numericInput("shape_trs", "shape_trs", val = 124, max = 127),
            numericInput("shape_ins", "shape_ins", val = 6, max = 127),
            numericInput("shape_del", "shape_del", val = 2, max = 127)
          )
        )
      )
    )
  )
)

