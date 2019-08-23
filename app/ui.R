shinyUI(
  navbarPage(
    "mutt",
    tabPanel(
      "taxa",
      wellPanel(shinyDirButton("import", "Import", "Please select a directory...", F)),
      verbatimTextOutput("root"),
      withSpinner(DT::DTOutput("taxa"), type = 8)
    ),
    tabPanel(
      "map",
      sidebarLayout(
        sidebarPanel(
          actionButton("run", "run", width = "100%"),
          h4("Map Controls"),
          sliderInput("range", "range", min = 1, max = 1, value = c(1, 1)),
          sliderInput("alpha", "alpha", min = 0, max = 1, value = 0.5),
          checkboxGroupInput("types", "types", choices = types, choiceValues = types, selected = types, inline = T),
          numericInput("height", "height (px)", 800, step = 100),
          fluidRow(
            column(6, numericInput("rel1", "rel1", 1, min = 1)),
            column(6, numericInput("rel2", "rel2", 4, min = 1))
          ),
          h4("CDS Controls"),
          selectizeInput("accession", "accession", NULL),
          selectizeInput("product", "product", NULL),
          fluidRow(
            column(6, checkboxGroupInput("labels", "label_display", choices = labels, choiceValues = labels, selected = labels[1], inline = T)),
            column(6, numericInput("label_size", "label_size", 4, min = 0))
          ),
          h4("Export Controls"),
          fluidRow(
            column(6, selectInput("ext", "format", choices = ext, selected = "pdf")),
            column(6, selectInput("units", "units", choices = units, selected = "in"))
          ),
          fluidRow(
            column(6, numericInput("exp_width", "width", 8, min = 0)),
            column(6, numericInput("exp_height", "height", 10, min = 0))
          ),
          fluidRow(
            column(6, numericInput("scale", "scale", 1, min = 0.1, step = 0.1)),
            column(6, numericInput("dpi", "dpi", 300, min = 72, step = 50))
          ),
          downloadButton("export", "export"),
          width = 3
        ), 
        mainPanel(
          withSpinner(plotOutput("map", height = 800), type = 8))
        )
    ),
    tabPanel(
      "calls",
      withSpinner(DT::DTOutput("calls"), type = 8)
    ),
    tabPanel(
      "conf",
      fluidPage(
        wellPanel(
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
              numericInput("shape_trv", "shape_trv", val = 124, min = 0, max = 127),
              numericInput("shape_trs", "shape_trs", val = 124, min = 0, max = 127),
              numericInput("shape_ins", "shape_ins", val = 6, min = 0, max = 127),
              numericInput("shape_del", "shape_del", val = 2, min = 0, max = 127)
            )
          )
        ),
        plotOutput("shape")
      )
    )
  )
)

