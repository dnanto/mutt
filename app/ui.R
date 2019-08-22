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
          selectizeInput("accession", "accession", NULL),
          selectizeInput("product", "product", NULL),
          sliderInput("range", "range", min = 1, max = 1, value = c(1, 1)),
          sliderInput("rel_heights", "rel_heights", min = 1, max = 10, value = c(1, 4), step = 0.125, dragRange = F),
          sliderInput("alpha", "alpha", min = 0, max = 1, value = 0.5),
          checkboxGroupInput("types", "types", choices = types, choiceValues = types, selected = types, inline = T),
          checkboxGroupInput("labels", "labels", choices = labels, choiceValues = labels, selected = labels[1], inline = T),
          numericInput("label_size", "label_size", 4, min = 0),
          actionButton("run", "run", width = "100%"),
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
              numericInput("shape_trv", "shape_trv", val = 124, max = 127),
              numericInput("shape_trs", "shape_trs", val = 124, max = 127),
              numericInput("shape_ins", "shape_ins", val = 6, max = 127),
              numericInput("shape_del", "shape_del", val = 2, max = 127)
            )
          )
        ),
        plotOutput("shape")
      )
    ),
    tabPanel(
      "export",
      fluidPage(
        wellPanel(
          selectInput("ext", "extension", choices = ext, selected = "pdf"),
          selectInput("units", "units", choices = units, selected = "in"),
          numericInput("width", "width", 8, min = 0),
          numericInput("height", "height", 10, min = 0)
        ),
        downloadButton("export", "export")
      )
    )
  )
)

