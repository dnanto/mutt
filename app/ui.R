library(shinydashboard)

dashboardPage(
  dashboardHeader(title = "mutt"),
  dashboardSidebar(
    sliderInput("range", "range", min = 1, max = 1, value = c(1, 1), width = "100%"),
    fluidRow(
      column(4, numericInput("height", "height", 1000, step = 100)),
      column(4, numericInput("rel_cds", "rel_cds", 1, min = 1)),
      column(4, numericInput("rel_msa", "rel_msa", 4, min = 1))
    ),
    checkboxGroupInput("types", "types", choices = types, choiceValues = types, selected = types, inline = T, width = "100%"),
    selectizeInput("product", "product", NULL, width = "100%"),
    fluidRow(
      column(4, numericInput("text_size_cds", "text_size_cds", 2, min = 0, step = 0.5)),
      column(4, numericInput("line_size_cds", "line_size_cds", 0.25, min = 0, step = 0.25)),
      column(4, checkboxInput("label", "display", value = T))
    ),
    fluidRow(
      column(4, numericInput("arrow_length", "arrow_length", 0.0025, min = 0, step = 0.0025)),
      column(4, selectInput("arrow_units", "arrow_units", arrow_units, arrow_units[1])),
      column(4, selectInput("arrow_type", "arrow_type", arrow_type, arrow_type[2]))
    ),
    fluidRow(
      column(4, selectInput("lineend", "line-end", lineend, lineend[2])),
      column(4, selectInput("linejoin", "line-join", linejoin, linejoin[2])),
      column(4, numericInput("segment_size", "segment_size", 2, min = 0, step = 1))
    ),
    numericInput("alpha", "alpha", 0.5, min = 0, max = 1, step = 0.1, width = "100%"),
    fluidRow(
      column(4, numericInput("size_trv", "size_trv", 10, min = 0)),
      column(4, numericInput("shape_trv", "shape_trv", 124, min = 0, max = 127)),
      column(4, spectrumInput("color_trv", "color_trv", selected = "magenta", update_on = "change"))
    ),
    fluidRow(
      column(4, numericInput("size_trs", "size_trs", 10, min = 0)),
      column(4, numericInput("shape_trs", "shape_trs", 124, min = 0, max = 127)),
      column(4, spectrumInput("color_trs", "color_trs", selected = "cyan", update_on = "change"))
    ),
    fluidRow(
      column(4, numericInput("size_ins", "size_ins", 4, min = 0)),
      column(4, numericInput("shape_ins", "shape_ins", 6, min = 0, max = 127)),
      column(4, spectrumInput("color_ins", "color_ins", selected = "black", update_on = "change"))
    ),
    fluidRow(
      column(4, numericInput("size_del", "size_del", 4, min = 0)),
      column(4, numericInput("shape_del", "shape_del", 2, min = 0, max = 127)),
      column(4, spectrumInput("color_del", "color_del", selected = "black", update_on = "change"))
    ),
    width = 400
  ),
  dashboardBody(
    box(
      title = "Taxa", width = 12, collapsible = T,
      wellPanel(
        shinyFilesButton(
          "import_msa", "Import MSA", "Please select a multiple sequence alignment FASTA file...", 
          multiple = F, icon = icon("file-import")
        ),
        verbatimTextOutput("path_msa", placeholder = T)
      ),
      DT::DTOutput("taxa")
    ),
    box(
      title = "CDS", width = 12, collapsible = T,
      wellPanel(
        shinyFilesButton(
          "import_gbk", "Import GenBank", "Please select a GenBank file...",
          multiple = F, icon = icon("file-import")
        ),
        verbatimTextOutput("path_gbk", placeholder = T)
      ),
      selectizeInput("accession", "accession", NULL, width = "100%"),
      DT::DTOutput("cds")
    ),
    box(
      title = "Export", width = 12, collapsible = T, collapsed = T,
      fluidRow(
        column(6, selectInput("ext", "format", choices = ext, selected = "pdf")),
        column(6, selectInput("units", "units", choices = units, selected = "in"))
      ),
      fluidRow(
        column(6, numericInput("exp_width", "width", 10, min = 0)),
        column(6, numericInput("exp_height", "height", 15, min = 0))
      ),
      fluidRow(
        column(6, numericInput("scale", "scale", 1, min = 0.1, step = 0.1)),
        column(6, numericInput("dpi", "dpi", 300, min = 72, step = 50))
      ),
      wellPanel(downloadButton("export", "Export"))
    ),
    actionButton("run", "Plot", width = "100%"),
    p(),
    withSpinner(plotOutput("map", height = 1000), type = 8)
  )
)