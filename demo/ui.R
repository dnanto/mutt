fluidPage(
  wellPanel(
    shinyFilesButton(
      "import", "Import MSA", "Please select a multiple sequence alignment FASTA file...", 
      multiple = F, icon = icon("file-import")
    ),
    verbatimTextOutput("path", placeholder = T)
  ),
  tagAppendAttributes(
    textAreaInput("msa", "msa", height = "200px", resize = "vertical"), 
    style = "width: 100%; font-family: monospace;"
  ),
  withSpinner(plotOutput("plot"), type = 8),
  withSpinner(DT::DTOutput("call"), type = 8)
)
