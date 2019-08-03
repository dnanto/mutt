shinyUI(
  tagList(
    tags$head(
      tags$style(
        "#log { font-family: monospace; overflow-y: scroll; }"
      )
    ),
    navbarPage(
      "mutt",
      tabPanel(
        "data",
        shinyDirButton("import", "Import", "Please select a directory...", F),
        verbatimTextOutput("root"),
        DT::DTOutput("taxa")
      ),
      tabPanel(
        "rtt",
        selectInput("objective", "objective", objective, objective[1]),
        selectInput("distmethod", "distmethod", distmethod, distmethod[1]),
        plotOutput("tree_rtt"),
        plotOutput("plot_rtt")
      )
    )
  )
)
