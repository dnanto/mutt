shinyServer(function(input, output, session) {
  
  ## msa
  
  shinyFileChoose(input, "import", roots = roots)
  
  path <- eventReactive(input$import, pull(parseFilePaths(roots, input$import), datapath))
  
  output$path <- eventReactive(path(), { req(nchar(path <- path()) > 0); path; })
  
  raw <- eventReactive(path(), { req(nchar(path <- path()) > 0); read_lines(path); })
  
  observeEvent(raw(), {
    req(lines <- raw())
    lines <- lines[grep("^[^>]", lines)] %>% paste(collapse = "\n")
    updateTextAreaInput(session, "box", value = lines)
  })

  output$plot <- renderPlot({
    req(input$box)
    str_split(input$box, "\n", simplify = T) %>% 
      str_split("", simplify = T) %>%
      plot_trace()
  })
  
})
