shinyServer(function(input, output, session) {
  
  ## msa
  
  shinyFileChoose(input, "import", roots = roots)
  
  path <- eventReactive(input$import, pull(parseFilePaths(roots, input$import), datapath))
  
  output$path <- eventReactive(path(), { req(nchar(path <- path()) > 0); path; })
  
  raw <- eventReactive(path(), { req(nchar(path <- path()) > 0); read_lines(path); })
  
  observeEvent(raw(), {
    req(lines <- raw())
    lines <- lines[grep("^[^>]", lines)] %>% paste(collapse = "\n")
    updateTextAreaInput(session, "msa", value = lines)
  })
  
  result <- eventReactive(input$msa, {
    req(input$msa) %>%
      str_split("\n", simplify = T) %>% 
      str_split("", simplify = T) %>%
      as.matrix() %>%
      set_rownames(c("ref", paste("alt", 1:(nrow(.)-1), sep = "-"))) %>%
      trace()
  })
  
  output$plot <- renderPlot({
    req(result <- result())
    result$plot
  })
  
  output$call <- DT::renderDT({
    req(result <- result()) 
    DT::datatable(
      select(result$call, id, mut, REF, ALT, POS),
      rownames = FALSE,
      style = "bootstrap",
      class = "table-bordered table-striped table-hover responsive",
      filter = list(position = "top")
    )
  })
  
})
