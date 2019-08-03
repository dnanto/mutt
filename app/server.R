shinyServer(function(input, output, session) {
  
  # data
  
  shinyDirChoose(input, "import", roots = roots)
  
  root <- eventReactive(input$import, {
    files <- parseDirPath(roots, input$import)
    if (length(files)) files
  })
  
  phy.ml <- reactive({
    path <- file.path(root(), "phy.treefile")
    req(file.exists(path))
    phy <- treeio::read.iqtree(path)
    phy@phylo$tip.date <- parse_tip_date(phy@phylo$tip.label)
    phy
  })
  
  taxa <- reactive({
    phy <- phy.ml(); req(phy);
    str_split_fixed(phy@phylo$tip.label, "_", 2) %>%
      as.data.frame() %>%
      setNames(c("acc", "date"))
  })
  
  output$taxa <- DT::renderDT({ 
    req(taxa <- taxa())
    DT::datatable(
      taxa,
      rownames = FALSE,
      style = "bootstrap",
      class = "table-bordered table-striped table-hover responsive",
      filter = list(position = "top")
    )
  })
  
  # rtt
  
  phy.rtt <- reactive({
    req(phy <- phy.ml())
    tip.date <- phy@phylo$tip.date
    phy <- rtt(phy@phylo, as.numeric(tip.date), objective = input$objective)
    phy$tip.date <- tip.date
    phy
  })

  output$tree_rtt <- renderPlot({ req(phy <- phy.rtt()); plot(phy) })

  output$plot_rtt <- renderPlot({
    req(phy <- phy.rtt())
    data.frame(date = phy$tip.date, rtt = adephylo::distRoot(phy, method = input$distmethod)) %>%
      ggplot(aes(date, rtt)) +
        geom_point() +
        geom_smooth(method = "lm") +
        theme_minimal()
  })
  
})
