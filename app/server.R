shinyServer(function(input, output, session) {
  
  # taxa
  
  shinyDirChoose(input, "import", roots = roots)
  
  root <- eventReactive(input$import, {
    files <- parseDirPath(roots, input$import)
    if (length(files)) files
  })
  
  taxa <- reactive({
    path <- file.path(root(), "rec.fas")
    req(file.exists(path))
    tokens <- names(ape::read.FASTA(path)) %>% str_split_fixed(" ", 2)
    data <- setNames(as.data.frame(tokens), c("id", "definition"))
    meta <- file.path(root(), "rec.tsv")
    if (file.exists(meta)) data <- merge(data, read_tsv(meta), all.x = T)
    data
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
  
  # gmm
  
  observeEvent(taxa(), {
    choices <- pull(taxa(), "id")
    updateSelectInput(session, "reference", choices = choices)
  })
  
  output$reference <- renderUI({
    req(taxa <- taxa())
    choices <- pull(taxa, "id")
    selectInput("reference", "reference", choices, choices[1])
  })

  data <- eventReactive(input$run, {
    root <- root()
    path <- file.path(root, "msa.fas")
    req(ref <- input$reference)
    
    plt <- NULL
    txt <- file.path(root, "gmm.txt")
    log <- file.path(root, "gmm.log")
    if (file.exists(txt)) file.remove(txt)
    if (file.exists(log)) file.remove(log)
    file.create(log)
    
    withProgress({
      cmd <- paste("./gmm.sh", path, ref, wait = F)
      system(cmd, wait = F)
      
      while(file.exists(log)) try(incProgress(message = read_file(log)), silent = T)
      
      incProgress(message = "plotting...")
      files <- pull(read_tsv(txt, col_names = "file"), "file")
      
      len <- 
        files[which.max(file.size(files))] %>%
        read_lines() %>% 
        .[grep("^##contig=", .)] %>%
        str_match(., pattern = "length=(\\d+)") %>%
        .[,2] %>%
        as.integer()
      
      data <- list(vcf = calls(files), len = len)
      path <- file.path(root, "rec.gbk")
      if (file.exists(path))
        data$gbk <- genbankr::readGenBank(text = read_gbk_acc(path, input$reference))
    })
    
    data
  })
  
  output$gmm <- renderPlot({
    isolate(input$reference)
    
    req(data <- data())
    
    gbk <- data$gbk
    cds <- genbankr::cds(gbk)
    
    lvl <- overlevels(cds)
    lvl <- factor(lvl, levels = rev(unique(lvl)))
    
    p1 <-
      data.frame(
        x = cds@ranges@start, 
        xend = cds@ranges@start + cds@ranges@width - 1, 
        product = cds$product,
        strand = cds@strand,
        size = 10,
        lvl = lvl
      ) %>%
      ggplot(aes(x = x, xend = xend, y = lvl, yend = lvl, size = size)) +
        geom_segment(aes(color = product)) +
        xlab("") + 
        ylab("") +
        xlim(1, data$len) +
        theme_minimal() +
        theme(legend.position = "none")
    
    types <- c("trv", "trs", "ins", "del")
    color <- setNames(c(input$color_trv, input$color_trs, input$color_ins, input$color_del), types)
    size <- setNames(c(input$size_trv, input$size_trv, input$size_ins, input$size_del), types)
    shape <- setNames(c(input$shape_trv, input$shape_trs, input$shape_ins, input$shape_del), types)
    
    p2 <-
      ggplot(data$vcf, aes(POS, id)) +
      geom_point(aes(color = type, size = type, shape = type)) +
      scale_color_manual(values = color) +
      scale_size_manual(values = size) +
      scale_shape_manual(values = shape) +
      xlim(1, data$len) +
      xlab("pos") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 4))
  })
  
  # calls
  
  output$calls <- DT::renderDT({
    req(data <- data())
    DT::datatable(
      data$vcf,
      rownames = FALSE,
      style = "bootstrap",
      class = "table-bordered table-striped table-hover responsive",
      filter = list(position = "top")
    )
  })
  
})
