shinyServer(function(input, output, session) {
  
  # taxa
  
  shinyDirChoose(input, "import", roots = roots)
  
  root <- eventReactive(input$import, {
    files <- parseDirPath(roots, input$import)
    if (length(files)) files
  })
  
  taxa <- eventReactive(root(), {
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
    req(taxa <- taxa())
    updateSelectInput(session, "reference", choices = pull(taxa, "id"))
  })
  
  cds <- eventReactive(input$reference, {
    path <- file.path(root(), "rec.gbk")
    req(file.exists(path), input$reference)
    genbankr::readGenBank(text = read_gbk_acc(path, input$reference)) %>% genbankr::cds()
  })

  observeEvent(cds(), {
    req(cds <- cds())
    updateSelectInput(session, "product", choices = cds$product)
  })
  
  observeEvent(input$product, {
    path <- file.path(root(), "msa.fas")
    
    req(cds <- cds(), input$product, file.exists(path))
    cds <- as.data.frame(cds)
    
    len <- as.integer(Biostrings::fasta.seqlengths(path, nrec = 1))
    idx <- which(cds$product == input$product)
    updateSliderInput(session, "range", max = len, value = c(min(cds$start[idx]), max(cds$end[idx])))
  })
  
  rv <- reactiveValues(vcf = NULL, plt.cds = NULL, plt.gmm = NULL)
  
  vcf <- eventReactive(input$run, {
    root <- root()
    msa <- file.path(root, "msa.fas")
    req(ref <- input$reference, file.exists(msa))

    txt <- file.path(root, "gmm.txt")
    log <- file.path(root, "gmm.log")
    if (file.exists(txt)) file.remove(txt)
    if (file.exists(log)) file.remove(log)
    file.create(log)

    withProgress({
      cmd <- paste("./gmm.sh", msa, ref, wait = F)
      system(cmd, wait = F)

      while(file.exists(log)) try(incProgress(message = read_file(log)), silent = T)

      incProgress(message = "process vcf files...")
      calls(pull(read_tsv(txt, col_names = "file"), "file"))
    })
  })

  output$gmm <- renderPlot({
    req(cds <- cds(), vcf <- vcf())
    print("renderPlot")
    plt.cds <- plot_cds(cds) + xlim(input$range[1], input$range[2])
    plt.gmm <- plot_gmm(vcf, input) + xlim(input$range[1], input$range[2])
    cowplot::plot_grid(plt.cds, plt.gmm, ncol = 1, align = "v", rel_heights = c(1, 4))
  })

  # calls

  output$calls <- DT::renderDT({
    req(vcf <- vcf())
    DT::datatable(
      vcf,
      rownames = FALSE,
      style = "bootstrap",
      class = "table-bordered table-striped table-hover responsive",
      filter = list(position = "top")
    )
  })
  
})
