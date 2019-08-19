shinyServer(function(input, output, session) {
  
  # taxa
  
  shinyDirChoose(input, "import", roots = roots)
  
  root <- eventReactive(input$import, {
    files <- parseDirPath(roots, input$import)
    if (length(files)) files
  })
  
  taxa <- eventReactive(root(), {
    path <- file.path(root(), "msa.fas")
    req(file.exists(path))
    tokens <- names(ape::read.FASTA(path)) %>% str_split_fixed(" ", 2)
    data <- setNames(as.data.frame(tokens), c("id", "definition"))
    meta <- file.path(root(), "rec.tsv")
    if (file.exists(meta)) data <- merge(data, read_tsv(meta), all.x = T)
    data
  })
  
  len <- eventReactive(root(), {
    path <- file.path(root(), "msa.fas")
    req(file.exists(path))
    as.integer(Biostrings::fasta.seqlengths(path, nrec = 1))
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
    req(len <- len(), taxa <- taxa())
    updateSelectizeInput(session, "reference", choices = pull(taxa, "id"), selected = NA)
    updateSliderInput(session, "range", max = len, value = c(1, len))
  })
  
  cds <- eventReactive(input$reference, {
    path <- file.path(root(), "rec.gbk")
    req(input$reference)
    if (file.exists(path))
      genbankr::readGenBank(text = read_gbk_acc(path, input$reference)) %>% genbankr::cds()
  })

  observeEvent(cds(), {
    req(cds <- cds())
    updateSelectizeInput(session, "product", choices = cds$product, selected = NA)
  })
  
  observeEvent(input$product, {
    path <- file.path(root(), "msa.fas")
    req(len <- len(), cds <- cds(), input$product, file.exists(path))
    cds <- as.data.frame(cds)
    idx <- which(cds$product == input$product)
    updateSliderInput(session, "range", max = len, value = c(min(cds$start[idx]), max(cds$end[idx])))
  })
  
  rv <- reactiveValues(vcf = NULL, plt.cds = NULL, plt.gmm = NULL)
  
  vcf <- eventReactive(input$run, {
    path <- file.path(root(), "msa.fas")
    req(ref <- input$reference, file.exists(path))
    withProgress({
      incProgress(1/4, "read multiple alignment...")
      msa <- Biostrings::readDNAMultipleAlignment(path)
      lab <- names(msa@unmasked) %>% str_split(" ", 2) %>% sapply(head, 1)
      mat <- str_split_fixed(msa, "", ncol(msa))
      idx <- grep(ref, lab) %>% c(., setdiff(1:nrow(mat), .))
      mat <- mat[idx, ]
      lab <- factor(lab[idx], levels = lab[idx])
      pos <- positions(mat[1, ])
      
      incProgress(1/4, "calculate indels...")
      ind <- call_ind(mat)
      incProgress(1/4, "calculate SNPs...")
      snp <- call_snp(mat) 
      bind_rows(ind, snp) %>% 
        mutate(id = factor(lab[idx], levels = rev(unique(lab))), POS = pos[POS]) %>%
        mutate_at("type", as.factor) %>%
        mutate_at("type", factor, levels = c("ins", "del", "trv", "trs"))
    })
  })

  output$gmm <- renderPlot({
    req(vcf <- vcf())
    print(!is_empty(cds))
    cds <- cds()
    plt <- plot_gmm(vcf, input) + xlim(input$range[1], input$range[2])
    if (!is_empty(cds))
      plt <- 
        (plot_cds(cds) + xlim(input$range[1], input$range[2])) %>% 
        cowplot::plot_grid(plt, ncol = 1, align = "v", rel_heights = c(1, 4))
    
    plt
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
