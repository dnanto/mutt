shinyServer(function(input, output, session) {
  
  # taxa
  
  shinyDirChoose(input, "import", roots = roots)
  
  root <- eventReactive(input$import, {
    files <- parseDirPath(roots, input$import)
    if (length(files)) files
  })
  
  taxa <- eventReactive(root(), {
    req(file.exists(path <- file.path(root(), "msa.fas")))
    tokens <- names(ape::read.FASTA(path)) %>% str_split_fixed(" ", 2)
    data <- setNames(as.data.frame(tokens), c("id", "definition"))
    meta <- file.path(root(), "rec.tsv")
    if (file.exists(meta)) data <- merge(data, read_tsv(meta), all.x = T)
    data
  })
  
  len <- eventReactive(root(), {
    req(file.exists(path <- file.path(root(), "msa.fas")))
    rec <- Biostrings::readDNAStringSet(path, nrec = 1)
    length(rec[[1]]) - Biostrings::countPattern("-", rec[[1]])
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
  
  rv <- reactiveValues(cds = NULL)
  
  observeEvent(root(), {
    file.exists(path <- file.path(root(), "rec.gbk"))
    loci <- c("")
    if (file.exists(path)) loci <- pull(read_gbk_loc(path), accn)
    updateSelectizeInput(session, "reference", choices = loci, selected = NA)
    updateSelectizeInput(session, "product", choices = c(""), selected = NA)
    updateSliderInput(session, "range", max = len(), value = c(1, len()))
  })
  
  observeEvent(input$reference, {
    rv$cds <- NULL
    if (file.exists(path <- file.path(root(), "rec.gbk"))) 
    {
      cds <- 
        genbankr::readGenBank(text = read_gbk_acc(path, input$reference)) %>% 
        genbankr::cds()
      updateSelectizeInput(session, "product", choices = cds$product, selected = NA)
      rv$cds <- cds
    }
  })
  
  observeEvent(input$product, {
    cds <- rv$cds %>% as.data.frame()
    idx <- which(cds$product == input$product)
    updateSliderInput(session, "range", value = c(min(cds$start[idx]), max(cds$end[idx])))
  })
  
  vcf <- eventReactive(input$run, {
    path <- file.path(root(), "msa.fas")
    withProgress({
      incProgress(1/5, "read multiple alignment...")
      msa <- Biostrings::readDNAMultipleAlignment(path)
      lab <- names(msa@unmasked) %>% str_split_fixed(" ", 2) %>% lapply(head, 1)
      mat <- str_split_fixed(msa, "", ncol(msa))
      pos <- positions(mat[1,])
      
      incProgress(1/5, "calling indels...")
      ind <- call_ind(mat)
      
      incProgress(1/5, "calling SNPs...")
      snp <- call_snp(mat)
      
      incProgress(1/5, "combine...")
      bind_rows(ind, snp) %>%
        mutate(
          id = factor(lab[idx], levels = rev(unique(lab))),
          POS = pos[POS]
        ) %>%
        mutate_at("type", as.factor) %>%
        mutate_at("type", factor, levels = c("ins", "del", "trv", "trs"))
    })
  })
  
  output$gmm <- renderPlot({
    vcf <- vcf()
    plt <- plot_gmm(vcf, input) + xlim(input$range[1], input$range[2])
    if (!is_empty(cds <- rv$cds))
      plt <- cowplot::plot_grid(
        plot_cds(cds) + xlim(input$range[1], input$range[2]), plt, 
        ncol = 1, align = "v", rel_heights = c(1, 4)
      )
    
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
