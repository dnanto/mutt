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
  
  # map
  
  rv <- reactiveValues(cds = NULL)
  
  observeEvent(root(), {
    file.exists(path <- file.path(root(), "rec.gbk"))
    loci <- c("")
    if (file.exists(path)) loci <- pull(read_gbk_loc(path), accn)
    updateSelectizeInput(session, "accession", choices = loci, selected = NA)
    updateSelectizeInput(session, "product", choices = c(""), selected = NA)
    updateSliderInput(session, "range", max = len(), value = c(1, len()))
  })
  
  observeEvent(input$accession, {
    rv$cds <- NULL
    if (file.exists(path <- file.path(root(), "rec.gbk"))) 
    {
      cds <- 
        genbankr::readGenBank(text = read_gbk_acc(path, input$accession)) %>% 
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
  
  plot_map_func <- function(vcf, input) {
    vcf <- vcf() %>% filter(type %in% input$types)
    
    title <- ggtitle(paste("Variant Calls: ", taxa()[1, "id"], taxa()[1, "definition"]))
    bold <- theme(plot.title = element_text(face = "bold"))
    plt <- plot_map(vcf, input) + xlim(input$range[1], input$range[2])
    
    if (is_empty(cds <- rv$cds))
    {
      plt <- plt + title + bold
    }
    else
    {
      p <- plot_cds(cds, input) + xlim(input$range[1], input$range[2]) + title + bold
      plt <- cowplot::plot_grid(p, plt, ncol = 1, align = "v", rel_heights = input$rel_heights)
    }
    
    plt
  }
  
  output$map <- renderPlot(plot_map_func(vcf(), input))
  
  output$export <- downloadHandler(
    filename = function() paste(taxa()[1, "id"], input$ext, sep="."),
    content = function(file) 
      ggsave(
        file, plot = plot_map_func(vcf(), input), 
        width = input$width, height = input$height, units = input$units
      )
  )
  
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
  
  # conf
  
  output$shape <- renderPlot({
    pch <- c(0:25, 32:127)
    l <- length(pch)
    n <- ceiling(sqrt(l))
    data.frame(pch = pch) %>% 
      rowid_to_column() %>% 
      mutate(col = ceiling(rowid / 12), row = rep(0:(n-1), n)[1:l] + 1) %>%
      ggplot(aes(row, col)) +
      geom_point(aes(shape = pch, size = 14)) + 
      geom_text(aes(label = pch), nudge_y = 0.5) + 
      scale_shape_identity() + 
      scale_y_reverse() + 
      theme_void() +
      theme(legend.position = "none", plot.title = element_text(face = "bold")) +
      ggtitle("shape key") 
  })
  
})
