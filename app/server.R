shinyServer(function(input, output, session) {
  
  ## msa
  
  shinyFileChoose(input, "import_msa", roots = roots)
  
  path_msa <- eventReactive(input$import_msa, pull(parseFilePaths(roots, input$import_msa), datapath))
  
  output$path_msa <- eventReactive(path_msa(), { req(nchar(path <- path_msa()) > 0); path; })
  
  observeEvent(path_msa(), updateSliderInput(session, "range", max = len(), value = c(1, 0.15 * len())))
  
  msa <- eventReactive(path_msa(), { req(nchar(path <- path_msa()) > 0); readDNAMultipleAlignment(path); })
  
  len <- eventReactive(msa(), { req(msa <- msa()); ncol(msa) - countPattern("-", msa@unmasked[[1]]); })
  
  output$taxa <- DT::renderDT({
    req(msa <- msa())
    names(msa@unmasked) %>% 
      enframe(name = NULL) %>% 
      separate(value, c("id", "definition"), sep = " ", extra = "merge", fill = "right") %>%
      DT::datatable(
        rownames = FALSE,
        style = "bootstrap",
        class = "table-bordered table-striped table-hover responsive",
        filter = list(position = "top")
      )
  })
  
  ## gbk
  
  shinyFileChoose(input, "import_gbk", roots = roots)
  
  path_gbk <- eventReactive(input$import_gbk, pull(parseFilePaths(roots, input$import_gbk), datapath))
  
  output$path_gbk <- eventReactive(path_gbk(), { req(nchar(path <- path_gbk()) > 0); path; })
  
  observeEvent(path_gbk(), {
    req(nchar(path <- path_gbk()) > 0)
    loci <- pull(read_gbk_loc(path), accn)
    updateSelectizeInput(session, "accession", choices = loci, selected = NA)
    updateSelectizeInput(session, "product", choices = c(""), selected = NA)
  })
  
  cds <- eventReactive(input$accession, {
    req(nchar(path <- path_gbk()) > 0)
    withProgress({
      incProgress(1/3, "reading genbank...")
      genbankr::readGenBank(text = read_gbk_acc(path, input$accession), ret.seq = F) %>% genbankr::cds()
    })
  })
  
  output$cds <- DT::renderDT(
    DT::datatable(
      select(as.data.frame(cds()), -translation),
      rownames = FALSE,
      style = "bootstrap",
      class = "table-bordered table-striped table-hover responsive",
      filter = list(position = "top"),
      options = list(scrollX = TRUE)
    )
  )
  
  observeEvent(cds(), updateSelectizeInput(session, "product", choices = cds()$product, selected = NA))
  
  observeEvent(input$product, {
    cds <- cds() %>% as.data.frame()
    idx <- which(cds$product == input$product)
    updateSliderInput(session, "range", value = c(min(cds$start[idx]), max(cds$end[idx])))
  })

  ## calls
  
  observeEvent(input$run, updateTabsetPanel(session, "tabset_panel", selected = "tab_map"))
  
  calls <- eventReactive(input$run, {
    withProgress({
      incProgress(1/5, "read multiple alignment...")
      msa <- msa()
      lab <- names(msa@unmasked) %>% str_split_fixed(" ", 2) %>% apply(1, head, 1)
      mat <- str_split_fixed(msa, "", ncol(msa))
      pos <- positions(mat[1,])
      
      incProgress(1/5, "calling indels...")
      ind <- call_ind(mat)
      
      incProgress(1/5, "calling SNPs...")
      snp <- call_snp(mat)
      
      incProgress(1/5, "combine...")
      mutate(
        bind_rows(ind, snp),
        POS = pos[POS],
        id = factor(lab[idx], levels = rev(unique(lab))),
        type = factor(as.factor(type), levels = c("ins", "del", "trv", "trs"))
      )
    })
  })
  
  plot_map <- function() 
  {
    calls <- filter(calls(), type %in% input$types & input$range[1] <= POS & POS <= input$range[2])
    
    title <- ggtitle(paste("Variant Calls [", names(msa()@unmasked)[1], "]"))
    
    color <- setNames(c(input$color_trv, input$color_trs, input$color_ins, input$color_del), types)
    alpha <- setNames(c(input$alpha_trv, input$alpha_trs, input$alpha_ins, input$alpha_del), types)
    shape <- setNames(c(input$shape_trv, input$shape_trs, input$shape_ins, input$shape_del), types)
    size <- setNames(c(input$size_trv, input$size_trv, input$size_ins, input$size_del), types)
    
    plt <-
      ggplot(calls, aes(POS, id)) +
      geom_point(aes(color = type, alpha = type, size = type, shape = type)) +
      scale_color_manual(values = color) +
      scale_alpha_manual(values = alpha) +
      scale_size_manual(values = size) +
      scale_shape_manual(values = shape) +
      xlab("pos") +
      theme_minimal() +
      theme(legend.position = "bottom", title = element_text(size = input$title_size))
    
    cds <- NULL
    
    if (nchar(input$accession) > 0)
    {
      cds <- cds()
      cds <- cds[input$range[1] <= start(cds) & start(cds) <= input$range[2]]
    }
    
    if (is_empty(cds)) 
    {
      plt <- plt + title
    }
    else
    {
      lvl <- overlevels(cds) %>% factor(., levels = rev(unique(.)))
      plt_cds <-
        data.frame(
          product = cds$product, strand = cds@strand, lvl = lvl,
          x = cds@ranges@start, xend = cds@ranges@start + cds@ranges@width - 1
        ) %>%
          transform(x = ifelse(strand == '-', xend, x), xend = ifelse(strand == '+', xend, x)) %>%
          ggplot(aes(x = x, xend = xend, y = lvl, yend = lvl)) +
            geom_segment(
              aes(color = product), 
              lineend = input$lineend, linejoin = input$linejoin, size = input$segment_size,
              arrow = arrow(length = unit(input$arrow_length, input$arrow_units), type = input$arrow_type)
            ) +
          {
            if (input$label)
              geom_text_repel(
                aes(label = product, x = (x + xend) / 2), 
                hjust = "left", vjust = "bottom", nudge_y = 0.25,
                size = input$text_size_cds, segment.size = input$line_size_cds
              )
          } +
          theme_minimal() +
          theme(legend.position = "none") +
          xlab("") + 
          ylab("product") +
          theme_minimal() +
          theme(
            title = element_text(size = input$title_size),
            legend.position = "none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank()
          ) +
          title
      
      plt <- cowplot::plot_grid(plt_cds, plt, ncol = 1, align = "v", rel_heights = c(input$rel_cds, input$rel_msa))
    }
    
    plt
  }
  
  observeEvent(input$height, output$map <- renderPlot(plot_map(), height = input$height))
  
  output$export <- downloadHandler(
    filename = function() paste(str_split_fixed(names(msa()@unmasked)[1], " ", 2)[,1], input$ext, sep="."),
    content = function(file) 
      ggsave(
        file, plot = plot_map(), 
        width = input$exp_width, height = input$exp_height, 
        units = input$units, scale = input$scale, dpi = input$dpi
      )
  )
  
  # calls
  
  output$calls <- DT::renderDT(
    DT::datatable(
      arrange(select(calls(), id, everything(), -idx), id, POS),
      rownames = FALSE,
      style = "bootstrap",
      class = "table-bordered table-striped table-hover responsive",
      filter = list(position = "top")
    )
  )
  
})
