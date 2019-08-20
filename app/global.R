library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinycssloaders)
library(tidyverse)

roots = c(home = "..")

read_gbk_loc <- function(path)
{
  read_lines(path) %>%
    .[grep("^LOCUS", .)] %>%
    enframe(name = NULL) %>%
    separate(value, c("LOCUS", "accn", "length", "bp", "biomol", "topology", "gbdiv", "modified"), "\\s+")
}

read_gbk_acc <- function(path, accession)
{
  lines <- read_lines(path)
  locus <- grep("^LOCUS", lines)
  version <- grep("^VERSION", lines)
  entry <- grep(accession, lines[version])[1]
  lines[locus[entry]:grep("^//", lines)[entry]]
}

positions <- function(ref)
{
  idx <- 0; pos <- c(); 
  for (chr in ref) { idx <- idx + (chr != "-"); pos <- c(pos, idx); }
  pos
}

call_ind <- function(msa)
{
  ref <- head(msa, 1)
  alt <- tail(msa, -1)
  
  mat <- (msa == "-") + 0
  gap <- (apply(mat, 2, sum) > 0) + 0
  run <- (mat[, 1:(ncol(mat)-1)] == mat[, 2:ncol(mat)]) %>% apply(2, all) %>% c(., last(.))
  
  paste(2 * gap + run, collapse = "") %>% 
    str_locate_all("(3+2?|2)") %>%
    as.data.frame() %>%
    apply(1, function(ele) {
      start <- ele["start"]
      end <- ele["end"]
      rng <- start:end
      data.frame(
        idx = 1 + 1:nrow(alt),
        REF = paste(ref[rng], collapse = ""), 
        ALT = apply(as.matrix(alt[, rng]), 1, paste, collapse = ""), 
        POS = start,
        stringsAsFactors = F
      )
    }) %>%
    bind_rows() %>%
    filter(REF != ALT & (str_detect(REF, "-") | str_detect(ALT, "-"))) %>%
    mutate_at("ALT", str_remove_all, "-") %>%
    mutate_at("REF", str_remove_all, "-") %>%
    mutate(len = nchar(ALT) - nchar(REF)) %>%
    mutate(type = c("ins", "del")[(len < 0) + 1])
}

call_snp <- function(msa)
{
  n <- 1
  lst <- list()
  for (pos in seq_along(msa[1, ]))
    for (idx in 2:nrow(msa))
      if (msa[1, pos] != "-" & msa[idx, pos] != "-" & msa[idx, pos] != msa[1, pos])
        lst[[n <- n + 1]] <- 
          data.frame(
            idx = idx, REF = msa[1, pos], ALT = msa[idx, pos], POS = pos, len = 1,
            stringsAsFactors = F
          ) %>% 
          mutate(type = ifelse(str_detect("AG GA CT TC", paste0(REF, ALT)), "trs", "trv"))
  bind_rows(lst)
}

overlevels <- function(ranges)
{
  data.frame(IRanges::findOverlaps(ranges)) %>% 
    setNames(c("x", "y")) %>% 
    filter(x <= y) %>% 
    pull("x") %>% 
    rle() %>% 
    .$lengths
}

plot_cds <- function(cds)
{
  lvl <- overlevels(cds)
  lvl <- factor(lvl, levels = rev(unique(lvl)))
  
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
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
}

plot_gmm <- function(vcf, input)
{
  types <- c("trv", "trs", "ins", "del")
  color <- setNames(c(input$color_trv, input$color_trs, input$color_ins, input$color_del), types)
  size <- setNames(c(input$size_trv, input$size_trv, input$size_ins, input$size_del), types)
  shape <- setNames(c(input$shape_trv, input$shape_trs, input$shape_ins, input$shape_del), types)
  
  ggplot(vcf, aes(POS, id)) +
  geom_point(aes(color = type, size = type, shape = type)) +
  scale_color_manual(values = color) +
  scale_size_manual(values = size) +
  scale_shape_manual(values = shape) +
  xlab("pos") +
  theme_minimal() +
  theme(legend.position = "bottom")
}
