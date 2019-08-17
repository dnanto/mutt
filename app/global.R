library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinycssloaders)
library(tidyverse)

roots = c(home = "..")

read_gbk_acc <- function(path, accession)
{
  lines <- read_lines(path)
  locus <- grep("^LOCUS", lines)
  version <- grep("^VERSION", lines)
  entry <- grep(accession, lines[version])[1]
  lines[locus[entry]:grep("^//", lines)[entry]]
}

read_vcf <- function(file)
{
  lines <- read_lines(file)
  fields <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  if (!is_empty(lines))
    read_tsv(lines, col_names = fields, col_types = "cicccccc", comment = "#") %>%
    mutate(id = basename(tools::file_path_sans_ext(file)))
}

calls <- function(files)
{
  vcf <-
    lapply(files, read_vcf) %>%
    bind_rows() %>%
    select(-CHROM, -ID, -QUAL, -FILTER, -INFO)
  
  ind <-
    filter(vcf, nchar(REF) > 1 | nchar(ALT) > 1) %>%
    mutate(type = c("ins", "del")[(nchar(REF) > 1) + 1], len = nchar(ALT) - nchar(REF))
  
  snp <-
    filter(vcf, nchar(REF) == 1 & nchar(ALT) == 1 & !(REF == "N" | ALT == "N")) %>%
    mutate(type = ifelse(str_detect("AG GA CT TC", paste0(REF, ALT)), "trs", "trv"))
  
  bind_rows(ind, snp) %>%
  select(id, type, POS, REF, ALT, len) %>%
  mutate_at("id", as.factor) %>%
  mutate_at("type", factor, levels = c("trv", "trs", "ins", "del"))
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
  
  gap <- 2 * apply(alt, 2, grepl, pattern = "-") %>% apply(2, any) %>% as.integer()
  run <- (alt[, 1:(ncol(alt)-1)] == alt[, 2:ncol(alt)]) %>% apply(2, all) %>% c(., last(.)) %>% as.integer()
  
  paste(gap + run, collapse = "") %>% 
    str_locate_all("(2|3+2?)") %>%
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
    mutate(type = c("ins", "del")[(nchar(REF) > 0) + 1], len = nchar(ALT) - nchar(REF))
}

call_snp <- function(msa)
{
  n <- 1
  lst <- list()
  for (pos in seq_along(msa[1, ]))
    for (idx in 2:nrow(msa))
      if (msa[1, pos] != "-" & msa[idx, pos] != "-" & msa[idx, pos] != msa[1, pos])
        lst[[n <- n + 1]] <- data.frame(
          idx = idx, REF = msa[1, pos], ALT = msa[idx, pos], POS = pos, 
          stringsAsFactors = F
        )
  bind_rows(lst) %>%
    mutate(type = ifelse(str_detect("AG GA CT TC", paste0(REF, ALT)), "trs", "trv"), len = 1)
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
