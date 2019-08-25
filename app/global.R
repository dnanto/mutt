library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinycssloaders)
library(shinydashboard)
library(Biostrings)
library(tidyverse)
library(ggrepel)

roots = c(home = "..")

ext <- c("eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg", "wmf")
units <- c("in", "cm", "mm")

types <- c("trv", "trs", "ins", "del")

lineend <- c("round", "butt", "square")
linejoin <- c("round", "mitre", "bevel")
arrow_units <- c("cm", "mm", "inches")
arrow_ends <- c("last", "first", "both")
arrow_type <- c("open", "closed")

read_gbk_loc <- function(path)
{
  read_lines(path) %>%
    .[grep("^LOCUS", .)] %>%
    enframe(name = NULL) %>%
    separate(
      value, 
      c("LOCUS", "accn", "length", "bp", "biomol", "topology", "gbdiv", "modified"), 
      "\\s+"
    )
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

call_snp <- function(mat)
{
  ref <- mat[1, ]
  pos <- seq_along(ref)
  lapply(2:nrow(mat), function(idx) {
    row <- mat[idx, ]
    data.frame(POS = pos[ref != '-' & row != '-' & ref != row]) %>% 
      mutate(idx = idx, ALT = row[POS])
  }) %>% 
    bind_rows() %>%
    mutate(REF = ref[POS], len = 1) %>%
    mutate(type = ifelse(str_detect("AG GA CT TC", paste0(REF, ALT)), "trs", "trv"))
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
