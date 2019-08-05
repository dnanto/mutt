library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(tidyverse)

roots = c(home = "..")
types <- c("trv", "trs", "ins", "del")
color <- setNames(c("magenta", "cyan", "black", "black"), types)
size <- setNames(c(10, 10, 2, 2), types)
shape <- setNames(c(124, 124, 6, 2), types)

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
  fields <- str_split(lines[grep("^#CHROM", lines)], "\t", simplify = T) %>% str_remove("^#")
  if (!is_empty(lines))
    read_tsv(lines, col_names = fields, comment = "#") %>%
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
    mutate_if(is.character, as.factor)
}

overlevels <- function(ranges)
{
  data.frame(findOverlaps(ranges)) %>% 
    setNames(c("x", "y")) %>% 
    filter(x <= y) %>% 
    pull("x") %>% 
    rle() %>% 
    .$lengths
}
