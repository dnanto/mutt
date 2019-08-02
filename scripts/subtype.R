#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)
file <- if (length(args) == 0 || args[1] == "-") file("stdin") else args[1]

read_tsv(file, col_names = c("AccessionVersion", "SubType", "SubName")) %>%
bind_cols(
  .,
  bind_rows(
    apply(., 1, function(row) {
      key <- str_split(row["SubType"], "\\|")[[1]]
      val <- str_split(row["SubName"], "\\|")[[1]]
      data.frame(as.list(setNames(val, key)), stringsAsFactors = F)
    })
  )
) %>%
select(-SubType, -SubName) %>%
write.table("", row.names = F, quote = F, sep = "\t")
