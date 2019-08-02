#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)
file <- if (args[1] == "-") file("stdin") else args[1]
cols <- args[2:length(args)]

read_tsv(file("stdin"), col_types = cols(.default = "c")) %>%
  mutate_at(cols, lubridate::parse_date_time, orders = c("dbY", "Ymd", "bY", "Y"), quiet = T) %>%
  mutate_at(cols, strftime, format = "%Y-%m-%d") %>%
  write.table("", row.names = F, quote = F, sep = "\t")
