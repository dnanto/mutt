library(shiny)
library(shinyFiles)
library(tidyverse)
library(lubridate)
library(ggtree)
library(ape)

roots = c(home = "..")
orders <- c("dbY", "Ymd", "bY", "Y")
objective <- c("correlation", "rsquared", "rms")
distmethod <- c("patristic", "nNodes", "Abouheif", "sumDD")
model <- c("poisson", "negbin", "strictgamma", "relaxedgamma", "mixedgamma")

parse_tip_date <- function(string, pattern = "_")
{
  ymd(lapply(str_split(string, pattern), tail, 1))
}
