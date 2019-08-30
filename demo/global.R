library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinycssloaders)
library(shinydashboard)
library(Biostrings)
library(tidyverse)
library(ggrepel)

roots = c(home = "..")

plot_trace <- function(msa)
{
  mat <- msa == "-"
  gap <- apply(mat, 2, any)
  run <- (mat[,1:(ncol(mat)-1)] == mat[,2:ncol(mat)]) %>% apply(2, all) %>% c(last(.))
  vec <- 2 * gap + run
  loc <- paste(vec, collapse = "") %>% str_locate_all("3+2?|2") %>% as.data.frame()
  
  data <- mutate(loc, xmin = start - 0.25, xmax = end + 0.25, mut = if_else(msa[1, start] == "-", "ins", "del"))
  
  df <- rbind(msa, gap = gap + 0, run = run + 0, vec = vec + 0) %>% as.data.frame()
  
  rownames_to_column(df, "id") %>%
    mutate_all(as.character) %>%
    gather(key, val, -id) %>% 
    mutate_at("id", factor, levels = rev(rownames(df))) %>%
    mutate_at("key", str_remove, "^V") %>%
    mutate_at("key", as.numeric) %>%
    ggplot(aes(key, id)) + 
    geom_text(aes(label = val)) + 
    geom_rect(
      data = data, ymin = -Inf, ymax = Inf, alpha = 0.25,
      aes(x = NULL, y = NULL, xmin = xmin, xmax = xmax, fill = mut)
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid = element_blank())
}
  
