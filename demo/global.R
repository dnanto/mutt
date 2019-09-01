library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(Biostrings)
library(tidyverse)
library(magrittr)

roots = c(home = "..")

plot_trace <- function(msa)
{
  mat <- msa == "-"
  gap <- apply(mat, 2, any) - apply(mat, 2, all)
  run <- (mat[,1:(ncol(mat)-1)] == mat[,2:ncol(mat)]) %>% apply(2, all) %>% c(last(.))
  vec <- 2 * gap + run
  loc <- paste(vec, collapse = "") %>% str_locate_all("3+2?|2") %>% as.data.frame()
  
  df <- rbind(msa, gap = gap + 0, run = run + 0, vec = vec + 0) %>% as.data.frame()
  
  calls <- NULL
  if (nrow(loc) > 0)
      calls <-
      apply(loc, 1, function(ele) {
        rng <- ele["start"]:ele["end"]
        data.frame(
          id = tail(rownames(msa), -1),
          REF = paste(msa[1, rng], collapse = ""), 
          ALT = apply(as.matrix(msa[2:nrow(msa), rng]), 1, paste, collapse = ""), 
          POS = ele["start"],
          stringsAsFactors = F
        )
      }) %>%
      bind_rows() %>%
      filter(REF != ALT & (str_detect(REF, "-") | str_detect(ALT, "-"))) %>%
      mutate_at("id", factor, levels = rev(rownames(df))) %>%
      mutate(
        xmin = POS - 0.45, xmax = POS + nchar(ALT) - 1 + 0.45, 
        ymin = as.numeric(id) - 0.45, ymax = as.numeric(id) + 0.45,
        mut = factor(if_else(str_detect(REF, "-"), "ins", "del"), levels = c("ins", "del"))
      )
  
  boxes <- 
    mutate(loc, xmin = start - 0.25, xmax = end + 0.25, mut = if_else(msa[1, start] == "-", "ins", "del")) %>%
    mutate_at("mut", factor, levels = c("ins", "del"))
  
  plot <-
    rownames_to_column(df, "id") %>%
      mutate_all(as.character) %>%
      gather(pos, val, -id) %>% 
      mutate_at("id", factor, levels = rev(rownames(df))) %>%
      mutate_at("pos", str_remove, "^V") %>%
      mutate_at("pos", as.numeric) %>%
      ggplot(aes(pos, id)) + 
      geom_text(aes(label = val)) + 
      geom_rect(
        data = boxes, ymin = -Inf, ymax = Inf, alpha = 0.25,
        aes(x = NULL, y = NULL, xmin = xmin, xmax = xmax, fill = mut)
      ) +
      { 
        if (!is.null(calls)) 
          geom_rect(
            data = calls, alpha = 0.25,
            aes(x = NULL, y = NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = mut)
          ) 
      } +
      scale_fill_manual(values = c("ins" = "blue", "del" = "red")) +
      theme_minimal() +
      theme(legend.position = "bottom", panel.grid = element_blank())
  
  list(plot = plot, calls = calls)
}
  
