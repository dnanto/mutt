library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(shinydashboard)
library(colourpicker)
library(tidyverse)
library(ggrepel)

roots <- c(home = "..")

types <- c("trv", "trs", "ins", "del")

ext <- c("eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg", "wmf")
units <- c("in", "cm", "mm")
lineend <- c("round", "butt", "square")
linejoin <- c("round", "mitre", "bevel")
arrow_units <- c("cm", "mm", "inches")
arrow_ends <- c("last", "first", "both")
arrow_type <- c("open", "closed")

parse_header <- function(string)
{
	str_split_fixed(string, " ", 2) %>% as.data.frame() %>% set_names(c("id", "definition"))
}

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

pos_map <- function(ref)
{
	j <- 0; pos <- rep(0, length(ref));
	for (i in seq_along(ref)) pos[i] <- j <- j + (ref[i] != "-")
	pos
}

call_snp <- function(msa)
{
	ref <- msa[1, ]
	pos <- seq_along(ref)
	bind_rows(lapply(2:nrow(msa), function(idx) {
		alt <- msa[idx, ]
		pos <- pos[ref != '-' & alt != '-' & ref != alt]
		data.frame(idx = idx, pos = pos, ref = ref[pos], alt = alt[pos], stringsAsFactors = F)
	})) %>%
		mutate_at(c("ref", "alt"), toupper) %>%
		mutate(call = case_when(
			ref == "A" & alt == "C" ~ "trv", ref == "C" & alt == "A" ~ "trv",
			ref == "A" & alt == "G" ~ "trs", ref == "G" & alt == "A" ~ "trs",
			ref == "A" & alt == "T" ~ "trv", ref == "T" & alt == "A" ~ "trv",
			ref == "C" & alt == "G" ~ "trv", ref == "G" & alt == "C" ~ "trv",
			ref == "C" & alt == "T" ~ "trs", ref == "T" & alt == "C" ~ "trs",
			ref == "G" & alt == "T" ~ "trv", ref == "T" & alt == "G" ~ "trv"
		)) %>%
		filter(!is.na(call))
}

call_ind <- function(msa)
{
	loc <-
		sweep(tail(msa, -1) == "-", 2, msa[1,] == "-") %>% abs() %>%
		apply(1, paste, collapse = "") %>% str_locate_all("1+")

	bind_rows(lapply(seq_along(loc), function(idx) {
		ele	<- loc[[idx]]
		bind_rows(apply(as.data.frame(ele), 1, function(obj) {
			idx <- idx + 1; start <- obj["start"]; end <- obj["end"]; rng <- start:end;
			data.frame(
				idx = idx, pos = start,
				ref = paste(msa[1, rng], collapse = ""), alt = paste(msa[idx, rng], collapse = ""),
				stringsAsFactors = F, row.names = NULL
			)
		}))
	})) %>%
		filter(ref != alt) %>%
		mutate(call = factor(if_else(str_detect(ref, "-"), "ins", "del"), c("ins", "del")))
}

call_ind_overlap <- function(msa)
{
	mat <- msa == "-"
	gap <- apply(mat, 2, any) - apply(mat, 2, all)
	run <- (mat[, 1:(ncol(mat)-1)] == mat[, 2:ncol(mat)]) %>% apply(2, all) %>% c(last(.))
	loc <- paste(2 * gap + run, collapse = "") %>% str_locate_all("3+2?|2") %>% as.data.frame()

	bind_rows(apply(loc, 1, function(ele) {
		rng <- ele["start"]:ele["end"]
		data.frame(
			idx = 2:nrow(msa), pos = ele["start"],
			ref = paste(msa[1, rng], collapse = ""),
			alt = apply(as.matrix(msa[2:nrow(msa), rng]), 1, paste, collapse = ""),
			stringsAsFactors = F, row.names = NULL
		)
	})) %>%
		filter(ref != alt) %>%
		mutate(call = factor(if_else(str_detect(ref, "-"), "ins", "del"), c("ins", "del")))
}

overlevels <- function(ranges)
{
	IRanges::findOverlaps(ranges) %>%
		as.data.frame() %>%
		setNames(c("x", "y")) %>%
		filter(x <= y) %>%
		pull("x") %>%
		rle() %>%
		.$lengths
}

plot_cds <- function(cds, input)
{
	cds <-
		mutate(cds, x = ifelse(strand == '-', end, start), xend = ifelse(strand == "+", end, start)) %>%
		filter(input$range[1] <= x & x <= input$range[2]) %>%
		mutate(xend = pmin(xend, input$range[2]))

	lvl <-
		IRanges::IRanges(cds$start, cds$end) %>%
		overlevels() %>%
		factor(., levels = rev(unique(.)))

	ggplot(cds, aes(x = x, xend = xend, y = lvl, yend = lvl)) +
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
		xlim(input$range[1], input$range[2]) +
		xlab("") +
		ylab("product") +
		ggtitle(paste("Coding [", input$accession, "]")) +
		theme_minimal() +
		theme(
			legend.position = "none",
			title = element_text(size = input$title_size),
			axis.text.x = element_blank(),
			axis.text.y = element_blank()
		)

}

plot_map <- function(calls, input)
{
	lab <- calls$lab
	ind <- filter(calls$ind, call %in% input$types & input$range[1] <= pos & pos <= input$range[2])
	snp <- filter(calls$snp, call %in% input$types & input$range[1] <= pos & pos <= input$range[2])

	trv <- filter(snp, call == "trv")
	trs <- filter(snp, call == "trs")
	ins <- filter(ind, call == "ins")
	del <- filter(ind, call == "del")

	y_trs <- position_nudge(y = input$y_trs)
	y_trv <- position_nudge(y = input$y_trv)
	y_ins <- position_nudge(y = input$y_ins)
	y_del <- position_nudge(y = input$y_del)

	color <- setNames(c(input$color_trv, input$color_trs, input$color_ins, input$color_del), types)
	shape <- setNames(c(input$shape_trv, input$shape_trs, input$shape_ins, input$shape_del), types)
	size <- setNames(c(input$size_trv, input$size_trv, input$size_ins, input$size_del), types)

	ggplot() +
		geom_point(data = trs, aes(pos, id, color = call, size = call, shape = call), position = y_trv) +
		geom_point(data = trv, aes(pos, id, color = call, size = call, shape = call), position = y_trs) +
		geom_point(data = ins, aes(pos, id, color = call, size = call, shape = call), position = y_ins) +
		geom_point(data = del, aes(pos, id, color = call, size = call, shape = call), position = y_del) +
		scale_color_manual(values = color) +
		scale_size_manual(values = size) +
		scale_shape_manual(values = shape) +
		xlim(input$range[1], input$range[2]) +
		ggtitle(paste("Variation [", lab[1], "]")) +
		theme_minimal() +
		theme(
			legend.position = "bottom",
			title = element_text(size = input$title_size)
		)
}
