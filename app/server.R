shinyServer(function(input, output, session) {

	## msa

	#### get the path to the msa file
	shinyFileChoose(input, "import_msa", roots = roots)
	path_msa <- eventReactive(input$import_msa, pull(parseFilePaths(roots, input$import_msa), datapath))
	output$path_msa <- eventReactive(path_msa(), path_msa())

	#### load the multiple alignment
	msa <- eventReactive(path_msa(), {
		req(path <- path_msa())
		ape::read.dna(path, format = "fasta", as.character = T, as.matrix = T)
	})

	#### calc alignment <-> ref position mapping
	ref2aln <- eventReactive(
		msa(),
		pos_map(msa()[1, ]) %>%
			enframe(name = "aln", value = "ref") %>%
			distinct(ref, .keep_all = T) %>%
			pull(aln)
	)

	#### output the taxa table
	output$taxa <- DT::renderDT(
		DT::datatable(
			parse_header(rownames(msa())),
			rownames = FALSE,
			style = "bootstrap",
			class = "table-bordered table-striped table-hover responsive",
			filter = list(position = "top")
		)
	)

	#### update slider range based on msa length
	observeEvent(
		msa(),
		updateSliderInput(session, "range", max = ncol(msa()), value = c(1, 0.15 * ncol(msa())))
	)

	## gbk

	#### get the path to the genbank file
	shinyFileChoose(input, "import_gbk", roots = roots)
	path_gbk <- eventReactive(input$import_gbk, pull(parseFilePaths(roots, input$import_gbk), datapath))
	output$path_gbk <- eventReactive(path_gbk(), path_gbk())

	#### update the accession choices and clear the product choices
	observeEvent(path_gbk(), {
		req(path <- path_gbk())
		loci <- pull(read_gbk_loc(path), accn)
		updateSelectizeInput(session, "accession", choices = loci, selected = NA)
		updateSelectizeInput(session, "product", choices = c(""), selected = NA)
	})

	#### load the CDS objects
	cds <- eventReactive(input$accession, {
		req(path <- path_gbk(), pos <- ref2aln())
		genbankr::readGenBank(text = read_gbk_acc(path, input$accession), ret.seq = F) %>%
			genbankr::cds() %>%
			as.data.frame() %>%
			mutate(astart = pos[start], aend = pos[end]) %>%
			select(-translation)
	})

	#### output the CDS table
	output$cds <- DT::renderDT(
		DT::datatable(
			cds(),
			rownames = FALSE,
			style = "bootstrap",
			class = "table-bordered table-striped table-hover responsive",
			filter = list(position = "top"),
			options = list(scrollX = TRUE)
		)
	)

	#### update the product choices
	observeEvent(cds(), updateSelectizeInput(session, "product", choices = cds()$product, selected = NA))

	#### update the slider based on the selected product
	observeEvent(input$product, {
		cds <- cds()
		idx <- which(cds$product == input$product)
		updateSliderInput(session, "range", value = c(min(cds$astart[idx]), max(cds$aend[idx])))
	})

	## map

	#### click on the map tab on run button click
	observeEvent(input$run, updateTabsetPanel(session, "tabset_panel", selected = "tab_map"))

	#### calculate the variation calls
	calls <- eventReactive(input$run, {
		req(msa <- msa())

		lab <- parse_header(rownames(msa)) %>% pull(id)

		snp <-
			call_snp(msa) %>%
			mutate_at("call", factor, levels = types) %>%
			mutate(id = factor(lab[idx], levels = rev(lab)))

		ind <-
			call_ind(msa) %>%
			mutate_at("call", factor, levels = types) %>%
			mutate(id = factor(lab[idx], levels = rev(lab)))

		list(lab = lab, snp = snp, ind = ind)
	})

	#### plot

	plot_cow <- function()
	{
		plt <- plot_map(calls(), input)
		if (nchar(input$accession) > 0)
		{
			plt <- cowplot::plot_grid(
				plot_cds(cds(), input), plt, ncol = 1, align = "v",
				rel_heights = c(input$rel_height_cds, input$rel_height_map)
			)
		}
		plt
	}

	observe(output$map <- renderPlot(plot_cow(), height = input$height))

	## calls

	### output the calls table
	output$calls <- DT::renderDT({
		calls <- calls()
		DT::datatable(
			arrange(select(bind_rows(calls$snp, calls$ind), id, everything(), -idx), id, pos),
			rownames = FALSE,
			style = "bootstrap",
			class = "table-bordered table-striped table-hover responsive",
			filter = list(position = "top")
		)
	})

	output$export <- downloadHandler(
		filename = function() paste(calls()$lab[1], input$ext, sep="."),
		content = function(file)
			ggsave(
				file, plot = plot_cow(),
				width = input$exp_width, height = input$exp_height,
				units = input$units, scale = input$scale, dpi = input$dpi
			)
	)

})
