dashboardPage(
	dashboardHeader(title = "mutt"),
	dashboardSidebar(
		actionButton("run", "Run", width = "93%"),
		sliderInput("range", "range", min = 1, max = 1, value = c(1, 1), width = "100%"),

		fluidRow(
			column(4, numericInput("height", "height", 800, step = 100), style = "padding-right: 0px;"),
			column(4, numericInput("rel_height_cds", "rel_height_cds", 1, step = 1, min = 1), style = "padding: 0px;"),
			column(4, numericInput("rel_height_map", "rel_height_map", 4, step = 1, min = 1), style = "padding-left: 0px;")
		),
		checkboxGroupInput("types", "types", choices = types, selected = types, inline = T, width = "100%"),
		fluidRow(
			column(6, selectizeInput("accession", "accession", NULL), style = "padding-right: 0px;"),
			column(6, selectizeInput("product", "product", NULL), style = "padding-left: 0px;")
		),
		fluidRow(
			column(3, numericInput("title_size", "title_size", 10, min = 0, width = "100%"), style = "padding-right: 0px;"),
			column(3, numericInput("text_size_cds", "text_size_cds", 2, min = 0, step = 0.5), style = "padding: 0px;"),
			column(3, numericInput("line_size_cds", "line_size_cds", 0.25, min = 0, step = 0.25), style = "padding: 0px;"),
			column(3, checkboxInput("label", "label", value = T), style = "padding-left: 0px;")
		),
		fluidRow(
			column(4, numericInput("arrow_length", "arrow_length", 1, min = 0, step = 1), style = "padding-right: 0px;"),
			column(4, selectInput("arrow_units", "arrow_units", arrow_units, arrow_units[2]), style = "padding: 0px;"),
			column(4, selectInput("arrow_type", "arrow_type", arrow_type, arrow_type[2]), style = "padding-left: 0px;")
		),
		fluidRow(
			column(4, selectInput("lineend", "line-end", lineend, lineend[2]), style = "padding-right: 0px;"),
			column(4, selectInput("linejoin", "line-join", linejoin, linejoin[2])), style = "padding: 0px;",
			column(4, numericInput("segment_size", "segment_size", 2, min = 0, step = 1), style = "padding-left: 0px;")
		),
		fluidRow(
			column(3, colourInput("color_trv", "color_trv", value = "magenta", allowTransparent = T), style = "padding-right: 0px;"),
			column(3, numericInput("shape_trv", "shape_trv", 124, min = 0, max = 127), style = "padding: 0px;"),
			column(3, numericInput("size_trv", "size_trv", 10, min = 0), style = "padding: 0px;"),
			column(3, numericInput("y_trv", "y_trv", 0, step = 0.05), style = "padding-left: 0px;")
		),
		fluidRow(
			column(3, colourInput("color_trs", "color_trs", value = "cyan", allowTransparent = T), style = "padding-right: 0px;"),
			column(3, numericInput("shape_trs", "shape_trs", 124, min = 0, max = 127), style = "padding: 0px;"),
			column(3, numericInput("size_trs", "size_trs", 10, min = 0), style = "padding: 0px;"),
			column(3, numericInput("y_trs", "y_trs", 0, step = 0.05), style = "padding-left: 0px;")
		),
		fluidRow(
			column(3, colourInput("color_ins", "color_ins", value = "black", allowTransparent = T), style = "padding-right: 0px;"),
			column(3, numericInput("shape_ins", "shape_ins", 6, min = 0, max = 127), style = "padding: 0px;"),
			column(3, numericInput("size_ins", "size_ins", 4, min = 0), style = "padding: 0px;"),
			column(3, numericInput("y_ins", "y_ins", 0, step = 0.05), style = "padding-left: 0px;")
		),
		fluidRow(
			column(3, colourInput("color_del", "color_del", value = "black", allowTransparent = T), style = "padding-right: 0px;"),
			column(3, numericInput("shape_del", "shape_del", 2, min = 0, max = 127), style = "padding: 0px;"),
			column(3, numericInput("size_del", "size_del", 4, min = 0), style = "padding: 0px;"),
			column(3, numericInput("y_del", "y_del", 0, step = 0.05), style = "padding-left: 0px;")
		),
		width = 450
	),
	dashboardBody(
		tabsetPanel(
			id = "tabset_panel",
			tabPanel(
				"Taxa", value = "tab_taxa",
				box(
					width = 12,
					wellPanel(
						shinyFilesButton(
							"import_msa", "Import MSA", "Please select a multiple sequence alignment FASTA file...",
							multiple = F, icon = icon("file-import")
						),
						verbatimTextOutput("path_msa", placeholder = T)
					),
					DT::DTOutput("taxa")
				)
			),
			tabPanel(
				"CDS", value = "tab_cds",
				box(
					width = 12,
					wellPanel(
						shinyFilesButton(
							"import_gbk", "Import GenBank", "Please select a GenBank file...",
							multiple = F, icon = icon("file-import")
						),
						verbatimTextOutput("path_gbk", placeholder = T)
					),
					DT::DTOutput("cds")
				)
			),
			tabPanel(
				"Map", value = "tab_map", withSpinner(plotOutput("map"), type = 8)
			),
			tabPanel(
				"Calls", value = "tab_calls",
				box(width = 12, withSpinner(DT::DTOutput("calls"), type = 8))
			),
			tabPanel(
				"Export", value = "tab_export",
				box(
					width = 12,
					wellPanel(
						fluidRow(
							column(6, selectInput("ext", "format", choices = ext, selected = "pdf")),
							column(6, selectInput("units", "units", choices = units, selected = "in"))
						),
						fluidRow(
							column(6, numericInput("exp_width", "width", 10, min = 0)),
							column(6, numericInput("exp_height", "height", 15, min = 0))
						),
						fluidRow(
							column(6, numericInput("scale", "scale", 1, min = 0.1, step = 0.1)),
							column(6, numericInput("dpi", "dpi", 300, min = 72, step = 50))
						)
					),
					wellPanel(downloadButton("export", "Export"))
				)
			)
		)
	),
	skin = "green"
)
