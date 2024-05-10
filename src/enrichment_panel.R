## enrichment_panel.R ##

library(ggplot2)
library(dplyr)
library(reshape)

source(here::here("src/utils.R"))

Enrichment.ui <- function(id) {
  tabItem(
    tabName = "Enrichment",
    fluidRow(
      box(plotOutput(NS(id, "enrichment"), height = 400, hover = hoverOpts(NS(id, "plot.hover"))),
        verbatimTextOutput(NS(id, "hover.info")),
        width = 8, title = "Scatter of enrichmnet analysis results", footer = "Above we display the overall distribution of adjusted p-values across all enrichment pathways that are filtered out by the set p-value threshold on the slider. Here enrichmnet pathways are labelled by their index on the x-axis and their corresponding adjusted p-values on the y-axis. You may hover over the plotted points to see what enrichment pathway they represent, the names will be displayed as they appear in the pathway file you have provided. When this page is loaded, it performs Gene Set Enrichment Analysis using a ranked list of genes from your expression data. The genes are ranked by their PCA loadings. The first instance of the page diplays the results of the enrichmnent analysis using the default pricincipal component (the first one) to rank the genes."
      ),
      box(selectInput(inputId = NS(id, "list.pcas.enrichment"), label = "Choose your principal component from the list", choices = c("PC1", "PC2", "PC3"), multiple = FALSE, selected = c("PC1")), width = 4),
      box(sliderInput(NS(id, "enrichment.p.threshold"),
        label = "Select p-value threshold",
        min = 0.001, value = 0.05, max = 0.05, step = 0.01, ticks = TRUE
      ), width = 3), # box(
      #   sliderInput("enrichment.n.abundant",
      #     label = "Select top number of genes ordered by PCA results",
      #     min = 10, value = length(rownames(expression.matrix)), max = length(rownames(expression.matrix)), step = 500, ticks = TRUE
      #   ),
      #   width = 4
      # ),
      box(plotOutput(NS(id, "upgenes"), height = 500), width = 6, title = "Bar plot of top 10 enrichmnent pathways with positive enrichment score", footer = "A positive enrichment score means that a gene set belonging to a particular pathway is more enriched in the positively-regulated genes."),
      box(plotOutput(NS(id, "downgenes"), height = 500), width = 6, title = "Bar plot of top 10 enrichmnent pathways with negative enrichment score", footer = "A negative enrichment score corresponds to enrichment in the negatively regulated genes."),
      box(DT::dataTableOutput(NS(id, "entable"), fill = TRUE),
        width = 12, title = "Summary table of enrichment analysis", footer = "You can click on the rows of the table to extract more information on the pathways. In particular, a list and a heat-map of leading edge genes belogning to a chosen pathway will be displayed below."
      ),
      box(verbatimTextOutput(NS(id, "engenes")),
        width = 12, title = "Leading edge genes", footer = "Leading edge genes are the subset of genes found in the ranking at the maximal ES."
      ),
      box(plotOutput(NS(id, "heatmap")),
        width = 12, title = "Heatmap of leading edge genes"
      )
    )
  )
}

Enrichment.server <- function(id, expression.matrix, pathways.gmt.file) {
  moduleServer(
    id,
    function(input, output, session) {
      fgseaRes <- reactive({
        frame <- memo_do_enrichment(expression.matrix = expression.matrix, pca = input$list.pcas.enrichment, pathways.gmt.file = pathways.gmt.file)
        frame[frame$padj < input$enrichment.p.threshold, ]
      })


      output$enrichment <- renderPlot({
        ggplot(data = fgseaRes(), aes(x = as.numeric(row.names(fgseaRes())), y = -log10(padj))) +
          geom_point(color = "purple", size = 5) +
          labs(x = "Pathway index", y = "-log10(padj)") +
          geom_hline(yintercept = -log10(0.05)) +
          scale_shape_discrete(
            name = "Legend",
            labels = c("padj=0.05")
          ) +
          theme_minimal()
      })

      displayed_text <- reactive({
        req(input$plot.hover)
        hover <- input$plot.hover
        dist <- sqrt((hover$x - as.numeric(row.names(fgseaRes())))^2 + (hover$y + log10(fgseaRes()$padj))^2)

        if (min(dist) < 0.3) {
          fgseaRes()$pathway[which.min(dist)]
        } else {
          NULL
        }
      })

      output$hover.info <- renderPrint({
        req(displayed_text())

        cat("Name\n")
        displayed_text()
      })

      enrichment.up.plot <- reactive({
        myplot <- bar_enrichment_up(
          fgseaRes = fgseaRes(),
          pathways.gmt.file = pathways.gmt.file
        )
        myplot
      })
      output[["upgenes"]] <- renderPlot(enrichment.up.plot())

      enrichment.down.plot <- reactive({
        myplot <- bar_enrichment_down(
          fgseaRes = fgseaRes(),
          pathways.gmt.file = pathways.gmt.file
        )
        myplot
      })
      output[["downgenes"]] <- renderPlot(enrichment.down.plot())

      enrichment.table <- reactive({
        mytable <- table_enrichment(
          fgseaRes = fgseaRes(),
          pathways.gmt.file = pathways.gmt.file
        )
        mytable
      })
      output[["entable"]] <- DT::renderDataTable(enrichment.table(),
        server = FALSE, rownames = FALSE, extensions = "Buttons", selection = "single",
        options = list(autoWidth = TRUE, buttons = c("csv", "excel"), dom = "Bfrtip", columnDefs = list(list(visible = FALSE, targets = 4)))
      )

      s <- reactiveVal(-1)
      observe({
        s <- input$entable_rows_selected
      })

      output$engenes <- renderPrint({
        s <- input$entable_rows_selected
        if (length(s)) {
          cat("The following 'leading edge' genes fall under the selected enrichment category:\n\n")
          selected.genes <- get_genes_from_column(enrichment.table, s, 5)
          # enrichment.table()[s, 5]
          selected.genes <- as.vector(strsplit(selected.genes, ",")[[1]])

          find.labels <- expression.matrix.symbols[expression.matrix.symbols$ensembl %in% selected.genes, ]
          labels <- find.labels$NAME
          cat(labels)
        } else {
          cat("Click on rows to print the 'leading edge' genes belonging to the selected category\n\n")
        }
      })

      heatmap.plot <- reactive({
        s <- input$entable_rows_selected
        if (length(s)) {
          myplot <- heatmap_enrichment(
            expression.matrix = expression.matrix,
            genes.list = get_genes_from_column(enrichment.table, s, 5)
          )
          myplot
        }
      })
      output[["heatmap"]] <- renderPlot(heatmap.plot())

      rows.selected <- reactive({
        s <- input$entable_rows_selected
        enrichment.table()[s, ]
      })

      return(rows.selected)
    }
  )
}

# Function for the Manhattan plot enrichment analysis results. The user is prompted to select top abundant genes from the slider, these genes are then used for plotting.

plot_enrichment <- function(
    fgseaRes,
    dataset.organism = "mmusculus") {
  enrichment.plot <- ggplot(data = fgseaRes, aes(x = as.numeric(row.names(fgseaRes)), y = -log10(padj))) +
    geom_point(color = "purple") +
    labs(x = "Pathway index", y = "-log10(padj)") +
    geom_hline(yintercept = -log10(0.05)) +
    scale_shape_discrete(
      name = "Legend",
      labels = c("padj=0.05")
    )
  theme_minimal()

  enrichment.plot
}

# A dot plot visualisiting pathways with positive ESs. TODO: write one function for the dot plot that will show positive and negative ESs; need just one additional parameter to implement this.

bar_enrichment_up <- function(
    fgseaRes,
    dataset.organism = "mmusculus",
    pathways.gmt.file) {
  topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n = 10), ]

  genes.up <- topPathwaysUp[, c("pathway", "padj", "size", "ES")]
  genes.up <- genes.up[order(genes.up$ES), ]

  bar.plot.up <- ggplot(genes.up, aes(x = ES, y = pathway, size = size, colour = log10(padj))) +
    geom_point(aes(fill = log10(padj)), colour = "black", shape = 21) +
    viridis::scale_color_viridis(option = "plasma", breaks = scales::pretty_breaks(n = 5)) +
    viridis::scale_fill_viridis(option = "plasma") +
    theme_minimal() +
    labs(x = "Enrichment score", y = "Pathway") +
    scale_y_discrete(limits = genes.up$pathway)

  # rng <- range(genes.up$padj)

  # bar.plot.up <- ggplot(genes.up, aes(x = size, y = pathway, fill = padj)) +
  #   scale_fill_continuous(
  #     low = "purple", high = "blue", name = "padj",
  #     guide = guide_colorbar(reverse = TRUE)
  #   ) +
  #   scale_fill_gradient2(low = "blue", mid = "cyan", high = "purple",
  #                      midpoint = mean(rng),  # Same midpoint for both plots
  #                      breaks = seq(-100, 100, 4),  # Breaks in the scale bar
  #                      limits = c(floor(rng[1]), ceiling(rng[2]))) +
  #   geom_col() +
  #   labs(x = "Size", y = "Pathway") +
  #   scale_y_discrete(limits = genes.up$pathway)

  # bar.plot.up + coord_flip()

  return(bar.plot.up)
}

# A dot plot visualisiting pathways with negative ESs.

bar_enrichment_down <- function(
    fgseaRes,
    dataset.organism = "mmusculus",
    pathways.gmt.file) {
  topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n = 10), ]

  genes.down <- topPathwaysDown[, c("pathway", "padj", "size", "ES")]
  genes.down <- genes.down[order(genes.down$ES), ]

  bar.plot.down <- ggplot(genes.down, aes(x = ES, y = pathway, size = size, colour = log10(padj))) +
    geom_point(aes(fill = log10(padj)), colour = "black", shape = 21) +
    viridis::scale_color_viridis(option = "plasma", breaks = scales::pretty_breaks(n = 5)) +
    viridis::scale_fill_viridis(option = "plasma") +
    theme_minimal() +
    labs(x = "Enrichment score", y = "Pathway") +
    scale_y_discrete(limits = genes.down$pathway)

  return(bar.plot.down)
}

# Table with the summary from the enrichment analysis, based off the top abundant genes selected by the user

table_enrichment <- function(
    fgseaRes,
    dataset.organism = "mmusculus",
    pathways.gmt.file) {
  pathways <- fgseaRes$pathway
  pval <- fgseaRes$pval
  padj <- fgseaRes$padj
  size <- fgseaRes$size
  leadingEdge <- fgseaRes$leadingEdge

  enrichment.table <- data.frame(Pathways = pathways, pval = pval, padj = padj, Size = size, Genes = sapply(leadingEdge, paste, collapse = ","))

  return(enrichment.table)
}

# Heatmap of genes wiht normalised expression z-scores for a given pathway index

heatmap_enrichment <- function(
    expression.matrix,
    genes.list) {
  genes.list <- as.vector(strsplit(genes.list, ",")[[1]])

  expression.matrix.sorted <- as.data.frame(expression.matrix) %>% filter(row.names(expression.matrix) %in% genes.list)

  find.labels <- expression.matrix.symbols[expression.matrix.symbols$ensembl %in% rownames(expression.matrix.sorted), ]
  labels <- find.labels$NAME

  scaled.expression.matrix <- as.data.frame(scale(expression.matrix.sorted))

  data.genes <- cbind(genes = labels, data.frame(scaled.expression.matrix, row.names = NULL))

  expression.marix.melt <- melt(data.genes)

  enrichment.heatmap <- ggplot(expression.marix.melt, aes(genes, variable)) + # Create heatmap with ggplot2
    geom_tile(aes(fill = value)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    viridis::scale_color_viridis(discrete = FALSE) +
    viridis::scale_fill_viridis(discrete = FALSE) +
    labs(x = "", y = "", fill = "z-score of expression")

  return(enrichment.heatmap)
}
