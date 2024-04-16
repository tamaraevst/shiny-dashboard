## app.R ##
library(shiny)
library(shinydashboard)
library(here)
library(DT)

source(here::here("src/data.R"))
source(here::here("src/pca_genes_panel.R"))
source(here::here("src/pca_samples_panel.R"))
source(here::here("src/enrichment_panel.R"))
source(here::here("src/enrichment.R"))
source(here::here("src/pathview_panel.R"))
source(here::here("src/ui.R"))

# Main function for the Shiny server

server <- function(input, output, session) {
  # myData <- reactive({
  #   inFile <- input$file1
  #   if (is.null(inFile)) return(NULL)
  #   data <- read.csv(inFile$datapath, header = TRUE)
  #   data
  # })

  output$contents <- DT::renderDataTable({
    DT::datatable(expression.matrix.freeze)
  })

  matrix <- memo_do_pca(expression.matrix.freeze)

  output$summary.pca <- DT::renderDataTable({
    DT::datatable(plot_pca_summary(matrix))
  })

  pca.plot <- reactive({
    myplot <- plot_pca(
      matrix,
      pcas = input$list_pcas,
      n.abundant = input$pca.n.abundant.genes,
      show.labels = input$enable_labels,
      label.size = input$pca.labels.size
    )
    myplot
  })
  output[["pca"]] <- renderPlot(pca.plot())

  pca.variance.plot <- reactive({
    myplot <- plot_pca_variance(
      matrix
    )
    myplot
  })
  output[["variance.pca"]] <- renderPlot(pca.variance.plot())

  # output$States_List <- renderUI({
  #   list_data <- colnames(meta)
  #   selectInput("group", "Select Meta Data Group", choices = list_data)
  # })

  # output$Enable_labels <- renderUI({
  #   list_data <- c(TRUE, FALSE)
  #   selectInput("labels", "Select Meta Data Group", choices = list_data)
  # })

  pca.samples.plot <- reactive({
    myplot <- plot_pca_samples(
      expression.matrix = expression.matrix.freeze,
      metadata = meta,
      pcas = input$list.pcas.samples,
      annotation.id = input$States_List,
      n.abundant = input$pca.n.abundant,
      show.labels = input$Enable_labels
    )
    myplot
  })
  output[["pca.samples"]] <- renderPlot(pca.samples.plot())

  compute.pca <- do_pca_reverse(expression.matrix.freeze)
  mytable <- plot_pca_summary(compute.pca)
  output$summary.pca.samples <- DT::renderDataTable({
    DT::datatable(mytable)
  })

  pca.samples.variance.plot <- reactive({
    myplot <- plot_pca_variance(
      compute.pca
    )
    myplot
  })
  output[["variance.pca.samples"]] <- renderPlot(pca.samples.variance.plot())

  # enrichment.plot <- reactive({
  #   myplot <- plot_enrichment(
  #     expression.matrix = expression.matrix.freeze,
  #     pca = input$list.pcas.enrichment
  #   )
  #   myplot
  # })
  # output[["enrichment"]] <- renderPlot(enrichment.plot())

  fgseaRes <- reactive({
    frame <- memo_do_enrichment(expression.matrix = expression.matrix.freeze, pca = input$list.pcas.enrichment)
    frame[frame$padj < 0.05, ]
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
  # scale_color_manual(
  # breaks = c('pval = 0.05'),
  # values = c('pval = 0.05' = 'purple')) +

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
      expression.matrix = expression.matrix.freeze,
      pcas = input$list.pcas.enrichment
    )
    myplot
  })
  output[["upgenes"]] <- renderPlot(enrichment.up.plot())

  enrichment.down.plot <- reactive({
    myplot <- bar_enrichment_down(
      expression.matrix = expression.matrix.freeze,
      pcas = input$list.pcas.enrichment
    )
    myplot
  })
  output[["downgenes"]] <- renderPlot(enrichment.down.plot())

  enrichment.table <- reactive({
    mytable <- table_enrichment(
      expression.matrix = expression.matrix.freeze,
      pcas = input$list.pcas.enrichment
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
        expression.matrix = expression.matrix.freeze,
        genes.list = get_genes_from_column(enrichment.table, s, 5)
      )
      myplot
    }
  })
  output[["heatmap"]] <- renderPlot(heatmap.plot())

  #  output[["entable.pathview"]] <- DT::renderDataTable(enrichment.table(),
  #     server = FALSE, rownames = FALSE, columnDefs = list(list(visible = FALSE, targets = 4))
  #     )

  observe(updateSelectInput(session, "condition1.pathview", choices = unique(meta[, input$metaname.pathview]), selected = unique(meta[, input$metaname.pathview])[1]))
  observe(updateSelectInput(session, "condition2.pathview", choices = unique(meta[, input$metaname.pathview]), selected = unique(meta[, input$metaname.pathview])[2]))
  observe(updateSelectInput(session, "control.pathview", choices = unique(meta[, input$metaname.pathview])))

  pv.out <- reactive({
    s <- input$entable_rows_selected
    id <- sub(".*?mmu(.*?)_.*", "\\1", enrichment.table()[s, 1])
    id.name <- enrichment.table()[s, 1]
    do_pathview(id, id.name, meta, input$condition1.pathview, input$condition2.pathview, input$metaname.pathview)
  })

  observeEvent(
    {
      input$condition1.pathview
      input$condition2.pathview
      input$entable_rows_selected
    },
    output$pathview <- renderImage(
      {
        s <- input$entable_rows_selected
        id <- sub(".*?mmu(.*?)_.*", "\\1", enrichment.table()[s, 1])

        # A temp file to save the output.
        # This file will be removed later by renderImage
        outfile <- paste("mmu", as.character(id), ".pathview.png", sep = "")

        # Return a list containing the filename
        list(
          src = outfile,
          contentType = "image/png"
          # width = 800,
          # height = 1000
        )
      },
      deleteFile = FALSE
    )
  )

  output$info <- renderPrint({
    req(input$plot.click)
    x <- round(input$plot.click$x, 3)
    y <- round(input$plot.click$y, 3)
    cat("[", x, ", ", y, "] \n", sep = "")

    data <- get_data_from_the_click(pv.out(), x, y)

    if (length(data[1]) != 0) {
      cat("You have clicked on the node with KEGG name", data[1], "and label", data[2])
    } else {
      cat("This node is not supported here.")
    }
  })

  pathview.density <- reactive({
    myplot <- plot_density(
      pv.out()$plot.data.gene
    )

    myplot
  })

  observeEvent(
    {
      input$condition1.pathview
      input$condition2.pathview
      input$entable_rows_selected
    },
    output$pathview.density <- renderPlot(pathview.density())
  )

  bar.plot.foldchange <- reactive({
    req(input$plot.click)
    x <- round(input$plot.click$x, 3)
    y <- round(input$plot.click$y, 3)
    data <- get_data_from_the_click(pv.out(), x, y)

    plot <- plot_fc_of_gene(
      gene.id = data[1],
      metadata = meta,
      control = input$control.pathview,
      annotation.id = input$metaname.pathview
    )

    plot
  })
  output$bar.plot.foldchange <- renderPlot(bar.plot.foldchange())
}

start.time <- Sys.time()
analysis.erichment <- memo_do_enrichment(expression.matrix.freeze, c("PCA1"), length(rownames(expression.matrix.freeze)))
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

shinyApp(ui, server)
