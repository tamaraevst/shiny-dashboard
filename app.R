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

  enrichment.plot <- reactive({
    myplot <- plot_enrichment(
      expression.matrix = expression.matrix.freeze,
      pcas = input$list.pcas.enrichment
    )
    myplot
  })
  output[["enrichment"]] <- renderPlot(enrichment.plot())

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

  observeEvent(input$showPathview.enrichment, {
    output$pathview.enrichment <- renderImage(
      {
        s <- input$entable_rows_selected
        if (length(s)) {
          pathway <- sub(".*?mmu(.*?)_.*", "\\1", enrichment.table()[s, 1])

          start.time <- Sys.time()

          do_pathview(pathway, enrichment.table()[s, 5])

          end.time <- Sys.time()
          time.taken <- end.time - start.time
          cat("Time taken for pathview: \n", time.taken) 

          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- paste("mmu", as.character(pathway), ".png", sep = "")

          # Return a list containing the filename
          list(
            src = outfile,
            contentType = "image/png",
            width = 700,
            height = 800,
            alt = "This is alternate text"
          )
        }
      },
      deleteFile = TRUE
    )
  })

  observeEvent(input$showPathview.enrichment, {
    output$pathview.enrichment2 <- renderImage(
      {
        s <- input$entable_rows_selected
        if (length(s)) {
          pathway <- sub(".*?mmu(.*?)_.*", "\\1", enrichment.table()[s, 1])

          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- paste("mmu", as.character(pathway), ".pathview.png", sep = "")

          # Return a list containing the filename
          list(
            src = outfile,
            contentType = "image/png",
            width = 700,
            height = 800,
            alt = "This is alternate text"
          )
        }
      },
      deleteFile = TRUE
    )
  })

  output$pathview.genes <- renderText({
    paste0(input$selectedGenes, collapse = ", ")
  })

  output$pathview.id <- renderText({
    paste0(input$selectedPath, collapse = ", ")
  })

  observeEvent(input$showPathview, {
    output$pathview <- renderImage(
      {
        if (length(input$selectedGenes) && length(input$selectedPath)) {
          selected.genes <- as.vector(strsplit(input$selectedGenes, ",")[[1]])

          find.labels <- expression.matrix.symbols[expression.matrix.symbols$NAME %in% selected.genes, ]
          labels <- find.labels$ensembl

          do_pathview(input$selectedPath, labels)

          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- paste("mmu", as.character(input$selectedPath), ".png", sep = "")

          # Return a list containing the filename
          list(
            src = outfile,
            contentType = "image/png",
            width = 700,
            height = 800,
            alt = "This is alternate text"
          )
        }
      },
      deleteFile = TRUE
    )
  })

  observeEvent(input$showPathview, {output$pathview2 <- 
    renderImage(
      {
        if (length(input$selectedGenes) && length(input$selectedPath)) {
          selected.genes <- as.vector(strsplit(input$selectedGenes, ",")[[1]])

          find.labels <- expression.matrix.symbols[expression.matrix.symbols$NAME %in% selected.genes, ]
          labels <- find.labels$ensembl

          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- paste("mmu", as.character(input$selectedPath), ".pathview.png", sep = "")

          # Return a list containing the filename
          list(
            src = outfile,
            contentType = "image/png",
            width = 700,
            height = 800,
            alt = "This is alternate text"
          )
        }
      },
      deleteFile = TRUE
    )
  })
}

start.time <- Sys.time()
analysis.erichment <- memo_do_enrichment(expression.matrix.freeze, c("PCA1"), length(rownames(expression.matrix.freeze)))
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

shinyApp(ui, server)
