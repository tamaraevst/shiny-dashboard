## app.R ##
library(shiny)
library(shinydashboard)
library(here)
library(DT)

source(here::here("src/read.R"))
source(here::here("src/pca_ongenes.R"))
source(here::here("src/enrichment.R"))
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
    DT::datatable(expression.matrix)
  })

  matrix <- do_pca(expression.matrix)

  output$summary.pca <- DT::renderDataTable({
    DT::datatable(plot_pca_summary(matrix))
  })

  pca.plot <- reactive({
    myplot <- plot_pca(
      matrix,
      pcas = input$list_pcas,
      select.genes = input$pca.n.abundant,
      label.force = input$enable_labels,
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

  pca.groups.plot <- reactive({
    myplot <- plot_pca_groups(
      matrix
    )
    myplot
  })
  output[["groups.pca"]] <- renderPlot(pca.groups.plot())

  enrichment.plot <- reactive({
    myplot <- plot_enrichment(
      expression.matrix = expression.matrix,
      pcas = input$list.pcas.enrichment,
      n.abundant = input$enrichment.n.abundant,
      threshold = input$enrichment.threshold
    )
    myplot
  })
  output[["enrichment"]] <- renderPlot(enrichment.plot())

  enrichment.table <- reactive({
    mytable <- table_enrichment(
      expression.matrix = expression.matrix,
      pcas = input$list.pcas.enrichment,
      n.abundant = input$enrichment.n.abundant,
      threshold = input$enrichment.threshold
    )
    mytable
  })
  output[["entable"]] <- DT::renderDataTable(enrichment.table(),
    server = TRUE, rownames = FALSE, extensions = "Buttons", selection = "single",
    options = list(autoWidth = TRUE, buttons = c("csv", "excel"), dom = "Bfrtip", columnDefs = list(list(visible = FALSE, targets = 4)))
  )

  output$engenes <- renderPrint({
    s <- input$entable_rows_selected
    if (length(s)) {
      cat("The following genes fall under the selected enrichment category:\n\n")
      cat(enrichment.table()[s, 5])
    } else {
      cat("Click on rows to print the genes belonging to the selected category")
    }
  })
}

shinyApp(ui, server)
