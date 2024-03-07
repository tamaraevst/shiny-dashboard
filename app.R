## app.R ##
library(shiny)
library(shinydashboard)
library(here)
library(DT)

source(here::here("src/read.R"))
source(here::here("src/pca.R"))
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

  output$States_List <- renderUI({
    list_data <- colnames(meta)
    selectInput("group", "Select Meta Data Group", choices = list_data)
  })

  pca.plot <- reactive({
    myplot <- plot_pca(
      expression.matrix = expression.matrix,
      metadata = meta,
      annotation.id = input$States_List,
      n.abundant = input$pca.n.abundant,
      show.labels = input$Enable_labels,
      show.ellipses = TRUE,
    )
    myplot
  })
  output[["pca"]] <- renderPlot(pca.plot())

  enrichment.plot <- reactive({
    myplot <- plot_enrichment(
      expression.matrix = expression.matrix,
      n.abundant = input$enrichment.n.abundant,
    )
    myplot
  })
  output[["enrichment"]] <- renderPlot(enrichment.plot())

  enrichment.table <- reactive({
    mytable <- table_enrichment(
      expression.matrix = expression.matrix,
      n.abundant = input$enrichment.n.abundant
    )
    mytable
  })
  output[["entable"]] <- DT::renderDataTable(enrichment.table())
}

shinyApp(ui, server)
