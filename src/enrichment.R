## enrichment.R ##

library(ggplot2)
library(gprofiler2)
library(plotly)

source(here::here("src/utils.R"))

# Function for the Manhattan plot enrichment analysis results. The user is prompted to select top abundant genes from the slider, these genes are then used for plotting.

plot_enrichment <- function(
    expression.matrix,
    n.abundant = 50,
    organism = "mmusculus",
    threshold = 0.05) {
  expr.list <- sort_genes(expression.matrix, n.abundant)

  genes <- rownames(expr.list)

  gostres <- gost(genes, organism = organism)

  graph <- gostplot(gostres, capped = TRUE, interactive = FALSE)

  graph
}


# Table with the summary from the enrichment analysis, based off the top abundant genes selected by the user

table_enrichment <- function(
    expression.matrix,
    n.abundant = 50,
    organism = "mmusculus") {
  expr.list <- sort_genes(expression.matrix, n.abundant)

  genes <- rownames(expr.list)

  gostres <- gost(genes, organism = organism)

  enrichment.table <- gostres$result[c("source", "term_id", "term_name", "p_value", "significant")]

  return(enrichment.table)
}


# Ui part for the Enrichmnet window

Enrichment.ui <- tabItem(
  "Enrichment",
  fluidRow(
    box(
      title = "Select number of genes",
      sliderInput("enrichment.n.abundant",
        label = "Number of genes",
        min = 10, value = 50, max = 500, step = 10, ticks = TRUE
      ), width = 4
    ),
    box(plotOutput("enrichment", height = 400), width = 8),
    box(DT::dataTableOutput("entable"),
      height = 100,
      width = 10
    )
  )
)
