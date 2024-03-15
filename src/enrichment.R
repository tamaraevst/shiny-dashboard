## enrichment.R ##

library(ggplot2)
library(plotly)

source(here::here("src/utils.R"))

# Function for the Manhattan plot enrichment analysis results. The user is prompted to select top abundant genes from the slider, these genes are then used for plotting.

plot_enrichment <- function(
    expression.matrix,
    pcas = NULL,
    n.abundant = 50,
    dataset.organism = "mmusculus",
    threshold = 0.001) {
  genes.list <- choose_pca_genes(expression.matrix, pcas = pcas, n.abundant = n.abundant)

  genes <- rownames(genes.list)

  gostres <- gost(genes, ordered_query = TRUE, organism = dataset.organism, user_threshold = threshold)

  graph <- gostplot(gostres, capped = FALSE, interactive = FALSE)

  graph
}


# Table with the summary from the enrichment analysis, based off the top abundant genes selected by the user

table_enrichment <- function(
    expression.matrix,
    pcas = NULL,
    n.abundant = 50,
    dataset.organism = "mmusculus",
    threshold = 0.05) {
  genes.list <- choose_pca_genes(expression.matrix, pcas = pcas, n.abundant = n.abundant)

  genes <- rownames(genes.list)

  gostres <- gost(genes, ordered_query = TRUE, organism = dataset.organism, user_threshold = threshold, evcodes = TRUE)

  enrichment.table.sorted <- gostres$result[c("term_id", "term_name", "p_value", "intersection_size", "intersection")]

  is.num <- sapply(enrichment.table.sorted, is.numeric)
  enrichment.table.sorted[is.num] <- lapply(enrichment.table.sorted[is.num], format, digits = 3)

  enrichment.table <- enrichment.table.sorted

  return(enrichment.table)
}

# Ui part for the Enrichmnet window

Enrichment.ui <- tabItem(
  "Enrichment",
  fluidRow(
    box(plotOutput("enrichment", height = 400), width = 12),
    box(selectInput(inputId = "list.pcas.enrichment", label = "Choose your principal components from the list", choices = c("PCA1", "PCA2", "PCA3"), multiple = TRUE, selected = c("PCA1")), width = 4),
    box(
      sliderInput("enrichment.n.abundant",
        label = "Select top number of genes ordered by PCA results",
        min = 10, value = 50, max = 1000, step = 10, ticks = TRUE
      ),
      width = 4
    ),
    box(
      sliderInput("enrichment.threshold",
        label = "Select p-value threshold",
        min = 0.001, value = 0.05, max = 0.08, step = 0.005, ticks = TRUE
      ),
      width = 4
    ),
    box(DT::dataTableOutput("entable", fill = TRUE),
      width = 12
    ),
    box(verbatimTextOutput("engenes"),
      width = 12
    )
  )
)
