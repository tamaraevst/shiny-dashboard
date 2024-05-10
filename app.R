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
source(here::here("src/differential_exp_panel.R"))

expression.matrix.freeze <- as.matrix(read.csv(
  "example_data/expression_matrix_preprocessed.csv",
  sep = ",", row.names = 1
))

meta <- data.frame(
  srr = colnames(expression.matrix.freeze),
  timepoint = rep(c("0h", "12h", "36h"), each = 2)
)

meta.groups <- subset(meta, select = -1)

expression.matrix.symbols <- convert_from_esembl_to_symbol(expression.matrix.freeze)

gmt.file <- file.path("example_data/20221221_kegg_mmu.gmt")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Your Data", tabName = "yourdata", icon = icon("dashboard")),
    menuItem("PCA (genes)", tabName = "PCAgenes"),
    menuItem("PCA (samples)", tabName = "PCAsamples"),
    menuItem("Enrichment Analysis", tabName = "Enrichment"),
    menuItem("Differential Expression", tabName = "DE"),
    menuItem("Pathway Visualisation", tabName = "Pathview")
  )
)

body <- dashboardBody(tabItems(
  # Data page
  Data.ui("contents"),
  # PCA on samples page
  PCAsamples.ui("summary.pca.samples", expression.matrix.freeze, meta),
  # PCA on genes page
  PCA.ui("summary.pca", expression.matrix.freeze),
  # Enrichment analysis and heatmap page
  Enrichment.ui("enrichment"),
  # Differential expression analysis page
  DE.ui("summary.deseq2", meta.groups),
  # Pathview analysis page
  Pathview.ui("pathview", meta.groups)
))

ui <- dashboardPage(
  dashboardHeader(title = "Dashboard for omics data"),
  sidebar,
  body,
  skin = "purple"
)

server <- function(input, output, session) {
  Data.server("contents", expression.matrix.symbols)
  PCAsamples.server("summary.pca.samples", expression.matrix.freeze, meta)
  PCA.server("summary.pca", expression.matrix.freeze)
  rows.selected <- Enrichment.server("enrichment", expression.matrix.freeze, gmt.file)
  deseq2.results <- DE.server("summary.deseq2", expression.matrix.freeze, meta, meta.groups, rows.selected)
  Pathview.server("pathview", expression.matrix.freeze, meta, meta.groups, rows.selected, deseq2.results[[1]], deseq2.results[[2]])
}

shinyApp(ui, server)
