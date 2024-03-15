## pca.R ##

library(ggplot2)
library(dplyr)
library(factoextra)

source(here::here("src/utils.R"))

# Plotting function of the PCA based off the selected genes and samples chosen for grouping

plot_pca <- function(
    expr.PCA.list,
    pcas = c("PCA1", "PCA2"),
    select.genes = 10,
    label.force = FALSE,
    label.size = 3) {
  if (label.force == FALSE) {
    geom <- "point"
  } else {
    geom <- c("point", "text")
  }

  pca.array <- c()
  for (i in (1:length(pcas)))
  {
    if (pcas[i] == "PCA1") {
      pca.array <- c(pca.array, 1)
    }

    if (pcas[i] == "PCA2") {
      pca.array <- c(pca.array, 2)
    }

    if (pcas[i] == "PCA3") {
      pca.array <- c(pca.array, 3)
    }
  }

  res.ind <- get_pca_ind(expr.PCA.list)
  # print(res.ind$coord)          # Coordinates
  contrib.ind <- res.ind$contrib # Contributions to the PCs
  cos2.ind <- res.ind$cos2

  full.pca <- fviz_pca_ind(expr.PCA.list, col.ind = "contrib", axes = pca.array, geom = geom, label = "all", repel = TRUE, select.ind = list(contrib = select.genes), title = "", labelsize = label.size) +
    scale_color_gradient2(
      low = "yellow", mid = "pink",
      high = "blue"
    ) +
    theme_minimal()

  return(full.pca)
}

plot_pca_groups <- function(
    expr.PCA.list) {
  res.var <- get_pca_var(expr.PCA.list)

  groups.pca <- fviz_pca_var(expr.PCA.list,
    col.var = "contrib", # Color by contributions to the PC
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    repel = TRUE, # Avoid text overlapping
    title = ""
  ) + theme_minimal()

  return(groups.pca)
}

plot_pca_summary <- function(
    expr.PCA.list) {
  summary.pca <- summary(expr.PCA.list)
  table.summary <- summary.pca$importance

  return(table.summary)
}

plot_pca_variance <- function(expr.PCA.list) {
  variance.pca <- fviz_eig(expr.PCA.list, choice = "variance", geom = "line", linecolor = "blue", addlabels = TRUE, main = "") +
    theme_minimal()
  return(variance.pca)
}

# Ui for the PCA window

PCA.ui <- tabItem(
  "PCAplot",
  fluidRow(
    box(
      DT::dataTableOutput("summary.pca"),
      width = 15
    )
  ),
  fluidRow(
    box(
      plotOutput("variance.pca"),
      width = 15
    )
  ),
  fluidRow(
    title = "Principal Component Anaysis Plot", column(8, plotOutput("pca")),
    column(4, selectInput(inputId = "list_pcas", label = "Choose your principal components (note: please choose 2 components)", choices = c("PCA1", "PCA2", "PCA3"), multiple = TRUE, selected = c("PCA1", "PCA2"))),
    column(4, numericInput("pca.n.abundant", label = "Select number of genes sorted by their contribution to PCA", length(rownames(expression.matrix)), min = 10, max = length(rownames(expression.matrix)))),
    column(4, sliderInput("pca.labels.size",
      label = "Select labels size",
      min = 1, value = 3, max = 5, step = 1, ticks = TRUE
    )),
    column(4, selectInput(inputId = "enable_labels", label = "Would you like labels on?", choices = c(TRUE, FALSE), selected = FALSE))
  ),
  fluidRow(
    box(
      plotOutput("groups.pca"),
      width = 15
    )
  )
)
