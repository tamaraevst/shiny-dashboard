## pca.R ##

library(ggplot2)
library(dplyr)

source(here::here("src/pca.R"))
source(here::here("src/data.R"))

get_eigenvalue <- function(
    expr) {
  eig <- (expr$sdev)^2

  return(eig)
}

# Plotting function of the PCA based off the selected genes and samples chosen for grouping

plot_pca <- function(
    expr.PCA.list,
    pcas = c("PCA1", "PCA2"),
    n.abundant = length(rownames(expression.matrix.freeze)),
    show.labels = FALSE,
    label.size = 0.1) {

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

  n.abundant <- min(n.abundant, nrow(expression.matrix.freeze))

  expr <- as.data.frame(expr.PCA.list$x)
  # cos2 <- expr * expr
  # data <- prcomp_reconstruct_data(expr.PCA.list)

  # ind <- get_pca_results(expr, data, expr.PCA.list$sdev *  expr.PCA.list$sdev)

  # contrib <- ind$contrib[, pca.array]
  # if(length(pca.array) > 1) {
  #   eig <- expr.PCA.list$sdev[pca.array] * expr.PCA.list$sdev[pca.array]
  #   # Adjust variable contributions by the Dimension eigenvalues
  #   contrib <- t(apply(contrib, 1,
  #                        function(var.contrib, pc.eig){var.contrib*pc.eig},
  #                        eig))
  #     contrib <-apply(contrib, 1, sum)/sum(eig)
  # }

  contrib <- sum_contrib_pca(expr.PCA.list, pcas)

  expr.sorted <- expr[, pca.array]
  names(expr.sorted) <- c("x", "y")

  names(contrib) <- NULL
  expr.sorted$contrib <- contrib

  expr.sorted <- tail(expr.sorted[order(expr.sorted$contrib), ], n.abundant)

  full.pca <- ggplot(expr.sorted, aes(x = .data$x, y = .data$y, colour = .data$contrib)) +
    theme_minimal() +
    viridis::scale_color_viridis() +
    viridis::scale_fill_viridis() +
    geom_point() +
    xlab(paste0(pcas[1])) +
    ylab(paste0(pcas[2]))

  if (show.labels) {
    find.labels <- expression.matrix.symbols[expression.matrix.symbols$ensembl %in% rownames(expr.sorted), ]
    labels <- find.labels$NAME

    full.pca <- full.pca +
      ggrepel::geom_label_repel(
        data = expr.sorted,
        mapping = aes(x = .data$x, y = .data$y, colour = .data$contrib, label = labels),
        max.overlaps = nrow(expr.sorted),
        # force = label.force,
        label.size = label.size,
        point.size = 1.0
      )
  }

  # full.pca <- fviz_pca_ind(expr.PCA.list, col.ind = "contrib", axes = pca.array, geom = geom, label = "all", repel = TRUE, select.ind = list(contrib = select.genes), title = "", labelsize = label.size) +
  #   scale_color_gradient2(
  #     low = "yellow", mid = "pink",
  #     high = "blue"
  #   ) +
  #   theme_minimal()

  return(full.pca)
}

plot_pca_biplot <- function(
    expr.PCA.list) {
  groups.pca <- fviz_pca_var(expr.PCA.list,
    col.var = "contrib", # Color by contributions to the PC
    gradient.cols = c("yellow", "pink", "blue"),
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

plot_pca_variance <- function(
    expr.PCA.list) {
  # variance.pca <- fviz_eig(expr.PCA.list, choice = "variance", geom = "line", linecolor = "blue", addlabels = TRUE, main = "") +
  #   theme_minimal()

  variance <- expr.PCA.list$sdev^2 / sum(expr.PCA.list$sdev^2)

  var.df <- data.frame(
    PC = paste0("PC", 1:6),
    variane = variance
  )

  variance.pca <- var.df %>% ggplot(aes(x = PC, y = variance, group = 1)) +
    theme_minimal() +
    geom_point(size = 4) +
    geom_col(fill = "slategray3") +
    geom_line(colour = "blue") +
    xlab("Principal Component") +
    ylab("Variance Explained")

  return(variance.pca)
}

# Ui for the PCA window

PCA.ui <- tabItem(
  "PCAplot",
  fluidRow(
    box(
      DT::dataTableOutput("summary.pca"),
      width = 15, title = "Summary of Princinpal Component analysis performed on the genes"
    )
  ),
  fluidRow(
    box(
      plotOutput("variance.pca"),
      width = 15, title = "Proportion of variance explained by each of the Principal Components"
    )
  ),
  fluidRow(
    box(plotOutput("pca", height = 500), width = 8, title = "PCA Plot", footer = "The plot represents the scatter of data across chosen principal components. The data-points (genes) are colour-coded by their contribution to the PCA. We note that the individual contributions to the PC are formally calculated as the
    ratio of the squared factor scores of a given data-point by the eigenvalue associated with that component. The contrbition appearing on the plot is a scalar and so takes on the form of a weighted sum of invidual contributions of the chosen principal components."),
    column(4, selectInput(inputId = "list_pcas", label = "Choose your principal components (note: please choose 2 components)", choices = c("PCA1", "PCA2", "PCA3"), multiple = TRUE, selected = c("PCA1", "PCA2"))),
    column(4, numericInput(inputId = "pca.n.abundant.genes", label = "Select number of genes sorted by their contribution to PCA", value = length(rownames(expression.matrix.freeze)), min = 10, max = length(rownames(expression.matrix.freeze)), step = 50)),
    column(4, sliderInput("pca.labels.size",
      label = "Select labels size",
      min = 0.1, value = 1, max = 1.5, step = 1, ticks = TRUE
    )),
    column(4, selectInput(inputId = "enable_labels", label = "Would you like labels on?", choices = c(TRUE, FALSE), selected = FALSE))
  )
)

# cat("---Analysed matrix---\n")

# matrix <- do_pca(expression.matrix.freeze)

# print(head(expression.matrix, 50))

# cat("---Original matrix---\n")

# og <- prcomp_reconstruct_data(matrix)

# print(head(og, 50))

# plot_pca(matrix)
