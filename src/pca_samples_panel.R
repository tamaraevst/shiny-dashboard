## pca.R ##

library(ggplot2)
library(dplyr)

source(here::here("src/pca.R"))

# Plotting function of the PCA based off the selected genes and samples chosen for grouping

plot_pca_samples <- function(
    expression.matrix,
    metadata,
    pcas = c("PCA1", "PCA2"),
    annotation.id = colnames(metadata)[1],
    n.abundant = NULL,
    show.labels = FALSE) {
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

  annotation.name <- annotation.id
  n.abundant <- min(n.abundant, nrow(expression.matrix))

  expr.PCA.list <- sort_genes(expression.matrix, n.abundant)
  expr.PCA.list <- do_pca_reverse(expr.PCA.list)

  expr.PCA.result <- as.data.frame(expr.PCA.list$x)
  expr.sorted <- expr.PCA.result[, pca.array]
  names(expr.sorted) <- c("x", "y")

  expr.PCA <- dplyr::mutate(
    as.data.frame(expr.sorted),
    name = factor(metadata[, 1], levels = metadata[, 1]),
    condition = if (!is.factor(metadata[, annotation.id])) {
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
    } else {
      metadata[, annotation.id]
    }
  )

  pca.plot <- ggplot(expr.PCA, aes(x = .data$x, y = .data$y, colour = .data$condition)) +
    theme_light() +
    viridis::scale_color_viridis(discrete = TRUE) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    geom_point() +
    ggforce::geom_mark_ellipse(aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE) +
    labs(
      x = paste0(pcas[1]),
      y = paste0(pcas[2]),
      colour = annotation.name
    )
  if (show.labels) {
    pca.plot <- pca.plot +
      ggrepel::geom_label_repel(
        data = expr.PCA,
        mapping = aes(x = .data$x, y = .data$y, colour = .data$condition, label = .data$name),
        max.overlaps = nrow(expr.PCA),
        label.size = 0.1,
        point.size = NA
      )
  }

  return(pca.plot)
}

# Ui for the PCA window

PCA.samples.ui <- tabItem(
  "PCAsamples",
  fluidRow(
    box(
      DT::dataTableOutput("summary.pca.samples"),
      width = 15, title = "Summary of Princinpal Component analysis performed on the samples"
    )
  ),
  fluidRow(
    box(
      plotOutput("variance.pca.samples"),
      width = 15, title = "Proportion of variance explained by each of the Principal Components"
    )
  ),
  fluidRow(
    box(plotOutput("pca.samples", height = 500), width = 8, title = "PCA plot", footer = "The plot represents the scatter of data across chosen principal components. Note that a negative value means the component is not significant, whilst a positive value means the component is significant. The ellipses are added for better visualisation and represent here the smallest possible ellipses enclosing the given set of data-points. These ellipses are determined using the Khachiyan algorithm and are expanded by 5mm (by default)."),
    box(
      numericInput("pca.n.abundant",
        label = "Select number of genes",
        min = 10, max = length(rownames(expression.matrix.freeze)), value = 500, step = 50
      ),
      width = 4
    ),
    box(selectInput(inputId = "States_List", label = "Choose your group", choices = colnames(meta)), width = 4),
    box(selectInput(inputId = "list.pcas.samples", label = "Choose your principal components (note: please choose 2 components)", choices = c("PCA1", "PCA2", "PCA3"), multiple = TRUE, selected = c("PCA1", "PCA2")), width = 4),
    box(selectInput(inputId = "Enable_labels", label = "Sample labels on?", choices = c(TRUE, FALSE)), width = 3, height = 2)
  )
)
