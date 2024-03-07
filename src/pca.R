## pca.R ##

library(ggplot2)
library(dplyr)

source(here::here("src/utils.R"))

# Plotting function of the PCA based off the selected genes and samples chosen for grouping

plot_pca <- function(
    expression.matrix,
    metadata,
    annotation.id = colnames(metadata)[1],
    n.abundant = NULL,
    show.labels = FALSE,
    show.ellipses = TRUE,
    label.force = 1) {
  annotation.name <- annotation.id
  n.abundant <- min(n.abundant, nrow(expression.matrix))

  expr.PCA.list <- sort_genes(expression.matrix, n.abundant) %>%
    t()
  expr.PCA.list <- expr.PCA.list[, apply(expr.PCA.list, 2, function(x) max(x) != min(x))] %>%
    stats::prcomp(center = TRUE, scale = TRUE)

  expr.PCA <- dplyr::mutate(
    as.data.frame(expr.PCA.list$x),
    name = factor(metadata[, 1], levels = metadata[, 1]),
    condition = if (!is.factor(metadata[, annotation.id])) {
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
    } else {
      metadata[, annotation.id]
    }
  )
  if (min(table(metadata[, annotation.id])) <= 2) {
    expr.PCA.2 <- expr.PCA
    expr.PCA.2$PC1 <- expr.PCA.2$PC1 * 1.001
    expr.PCA.2$PC2 <- expr.PCA.2$PC2 * 1.001
    expr.PCA.full <- rbind(expr.PCA, expr.PCA.2)
  } else {
    expr.PCA.full <- expr.PCA
  }
  pca.plot <- ggplot(expr.PCA.full, aes(x = .data$PC1, y = .data$PC2, colour = .data$condition)) +
    theme_light() +
    geom_point() +
    labs(
      x = paste0("PC1 (proportion of variance = ", summary(expr.PCA.list)$importance[2, 1] * 100, "%)"),
      y = paste0("PC2 (proportion of variance = ", summary(expr.PCA.list)$importance[2, 2] * 100, "%)"),
      colour = annotation.name
    )
  if (show.ellipses) {
    pca.plot <- pca.plot +
      ggforce::geom_mark_ellipse(aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
  }
  if (show.labels) {
    pca.plot <- pca.plot +
      ggrepel::geom_label_repel(
        data = expr.PCA,
        mapping = aes(x = .data$PC1, y = .data$PC2, colour = .data$condition, label = .data$name),
        max.overlaps = nrow(expr.PCA),
        force = label.force,
        label.size = 0.1,
        point.size = NA
      )
  }

  return(pca.plot)
}

# Ui for the PCA window

PCA.ui <- tabItem(
  "PCAplot",
  fluidRow(
    box(plotOutput("pca", height = 400), width = 8),
    box(
      title = "Select number of genes",
      sliderInput("pca.n.abundant",
        label = "Number of genes",
        min = 50, value = 500, max = 5000, step = 50, ticks = TRUE
      ), width = 4
    ),
    box(selectInput(inputId = "States_List", label = "Choose your group", choices = colnames(meta)), width = 4),
    box(selectInput(inputId = "Enable_labels", label = "Sample labels on", choices = c(TRUE, FALSE)), width = 3, height = 2)
  )
)

# expression.matrix <- as.matrix(read.csv(
#    "expression_matrix_preprocessed.csv", sep=",", row.names=1
# ))

# meta <- data.frame(
#   srr = colnames(expression.matrix),
#   timepoint = rep(c("0h", "12h", "36h"), each = 2)
# )

# plot_pca(
#         expression.matrix = expression.matrix,
#         metadata = meta,
#         # annotation.id = match(input[['pca.annotation']], colnames(meta)),
#         n.abundant = 200,
#         show.labels = TRUE,
#         show.ellipses = TRUE,
#       )
