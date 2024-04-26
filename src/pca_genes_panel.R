## pca_genes_panel.R ##

library(ggplot2)
library(dplyr)

source(here::here("src/pca.R"))

# Ui for the PCA window

PCA.ui <- function(id, expression.matrix) {
  tabItem(
    "PCAplot",
    fluidRow(
      box(
        DT::dataTableOutput(NS(id, "summary.pca")),
        width = 15, title = "Summary of Princinpal Component analysis performed on the genes"
      )
    ),
    fluidRow(
      box(
        plotOutput(NS(id, "variance.pca")),
        width = 15, title = "Proportion of variance explained by each of the Principal Components"
      )
    ),
    fluidRow(
      box(plotOutput(NS(id, "pca"), height = 500), width = 8, title = "PCA Plot", footer = "The plot represents the scatter of data across chosen principal components. The data-points (genes) are colour-coded by their contribution to the PCA. We note that the individual contributions to the PC are formally calculated as the
    ratio of the squared factor scores of a given data-point by the eigenvalue associated with that component. The contrbition appearing on the plot is a scalar and so takes on the form of a weighted sum of invidual contributions of the chosen principal components."),
      column(4, selectInput(inputId = NS(id, "list_pcas"), label = "Choose your principal components (note: please choose 2 components)", choices = c("PC1", "PC2", "PC3"), multiple = TRUE, selected = c("PC1", "PC2"))),
      column(4, numericInput(inputId = NS(id, "pca.n.abundant.genes"), label = "Select number of genes sorted by their contribution to PCA", value = length(rownames(expression.matrix)), min = 10, max = length(rownames(expression.matrix)), step = 50)),
      column(4, sliderInput(NS(id, "pca.labels.size"),
        label = "Select labels size",
        min = 0.1, value = 1, max = 1.5, step = 1, ticks = TRUE
      )),
      column(4, selectInput(inputId = NS(id, "enable_labels"), label = "Would you like labels on?", choices = c(TRUE, FALSE), selected = FALSE))
    )
  )
}

PCA.server <- function(id, expression.matrix) {
  moduleServer(
    id,
    function(input, output, session) {
      matrix.pca <- memo_do_pca(expression.matrix)

      output$summary.pca <- DT::renderDataTable({
        DT::datatable(plot_pca_summary(matrix.pca))
      })

      pca.plot <- reactive({
        myplot <- plot_pca(
          matrix.pca,
          expression.matrix = expression.matrix,
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
          matrix.pca
        )
        myplot
      })
      output[["variance.pca"]] <- renderPlot(pca.variance.plot())
    }
  )
}
# Plotting function of the PCA based off the selected genes and samples chosen for grouping

plot_pca <- function(
    expr.PCA.list,
    expression.matrix,
    pcas = c("PC1", "PC2"),
    n.abundant = length(rownames(expression.matrix)),
    show.labels = FALSE,
    label.size = 0.1) {
  pca.array <- c()
  for (i in (1:length(pcas)))
  {
    if (pcas[i] == "PC1") {
      pca.array <- c(pca.array, 1)
    }

    if (pcas[i] == "PC2") {
      pca.array <- c(pca.array, 2)
    }

    if (pcas[i] == "PC3") {
      pca.array <- c(pca.array, 3)
    }
  }

  n.abundant <- min(n.abundant, nrow(expression.matrix))

  expr <- as.data.frame(expr.PCA.list$x)

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


# cat("---Analysed matrix---\n")

# matrix <- do_pca(expression.matrix.freeze)

# print(head(expression.matrix, 50))

# cat("---Original matrix---\n")

# og <- prcomp_reconstruct_data(matrix)

# print(head(og, 50))

# plot_pca(matrix)
