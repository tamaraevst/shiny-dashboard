## pca_samples_panel.R ##

library(ggplot2)
library(dplyr)

source(here::here("src/pca.R"))

PCAsamples.ui <- function(id, expression.matrix, metadata) {
  tabItem(
    tabName = "PCAgenes",
    fluidRow(
      box(
        DT::dataTableOutput(NS(id, "summary.pca.samples")),
        width = 15, title = "Summary of Princinpal Component Analysis performed on the samples"
      )
    ),
    fluidRow(
      box(
        plotOutput(NS(id, "variance.pca.samples")),
        width = 15, title = "Proportion of variance explained by each of the Principal Components"
      )
    ),
    fluidRow(
      box(plotOutput(NS(id, "pca.samples"), height = 500), width = 8, title = "PCA plot", footer = "The plot represents the scatter of data across chosen principal components. Note that a negative value means the component is not significant, whilst a positive value means the component is significant. The ellipses are added for better visualisation and represent here the smallest possible ellipses enclosing the given set of data-points. These ellipses are determined using the Khachiyan algorithm and are expanded by 5mm (by default)."),
      box(
        numericInput(NS(id, "pca.n.abundant"),
          label = "Select number of genes",
          min = 10, max = length(rownames(expression.matrix)), value = 500, step = 50
        ),
        width = 4
      ),
      box(selectInput(inputId = NS(id, "States_List"), label = "Choose your group", choices = colnames(metadata)), width = 4),
      box(selectInput(inputId = NS(id, "list.pcas.samples"), label = "Choose your principal components (note: please choose 2 components)", choices = c("PC1", "PC2", "PC3"), multiple = TRUE, selected = c("PC1", "PC2")), width = 4),
      box(selectInput(inputId = NS(id, "Enable_labels"), label = "Sample labels on?", choices = c("Yes", "No")), width = 3)
    )
  )
}

PCAsamples.server <- function(id, expression.matrix, metadata) {
  moduleServer(
    id,
    function(input, output, session) {
      pca.samples.plot <- reactive({
        myplot <- plot_pca_samples(
          expression.matrix = expression.matrix,
          metadata = metadata,
          pcas = input$list.pcas.samples,
          annotation.id = input$States_List,
          n.abundant = input$pca.n.abundant,
          show.labels = input$Enable_labels
        )
        myplot
      })
      output[["pca.samples"]] <- renderPlot(pca.samples.plot())

      compute.pca <- reactive({
        do_pca_reverse_order(expression.matrix, input$list.pcas.samples, input$pca.n.abundant)
      })
      mytable <- reactive({
        plot_pca_summary(compute.pca())
      })
      output$summary.pca.samples <- DT::renderDataTable({
        DT::datatable(mytable())
      })

      pca.samples.variance.plot <- reactive({
        myplot <- plot_pca_variance(
          compute.pca()
        )
        myplot
      })
      output[["variance.pca.samples"]] <- renderPlot(pca.samples.variance.plot())
    }
  )
}

# Plotting function of the PCA based off the selected genes and samples chosen for grouping

plot_pca_samples <- function(
    expression.matrix,
    metadata,
    pcas = c("PC1", "PC2"),
    annotation.id = colnames(metadata)[1],
    n.abundant = NULL,
    show.labels = "Yes") {
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

  annotation.name <- annotation.id
  n.abundant <- min(n.abundant, nrow(expression.matrix))

  expr.PCA.list <- do_pca_reverse_order(expression.matrix, pcas, n.abundant)

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
  if (show.labels == "Yes") {
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
