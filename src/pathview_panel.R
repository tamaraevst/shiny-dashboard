## pathview_panel.R ##

library(pathview)

Pathview.ui <- function(id, metadata.groups) {
  tabItem(
    tabName = "Pathview",
    p("In order to use Pathview Analysis page, please make sure you have selected the pathway of interest in the enrichment analysis results table. The figure is downloaded directly from KEGG and then annotated using a log2 fold change of groups you have chosen on the Differential Expression Analysis page. You may click on the nodes (e.g. coloured in genes) on the diagram to extract gene ENTREZID and symbol names; you will also see the information on the statistical significance of that gene."),
    # fluidRow(box(selectInput(inputId = NS(id,"metaname.pathview"), label = "Choose your group", choices = colnames(metadata.groups), selected=colnames(metadata.groups)[1]), width = 3),
    # box(selectInput(inputId = NS(id,"condition1.pathview"), label = "Choose comparison condition #1", choices=unique(metadata.groups[, 1]), selected=unique(metadata.groups[,1])[1]), width = 3),
    # box(selectInput(inputId = NS(id,"condition2.pathview"), label = "Choose comparison condition #2", choices=unique(metadata.groups[, 1]), selected=unique(metadata.groups[,1])[2]), width = 3)),
    fluidRow(box(
      style = "overflow-x: scroll;overflow-y: scroll;", imageOutput(NS(id, "pathview"), height = 1100, width = 1100, fill = TRUE, click = NS(id, "plot.click")), title = "Pathview diagram",
      verbatimTextOutput(NS(id, "info")), width = 12
    )),
    fluidRow(box(plotOutput(NS(id, "pathview.density")), title = "Density plot", footer = "This is a density plot of the fold change of all genes present in the selected pathway.", width = 12)),
    fluidRow(box(selectInput(inputId = NS(id, "control.pathview"), label = "Choose your control sample", choices = unique(metadata.groups[, 1]), selected = metadata.groups[1, 1]), width = 3)),
    fluidRow(box(plotOutput(NS(id, "bar.plot.foldchange")), title = "Bar plot", footer = "This plot shows the fold change of the selected gene across all samples belonging to the chosen group. You should change the 'Control' sample to adjust the fold change calculation to your needs.", width = 10))
  )
}

Pathview.server <- function(id, expression.matrix, metadata, metadata.groups, rows.selected, de.results, metaname.pathview) {
  moduleServer(
    id,
    function(input, output, session) {
      # observe(updateSelectInput(session, "condition1.pathview", choices = unique(metadata.groups[, input$metaname.pathview]), selected=unique(metadata.groups[,input$metaname.pathview])[1]))
      # observe(updateSelectInput(session, "condition2.pathview", choices = unique(metadata.groups[, input$metaname.pathview]), selected=unique(metadata.groups[,input$metaname.pathview])[2]))
      # observe(updateSelectInput(session, "control.pathview", choices = unique(metadata.groups[, input$metaname.pathview])))

      pv.out <- reactive({
        id <- sub(".*?mmu(.*?)_.*", "\\1", rows.selected()[1])
        # id.name <- rows.selected()[[1]]
        result <- do_pathview(id, de.results())
      })

      output$pathview <- renderImage(
        {
          id <- sub(".*?mmu(.*?)_.*", "\\1", rows.selected()[1])

          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- paste("mmu", as.character(id), ".pathview.png", sep = "")

          # Return a list containing the filename
          list(
            src = outfile,
            contentType = "image/png"
            # width = 800,
            # height = 1000
          )
        },
        deleteFile = FALSE
      )

      output$info <- renderPrint({
        req(input$plot.click)
        x <- round(input$plot.click$x, 3)
        y <- round(input$plot.click$y, 3)
        cat("[", x, ", ", y, "] \n", sep = "")

        data <- get_data_from_the_click(pv.out(), x, y)

        if (length(data[1]) != 0) {
          cat("You have clicked on the node with KEGG name", data[1], "and label", data[2], "\n")
          find.gene <- de.results()[de.results()$gene_name == data[2], ]
          if (length(find.gene) != 0) {
            cat("Significance of the selected gene is as follows: pval = ", find.gene$pval, ", pvalAdj = ", find.gene$pvalAdj)
          } else {
            cat("For some reason I could not find a significance of this gene")
          }
        } else {
          cat("This node is not supported here. \n")
        }
      })

      pathview.density <- reactive({
        myplot <- plot_density(
          pv.out()$plot.data.gene
        )

        myplot
      })
      output$pathview.density <- renderPlot(pathview.density())

      # observe(updateSelectInput(session, "control.pathview", choices = unique(metadata.groups[, de.results()$group_name[1]]), selected = unique(metadata.groups[, de.results()$group_name[1]])[1]))

      bar.plot.foldchange <- reactive({
        req(input$plot.click)
        x <- round(input$plot.click$x, 3)
        y <- round(input$plot.click$y, 3)
        data <- get_data_from_the_click(pv.out(), x, y)

        plot <- plot_fc_of_gene(
          expression.matrix = expression.matrix,
          gene.id = data[1],
          metadata = metadata,
          control = input$control.pathview,
          annotation.id = metaname.pathview()
        )

        plot
      })
      output$bar.plot.foldchange <- renderPlot(bar.plot.foldchange())
    }
  )
}

# Function to generate pathview figures from KEGG, see https://bioconductor.org/packages/release/bioc/html/pathview.html for further details!
do_pathview <- function(
    path.id,
    de.results) {
  fold.change <- de.results$log2FC

  # print(colnames(expr.new))
  names(fold.change) <- de.results$ensembl

  bar.lim <- round(max(abs(max(fold.change)), abs(min(fold.change))))
  pv.out <- pathview(gene.data = fold.change, pathway.id = as.character(path.id), gene.idtype = "ENSEMBL", species = "mmu", multi.state = TRUE, out.suffix = "pathview", limit = list(gene = bar.lim, cpd = bar.lim))

  return(pv.out)
}

# Pathview results with fold change calculations allowed for sample-vs-sample.

do_pathview_manual <- function(
    expression.matrix,
    pathways.gmt.file,
    path.id,
    path.id.name,
    metadata,
    condition1,
    condition2,
    annotation.id) {
  pathways <- gmtPathways(pathways.gmt.file)

  # Here 'grep' is sensitive to pathways that have '(' and/or ')'. Please remove brackets if there any in your pathways file.
  find.index <- names(pathways)[base::grep(path.id.name, names(pathways))]
  find.genes <- scan(text = pathways[find.index][[1]], what = "", sep = ",")

  sorted.genes <- intersect(find.genes, rownames(expression.matrix))

  expression.matrix.sorted <- expression.matrix[sorted.genes, ] %>% t()

  expr.new <- dplyr::mutate(
    as.data.frame(expression.matrix.sorted),
    condition = if (!is.factor(metadata[, annotation.id])) {
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
    } else {
      metadata[, annotation.id]
    }
  )

  expr.new.refined1 <- expr.new %>%
    filter(condition == condition1) %>%
    dplyr::select(-condition) %>%
    t()

  expr.new.refined2 <- expr.new %>%
    filter(condition == condition2) %>%
    dplyr::select(-condition) %>%
    t()

  mean.frame <- cbind(rowMeans(expr.new.refined1), rowMeans(expr.new.refined2))

  fold.change <- apply(mean.frame, 1, function(x) log2(x[2] / x[1]))

  # print(colnames(expr.new))
  names(fold.change) <- colnames(expression.matrix.sorted)

  bar.lim <- round(max(abs(max(fold.change)), abs(min(fold.change))))

  pv.out <- pathview(gene.data = fold.change, pathway.id = as.character(path.id), gene.idtype = "ENSEMBL", species = "mmu", multi.state = TRUE, out.suffix = "pathview", limit = list(gene = bar.lim, cpd = bar.lim))

  return(pv.out)
}

# Pathview results for fold change calculations form DEseq2 results page.

plot_density <- function(
    data.frame) {
  x.lims <- max(abs(max(data.frame$mol.data)), abs(min(data.frame$mol.data)))
  p <- ggplot(data.frame, aes(x = mol.data)) +
    geom_density(fill = "lightblue") +
    labs(x = "log2(fold change)", y = "Density") +
    xlim(-1.01 * x.lims, 1.01 * x.lims)

  return(p)
}

# A bar plot showing FC values across samples for a given gene; requires the control sample to be chosen, against which the fold change is found

plot_fc_of_gene <- function(
    expression.matrix,
    gene.id,
    metadata,
    control = metadata[1, 1],
    annotation.id = colnames(metadata)[1]) {
  ensembl.names <- mapIds(org.Mm.eg.db, keys = gene.id, column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first")

  expr.new <- dplyr::mutate(
    as.data.frame(t(expression.matrix)),
    condition = if (!is.factor(metadata[, annotation.id])) {
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
    } else {
      metadata[, annotation.id]
    }
  )

  expr.new <- expr.new[c(ensembl.names, "condition")]
  colnames(expr.new)[1] <- "gene"

  new.data <- expr.new %>%
    group_by(condition) %>%
    summarise(mean = mean(gene))

  mean.cntrl <- subset(new.data, condition == control)$mean

  expr.others <- subset(new.data, condition != control)

  fold.change <- sapply(expr.others$mean, function(x) log2(mean.cntrl / x), simplify = "array")

  df <- data.frame(samples = expr.others$condition, logfc = fold.change)

  bar.plot.of.gene <- ggplot(data = df, aes(x = samples, y = logfc)) +
    geom_bar(stat = "identity", color = "blue", fill = "grey") +
    geom_text(aes(label = sprintf("%0.4f", logfc)), size = 3, position = position_stack(vjust = 0.5)) +
    theme_minimal() +
    xlab("Samples") +
    ylab("log2(fold change)")

  return(bar.plot.of.gene)
}

# Function to read the data from the pathview figure upon a click by a user

get_data_from_the_click <- function(
    pv.out,
    x,
    y) {
  plot.data.gene <- as.data.frame(pv.out$plot.data.gene)
  plot.data.gene <- plot.data.gene[abs(plot.data.gene$x - x) <= 23 & abs(plot.data.gene$y - y) <= 8, ]

  if (length(plot.data.gene$kegg.names) != 0) {
    kegg.names <- plot.data.gene$kegg.names
    labels <- plot.data.gene$labels
    data <- c(kegg.names, labels, plot.data.gene$mol.data)
  } else {
    data <- NULL
  }

  return(data)
}
