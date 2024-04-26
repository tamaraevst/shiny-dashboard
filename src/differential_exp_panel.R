## differential_exp.R ##

library(DESeq2)

DE.ui <- function(id, metadata.groups) {
  tabItem(
    tabName = "DE",
    fluidRow(
      box(selectInput(inputId = NS(id, "metaname.de"), label = "Choose your group", choices = colnames(metadata.groups), selected = colnames(metadata.groups)[1]), width = 3),
      box(selectInput(inputId = NS(id, "condition1.de"), label = "Choose comparison condition #1", choices = unique(metadata.groups[, 1]), selected = unique(metadata.groups[, 1])[1]), width = 3),
      box(selectInput(inputId = NS(id, "condition2.de"), label = "Choose comparison condition #2", choices = unique(metadata.groups[, 1]), selected = unique(metadata.groups[, 1])[2]), width = 3),
      box(DT::dataTableOutput(NS(id, "summary.deseq2"), fill = TRUE),
        width = 12, title = "Summary table of DE analysis"
      )
    ),
    fluidRow(
      box(plotOutput(NS(id, "volcano.de"), height = 500, brush = NS(id, "brush.volcano")), width = 9, title = "Volcano plot"),
      box(sliderInput(NS(id, "p.threshold"),
        label = "Select p-value threshold",
        min = 0.001, value = 0.05, max = 0.05, step = 0.01, ticks = TRUE
      ), width = 3),
      box(sliderInput(NS(id, "fc.threshold"),
        label = "Select FC threshold",
        min = 0.5, value = 1, max = 10, step = 0.5, ticks = TRUE
      ), width = 3),
      box(selectInput(inputId = NS(id, "diffexp.capping"), label = "Do you want to cap log10pval values?", choices = c("Yes", "No"), selected = "Yes"), width = 3)
    ),
    fluidRow(box(DT::dataTableOutput(NS(id, "table.volcano"), fill = TRUE), width = 10)),
  )
}

DE.server <- function(id, expression.matrix, metadata, metadata.groups, rows.selected) {
  moduleServer(
    id,
    function(input, output, session) {
      expression.matrix.sorted <- reactive({
        id <- sub(".*?mmu(.*?)_.*", "\\1", rows.selected()[1])
        id.name <- rows.selected()[1]
        select_genes_for_pathway(expression.matrix, id, id.name, metadata, input$condition1.de, input$condition2.de, input$metaname.de)
      })

      metadata.sorted <- reactive({
        select_meta_for_pathway(metadata.groups, input$condition1.de, input$condition2.de, input$metaname.de)
      })

      deseq2.results <- reactive({
        DEanalysis_deseq2(
          expression.matrix = expression.matrix.sorted(),
          condition = metadata.sorted(),
          var1 = input$condition1.de,
          var2 = input$condition2.de
        )
      })
      output$summary.deseq2 <- DT::renderDataTable({
        DT::datatable(deseq2.results())
      })

      volcano.plot <- reactive({
        myplot <- volcano_plot(deseq2.results(),
          pval.threshold = input$p.threshold,
          lfc.threshold = input$fc.threshold,
          log10pval.cap = input$diffexp.capping
        )
        myplot
      })
      output$volcano.de <- renderPlot(volcano.plot())

      output$table.volcano <- DT::renderDataTable({
        display.table <- deseq2.results() %>%
          dplyr::mutate(gene = .data$gene_name, log10pvalAdj = -log10(.data$pvalAdj)) %>%
          dplyr::filter(!is.na(.data$log10pvalAdj)) %>%
          dplyr::select(gene_name, log2FC, pval, pvalAdj, log10pvalAdj)
        brushedPoints(display.table, input$brush.volcano)
      })
    }
  )
}

DEanalysis_deseq2 <- function(
    expression.matrix,
    condition,
    var1,
    var2) {
  expression.matrix <- round(expression.matrix)
  expression.matrix <- # Remove genes of constant expression
    expression.matrix[matrixStats::rowMins(expression.matrix) !=
      matrixStats::rowMaxs(expression.matrix), ]

  condition <- factor(condition, levels = unique(condition))
  design <- stats::model.matrix(~ 0 + condition)

  deseq <- DESeq2::DESeqDataSetFromMatrix(expression.matrix, data.frame(condition), design)
  DESeq2::sizeFactors(deseq) <- stats::setNames(
    rep(1, ncol(expression.matrix)),
    colnames(expression.matrix)
  )
  deseq <- DESeq2::DESeq(deseq)
  if (which(condition == var1)[1] < which(condition == var2)[1]) contrast <- c(-1, 1) else contrast <- c(1, -1)
  deseq.res <- DESeq2::results(deseq, contrast = contrast)

  find.labels <- expression.matrix.symbols[expression.matrix.symbols$ensembl %in% rownames(expression.matrix), ]
  labels <- find.labels$NAME

  gene_id <- NULL
  pval <- NULL
  deseq2.output <- tibble::tibble(
    gene_name = labels,
    log2exp = log2(rowMeans(expression.matrix)),
    log2FC = deseq.res$log2FoldChange,
    pval = deseq.res$pvalue,
    pvalAdj = stats::p.adjust(pval, method = "BH"),
  )

  return(deseq2.output)
}

volcano_plot <- function(
    genes.de.results,
    pval.threshold = 0.05,
    lfc.threshold = 1,
    alpha = 0.1,
    log10pval.cap = "Yes") {
  df <- genes.de.results %>%
    dplyr::mutate(gene = .data$gene_name, log10pvalAdj = -log10(.data$pvalAdj)) %>%
    dplyr::filter(!is.na(.data$log10pvalAdj))

  if (all(df$log10pvalAdj >= -10)) log10pval.cap <- "No"
  if (log10pval.cap == "Yes") df$log10pvalAdj[df$log10pvalAdj < -10] <- -10

  vp <- ggplot(data = df, aes(x = log2FC, y = log10pvalAdj, colour = log2exp)) +
    viridis::scale_color_viridis() +
    viridis::scale_fill_viridis() +
    geom_point(size = 3) +
    theme_minimal() +
    xlab("log2(FC)") +
    ylab("-log10(pvalAdj)") +
    geom_vline(xintercept = lfc.threshold, colour = "grey") +
    geom_vline(xintercept = -lfc.threshold, colour = "grey") +
    geom_hline(yintercept = -log10(pval.threshold), colour = "blue")

  max.abs.lfc <- max(abs(df[df$log10pvalAdj > -Inf, ]$log2FC))
  vp <- vp + xlim(-max.abs.lfc, max.abs.lfc)

  return(vp)
}

# DEtable <- DEanalysis_deseq2(
#           expression.matrix = expression.matrix.freeze[, c(1,2,3,4)],
#           condition = meta[c(1,2,3,4), 2],
#           var1 = meta[1,2],
#           var2 = meta[3,2]
#         )

# print(DEtable)
