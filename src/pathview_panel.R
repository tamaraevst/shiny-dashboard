library(pathview)

do_pathview <- function(
    path.id,
    path.id.name,
    metadata,
    condition1,
    condition2,
    annotation.id) {
  pathways <- gmtPathways(gmt.file)

  find.index <- names(pathways)[grep(path.id.name, names(pathways))]
  find.genes <- scan(text = pathways[find.index][[1]], what = "", sep = ",")

  sorted.genes <- intersect(find.genes, rownames(expression.matrix.freeze))

  expression.matrix <- expression.matrix.freeze[sorted.genes, ] %>% t()

  expr.new <- dplyr::mutate(
    as.data.frame(expression.matrix),
    condition = if (!is.factor(metadata[, annotation.id])) {
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
    } else {
      metadata[, annotation.id]
    }
  )

  expr.new.refined1 <- expr.new %>%
    filter(condition == condition1) %>%
    select(-condition) %>%
    t()

  expr.new.refined2 <- expr.new %>%
    filter(condition == condition2) %>%
    select(-condition) %>%
    t()

  mean.frame <- cbind(rowMeans(expr.new.refined1), rowMeans(expr.new.refined2))

  fold.change <- apply(mean.frame, 1, function(x) log2(x[1] / x[2]))

  # fold.change <- apply(sort.expression.matrix, 1, function(x) log2(x[1]/x[2]))

  # entrez.names <- id2eg(names(fold.change), org = "Mm", unique.map=FALSE, category='ENSEMBL')
  # entrez.names <- geneannot.map(names(fold.change), in.type="ENSEMBL", out.type="SYMBOL", org = "Mm", unique.map=F)
  # print(entrez.names)
  # names(fold.change) <- as.data.frame(entrez.names)$ENTREZID
  # print(fold.change)

  entrez.names <- mapIds(org.Mm.eg.db, keys = names(fold.change), column = "ENTREZID", keytype = "ENSEMBL", multiVals = list)
  # print(entrez.names)

  names(entrez.names) <- NULL
  names(fold.change) <- entrez.names
  bar.lim <- round(max(abs(max(fold.change)), abs(min(fold.change))))

  pv.out <- pathview(gene.data = fold.change, pathway.id = as.character(path.id), species = "mmu", multi.state = TRUE, out.suffix = "pathview", limit = list(gene = bar.lim, cpd = bar.lim))
  # print(pv.out$plot.data.gene$labels)

  return(pv.out)
}

plot_density <- function(
    data.frame) {
  p <- ggplot(data.frame, aes(x = mol.data)) +
    geom_density(fill = "lightblue") +
    labs(x = "log2(fold change)", y = "Density")

  return(p)
}

plot_fc_of_gene <- function(
    gene.id,
    metadata,
    control = meta[1, 1],
    annotation.id = colnames(meta)[1]) {
  ensembl.names <- mapIds(org.Mm.eg.db, keys = gene.id, column = "ENSEMBL", keytype = "ENTREZID", multiVals = first)

  expr.new <- dplyr::mutate(
    as.data.frame(t(expression.matrix.freeze)),
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

get_data_from_the_click <- function(
    pv.out,
    x,
    y) {
  plot.data.gene <- as.data.frame(pv.out$plot.data.gene)
  plot.data.gene <- plot.data.gene[abs(plot.data.gene$x - x) <= 23 & abs(plot.data.gene$y - y) <= 8, ]

  if (length(plot.data.gene$kegg.names) != 0) {
    kegg.names <- plot.data.gene$kegg.names
    # convert.genes <- mapIds(org.Mm.eg.db, keys = kegg.names,
    # column = "ENSEMBL", keytype = "ENTREZID")
    # find.gene <- as.data.frame(expression.matrix.freeze[c(convert.genes), ])
    # colnames(find.gene) <- c("Normalised expression")
    labels <- plot.data.gene$labels
    data <- c(kegg.names, labels, plot.data.gene$mol.data)
  } else {
    data <- NULL
  }

  return(data)
}

Pathview.ui <- tabItem(
  "Pathview",
  p("In order to use Pathview Analysis page, please make sure you have selected the pathway of interest in the enrichment analysis results table. Below, you are now able to choose a group from your data and two samples to compare. A log2 fold change will be calculated for these selected samples and annotated onto KEGG pathway diagrams. You may click on the nodes (e.g. coloured in genes) on the diagram to extract gene ENTREZID and symbol names. A bar-plot showing the log2 fold change of this gene across samples will also be produced below."),
  fluidRow(
    box(selectInput(inputId = "metaname.pathview", label = "Choose your group", choices = colnames(meta), selected = colnames(meta)[1]), width = 3),
    box(selectInput(inputId = "condition1.pathview", label = "Choose comparison condition #1", choices = unique(meta[, 1]), selected = meta[1, 1]), width = 3),
    box(selectInput(inputId = "condition2.pathview", label = "Choose comparison condition #2", choices = unique(meta[, 1]), selected = meta[1, 2]), width = 3)
  ),
  fluidRow(box(
    style = "overflow-x: scroll;overflow-y: scroll;", imageOutput("pathview", height = 1100, width = 1100, fill = TRUE, click = "plot.click"), title = "Pathview diagram",
    verbatimTextOutput("info"), width = 12
  )),
  fluidRow(box(plotOutput("pathview.density"), title = "Density plot", footer = "This is a density plot of the fold change of all genes present in the selected pathway.", width = 12)),
  fluidRow(box(selectInput(inputId = "control.pathview", label = "Choose your control sample", choices = unique(meta[, 1]), selected = meta[1, 1]), width = 3)),
  fluidRow(box(plotOutput("bar.plot.foldchange"), title = "Bar plot", footer = "This plot shows the fold change of the selected gene across all samples belonging to the chosen group. You should change the 'Control' sample to adjust the fold change calculation to your needs.", width = 10))
)
