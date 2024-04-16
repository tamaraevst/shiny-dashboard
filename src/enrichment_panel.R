## enrichment.R ##

library(ggplot2)
library(dplyr)
library(reshape)

source(here::here("src/data.R"))
source(here::here("src/utils.R"))

# Function for the Manhattan plot enrichment analysis results. The user is prompted to select top abundant genes from the slider, these genes are then used for plotting.

plot_enrichment <- function(
    expression.matrix,
    pca = c("PCA1"),
    n.abundant = length(rownames(expression.matrix.freeze)),
    dataset.organism = "mmusculus") {
  fgseaRes <- memo_do_enrichment(expression.matrix, pca, n.abundant)

  enrichment.plot <- ggplot(data = fgseaRes, aes(x = as.numeric(row.names(fgseaRes)), y = -log10(padj))) +
    geom_point(color = "purple") +
    labs(x = "Pathway index", y = "-log10(padj)") +
    geom_hline(yintercept = -log10(0.05)) +
    scale_shape_discrete(
      name = "Legend",
      labels = c("padj=0.05")
    )
  # scale_color_manual(
  # breaks = c('pval = 0.05'),
  # values = c('pval = 0.05' = 'purple')) +
  theme_minimal()

  enrichment.plot
}

bar_enrichment_up <- function(
    expression.matrix,
    pcas = c("PCA1"),
    n.abundant = length(rownames(expression.matrix.freeze)),
    dataset.organism = "mmusculus") {
  fgseaRes <- memo_do_enrichment(expression.matrix, pcas, n.abundant)

  topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n = 10), ]

  genes.up <- topPathwaysUp[, c("pathway", "padj", "size", "ES")]
  genes.up <- genes.up[order(genes.up$ES), ]

  bar.plot.up <- ggplot(genes.up, aes(x = ES, y = pathway, size = size, colour = log10(padj))) +
    geom_point(aes(fill = log10(padj)), colour = "black", shape = 21) +
    viridis::scale_color_viridis(option = "plasma", breaks = scales::pretty_breaks(n = 5)) +
    viridis::scale_fill_viridis(option = "plasma") +
    theme_minimal() +
    # scale_colour_continuous(low="purple", high="blue", name = "log10(padj)", limits = c(min(log10(genes.up$padj)), max(log10(genes.up$padj))),
    #         guide=guide_colorbar(reverse=TRUE)) +
    labs(x = "Enrichment score", y = "Pathway") +
    scale_y_discrete(limits = genes.up$pathway)

  # rng <- range(genes.up$padj)

  # bar.plot.up <- ggplot(genes.up, aes(x = size, y = pathway, fill = padj)) +
  #   scale_fill_continuous(
  #     low = "purple", high = "blue", name = "padj",
  #     guide = guide_colorbar(reverse = TRUE)
  #   ) +
  #   scale_fill_gradient2(low = "blue", mid = "cyan", high = "purple",
  #                      midpoint = mean(rng),  # Same midpoint for both plots
  #                      breaks = seq(-100, 100, 4),  # Breaks in the scale bar
  #                      limits = c(floor(rng[1]), ceiling(rng[2]))) +
  #   geom_col() +
  #   labs(x = "Size", y = "Pathway") +
  #   scale_y_discrete(limits = genes.up$pathway)

  # bar.plot.up + coord_flip()

  return(bar.plot.up)
}

bar_enrichment_down <- function(
    expression.matrix,
    pcas = c("PCA1"),
    n.abundant = length(rownames(expression.matrix.freeze)),
    dataset.organism = "mmusculus") {
  fgseaRes <- memo_do_enrichment(expression.matrix, pcas, n.abundant)

  topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n = 10), ]

  genes.down <- topPathwaysDown[, c("pathway", "padj", "size", "ES")]
  genes.down <- genes.down[order(genes.down$ES), ]

  bar.plot.down <- ggplot(genes.down, aes(x = ES, y = pathway, size = size, colour = log10(padj))) +
    geom_point(aes(fill = log10(padj)), colour = "black", shape = 21) +
    viridis::scale_color_viridis(option = "plasma", breaks = scales::pretty_breaks(n = 5)) +
    viridis::scale_fill_viridis(option = "plasma") +
    theme_minimal() +
    labs(x = "Enrichment score", y = "Pathway") +
    scale_y_discrete(limits = genes.down$pathway)

  return(bar.plot.down)
}

# Table with the summary from the enrichment analysis, based off the top abundant genes selected by the user

table_enrichment <- function(
    expression.matrix,
    pcas = c("PCA1"),
    n.abundant = length(rownames(expression.matrix.freeze)),
    dataset.organism = "mmusculus") {
  fgseaRes <- memo_do_enrichment(expression.matrix, pcas, n.abundant)

  pathways <- fgseaRes$pathway
  pval <- fgseaRes$pval
  padj <- fgseaRes$padj
  size <- fgseaRes$size
  leadingEdge <- fgseaRes$leadingEdge

  enrichment.table <- data.frame(Pathways = pathways, pval = pval, padj = padj, Size = size, Genes = sapply(leadingEdge, paste, collapse = ","))

  # genes.list <- choose_pca_genes(expression.matrix, pcas = pcas, n.abundant = n.abundant)

  # genes <- rownames(genes.list)

  # gostres <- gost(genes, ordered_query = TRUE, organism = dataset.organism, user_threshold = threshold, evcodes = TRUE)

  # enrichment.table.sorted <- gostres$result[c("term_id", "term_name", "p_value", "intersection_size", "intersection")]

  # is.num <- sapply(enrichment.table.sorted, is.numeric)
  # enrichment.table.sorted[is.num] <- lapply(enrichment.table.sorted[is.num], format, digits = 3)

  # enrichment.table <- enrichment.table.sorted

  return(enrichment.table)
}

heatmap_enrichment <- function(
    expression.matrix,
    genes.list) {
  genes.list <- as.vector(strsplit(genes.list, ",")[[1]])

  expression.matrix.sorted <- as.data.frame(expression.matrix) %>% filter(row.names(expression.matrix) %in% genes.list)

  find.labels <- expression.matrix.symbols[expression.matrix.symbols$ensembl %in% rownames(expression.matrix.sorted), ]
  labels <- find.labels$NAME

  # convert.genes <- AnnotationDbi::select(getExportedValue('org.Mm.eg.db', 'org.Mm.eg.db'),
  #   keys = genes.list,
  #   columns = "SYMBOL", keytype = "ENSEMBL", multiVals="first"
  # )
  # print(convert.genes)

  scaled.expression.matrix <- as.data.frame(scale(expression.matrix.sorted))

  data.genes <- cbind(genes = labels, data.frame(scaled.expression.matrix, row.names = NULL))
  # print(data.genes)

  expression.marix.melt <- melt(data.genes)

  enrichment.heatmap <- ggplot(expression.marix.melt, aes(genes, variable)) + # Create heatmap with ggplot2
    geom_tile(aes(fill = value)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    viridis::scale_color_viridis(discrete = FALSE) +
    viridis::scale_fill_viridis(discrete = FALSE) +
    labs(x = "", y = "", fill = "z-score of expression")

  return(enrichment.heatmap)
}

# Ui part for the Enrichmnet window

Enrichment.ui <- tabItem(
  "Enrichment",
  fluidRow(
    box(plotOutput("enrichment", height = 400, hover = hoverOpts(id = "plot.hover")),
      verbatimTextOutput("hover.info"),
      width = 8, title = "Scatter of enrichmnet analysis results", footer = "Above we display the overall distribution of p-values across all enrichment pathways. Here enrichmnet pathways are labelled by their index on the x-axis and their corresponding p-values on the y-axis. When this page is loaded, it performs Gene Set Enrichment Analysis using a ranked list of genes from your privided data set. The genes are ranked by their PCA scores, or equivalently by the coordinates of the data projected onto the PCs. The first instance of the page diplays the results of the enrichmnent analysis using default pricincipal components to rank the genes (you can see what they are on the right-hand side panel). You should change the principal components if you are interested in different analysis."
    ),
    box(selectInput(inputId = "list.pcas.enrichment", label = "Choose your principal component from the list", choices = c("PCA1", "PCA2", "PCA3"), multiple = FALSE, selected = c("PCA1")), width = 4),
    # box(
    #   sliderInput("enrichment.n.abundant",
    #     label = "Select top number of genes ordered by PCA results",
    #     min = 10, value = length(rownames(expression.matrix)), max = length(rownames(expression.matrix)), step = 500, ticks = TRUE
    #   ),
    #   width = 4
    # ),
    box(plotOutput("upgenes", height = 500), width = 6, title = "Bar plot of top 10 enrichmnent pathways with positive enrichment score", footer = "A positive enrichment score means that a gene set belonging to a particular pathway is over-represented with respect to the list of ranked genes we performed enrichment analysis on."),
    box(plotOutput("downgenes", height = 500), width = 6, title = "Bar plot of top 10 enrichmnent pathways with negative enrichment score", footer = "A negative enrichment score means that a gene set belonging to a particular pathway is under-represented with respect to the list of ranked genes we performed enrichment analysis on."),
    box(DT::dataTableOutput("entable", fill = TRUE),
      width = 12, title = "Summary table of enrichment analysis", footer = "You can click on the rows of the table to extract more information on the pathways. In particular, a list and a heat-map of leading edge genes belogning to a chosen pathway will be displayed below."
    ),
    box(verbatimTextOutput("engenes"),
      width = 12, title = "Leading edge genes"
    ),
    box(plotOutput("heatmap"),
      width = 12, title = "Heatmap of leading edge genes"
    )
  )
  # fluidRow(box(actionButton("showPathview.enrichment", "Show Pathview", width = 200), footer = "Click on the button to generate the pathview results for a given set of genes and an enrichment pathway.")),
  # fluidRow(box(imageOutput("pathview.enrichment", height = 1000, fill = TRUE), width = 12)),
  # fluidRow(box(imageOutput("pathview.enrichment2", height = 1000, fill = TRUE), width = 12))
)

# expression.matrix <- as.matrix(read.csv(
#   "expression_matrix_preprocessed.csv",
#   sep = ",", row.names = 1
# ))

# pwys <- system.file("extdata/mmu_kegg_pwys.zip", package="EnrichmentBrowser")
# pathways <- getGenesets(org = "mmu", db = "kegg")

# print(head(pathways))

# genes.list = choose_pca_genes(
#     expression.matrix,
#     pcas = c("PCA1"),
#     n.abundant = 1000)

# convert.genes <- mapIds(org.Mm.eg.db, keys = rownames(genes.list),
#        column = "ENTREZID", keytype = "ENSEMBL")

# names(convert.genes) <- NULL

# genes.list <- cbind(convert.genes, data.frame(genes.list, row.names=NULL))

# gene.list <- as.data.frame(genes.list, colClasses = c("character", "numeric"), colnames=c('t', 'ID'))
# genes.list <- setNames(c(genes.list[,1], genes.list[,2]), c(genes.list$t, genes.list$ID))
# colnames(genes.list) <- c('t', 'ID')
# genes.list <- setNames(genes.list$ID, genes.list$t)

# print(head(genes.list))

# genes.list <- genes.list[!duplicated(names(genes.list))]

# str(genes.list)

# print(all(is.finite(genes.list)))

# pathways <- gmtPathways(gmt.file)
# str(head(pathways))

# fgseaRes <- fgsea(pathways, genes.list, minSize=15, maxSize=500, scoreType = "pos")

# print(dim(fgseaRes))

# p <- ggplot() +
#     geom_point(data = fgseaRes, aes(x = as.numeric(row.names(fgseaRes)), y = -log10(pval)), color='purple') +
#     labs(x='Pathway index', y = "-log10(pval)") +
#     geom_hline(yintercept=-log10(0.05)) +
#     theme_minimal()

# topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
# topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# print(topPathways)

# plot_enrichment(expression.matrix)

# bar_enrichment_down(expression.matrix)


# heatmap_enrichment(scaled.expression.matrix, "ENSMUSG00000065037,ENSMUSG00000037742")
