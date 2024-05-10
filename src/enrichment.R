## enrichment.R ##

library(fgsea)
library(memoise)

# Function to call enrichment analysis using Fast Gene Set Enrichment (fGSEA) implementation, see https://bioconductor.org/packages/release/bioc/html/fgsea.html for further details!

do_enrichment <- function(
    expression.matrix,
    pca,
    n.abundant = length(rownames(expression.matrix)),
    pathways.gmt.file) {
  pathways <- gmtPathways(pathways.gmt.file)

  genes.list <- sort_pca_genes_by_loading(
    expression.matrix,
    pca = pca,
    n.abundant = n.abundant
  )

  genes.list <- cbind(rownames(genes.list), data.frame(genes.list, row.names = NULL))

  colnames(genes.list) <- c("t", "ID")
  genes.list <- setNames(genes.list$ID, genes.list$t)

  # print(str(genes.list))
  # genes.list <- genes.list[!duplicated(names(genes.list))]

  minSize <- 15
  maxSize <- length(pathways)

  # nPermSimple = 10000
  cat("Begining Enrichment Analysis with the following parameters: \n")
  cat("- minSize = ", minSize, "\n")
  cat("- maxSize = ", maxSize, "\n")
  cat("- nPermSimple = 10000 \n")
  fgseaRes <- fgsea(pathways, genes.list, minSize = minSize, maxSize = maxSize, nPermSimple = 10000)
  cat("------------------------------\n")
  cat("Completed Enrichment Analysis \n")

  return(fgseaRes)
}

memo_do_enrichment <- memoise(do_enrichment)
