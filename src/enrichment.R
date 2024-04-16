library(fgsea)
library(memoise)

do_enrichment <- function(
    expression.matrix,
    pca,
    n.abundant = length(rownames(expression.matrix))) {
  # pathways <- getGenesets(org = "mmu", db = "kegg")
  # print(names(pathways))

  pathways <- gmtPathways(gmt.file)

  genes.list <- sort_pca_genes_by_scores(
    expression.matrix,
    pca = pca,
    n.abundant = n.abundant
  )

  # convert.genes <- mapIds(org.Mm.eg.db,
  #   keys = rownames(genes.list),
  #   column = "ENTREZID", keytype = "ENSEMBL"
  # )

  # names(convert.genes) <- NULL

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
  # fgseaResMain <- fgseaRes[, leadingEdge := mapIdsList(
  #                                      x=org.Mm.eg.db,
  #                                      keys=leadingEdge,
  #                                      keytype="ENTREZID",
  #                                      column="SYMBOL")]
  return(fgseaRes)
}

memo_do_enrichment <- memoise(do_enrichment)
