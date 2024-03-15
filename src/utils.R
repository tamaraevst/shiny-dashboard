## utils.R ##

# Sort by the most abundant genes

library(gprofiler2)

sort_genes <- function(
    expression.matrix,
    n.abundant) {
  n.abundant <- min(n.abundant, nrow(expression.matrix))

  expr.list <- expression.matrix %>%
    as.data.frame() %>%
    dplyr::filter(seq_len(nrow(expression.matrix)) %in%
      utils::tail(order(rowSums(expression.matrix)), n.abundant))

  expr.list
}

get_genes_from_enrichment <- function(
    table,
    row) {
  row.selected <- table()[row, 5]
  return(row.selected)
}

do_pca <- function(
    expression.matrix) {
  expr.PCA.list <- stats::prcomp(expression.matrix, center = TRUE, scale = TRUE)

  return(expr.PCA.list)
}

choose_pca_genes <- function(
    expression.matrix,
    pcas = c("PCA2"),
    n.abundant = 50) {
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

  matrix <- do_pca(expression.matrix)
  pca.contrib <- matrix$x

  if (length(pca.array) == 1) {
    ordered.pca.contrib <- pca.contrib[order(pca.contrib[, pca.array], decreasing = TRUE), pca.array]
    sort.pca.contrib <- as.data.frame(ordered.pca.contrib[1:n.abundant])
  } else {
    pca.contrib.abs.cropped <- pca.contrib[, pca.array]
    sort.pca.contrib <- sort_genes(pca.contrib.abs.cropped, n.abundant)
  }

  return(sort.pca.contrib)
}
