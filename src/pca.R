## pca.R ##

library(dplyr)
library(memoise)

# PCA on X

do_pca <- function(
    expression.matrix) {
  expr.PCA.list <- stats::prcomp(expression.matrix, center = TRUE, scale = TRUE)

  return(expr.PCA.list)
}

# PCA on X^T

do_pca_reverse <- function(
    expression.matrix) {
  result <- expression.matrix %>%
    t()
  result <- result[, apply(result, 2, function(x) max(x) != min(x))] %>%
    stats::prcomp(center = TRUE, scale = TRUE)

  return(result)
}

# PCA on X^T ordered by the expression

do_pca_reverse_order <- function(
    expression.matrix,
    pcas = c("PCA1", "PCA2"),
    n.abundant = NULL) {
  n.abundant <- min(n.abundant, nrow(expression.matrix))

  expr.PCA.list <- sort_genes(expression.matrix, n.abundant)
  expr.PCA.list <- do_pca_reverse(expr.PCA.list)

  return(expr.PCA.list)
}

# Get original matrix X from the PCA results

prcomp_reconstruct_data <- function(
    expr) {
  result <- t(t(expr$x %*% t(expr$rotation)) * expr$scale + expr$center)

  return(result)
}

# Useful diagnostics from the PCA: distances from the center, cos2 and contributions

get_pca_results <- function(ind.coord, data, eigenvalues) {
  eigenvalues <- eigenvalues[1:ncol(ind.coord)]

  pca.center <- rep(0, ncol(data))
  pca.scale <- rep(1, ncol(data))

  # Compute the square of the distance between an individual and the
  # center of gravity
  getdistance <- function(ind_row, center, scale) {
    return(sum(((ind_row - center) / scale)^2))
  }
  d2 <- apply(data, 1, getdistance, pca.center, pca.scale)

  # Compute the cos2
  cos2 <- function(ind.coord, d2) {
    return(ind.coord * ind.coord / d2)
  }
  ind.cos2 <- apply(ind.coord, 2, cos2, d2)

  # Individual contributions
  contrib <- function(ind.coord, eigenvalues, n.ind) {
    100 * (1 / n.ind) * (ind.coord * ind.coord / eigenvalues)
  }
  ind.contrib <- t(apply(ind.coord, 1, contrib, eigenvalues, nrow(ind.coord)))

  colnames(ind.coord) <- colnames(ind.cos2) <-
    colnames(ind.contrib) <- paste0("Dim.", 1:ncol(ind.coord))

  rnames <- rownames(ind.coord)
  if (is.null(rnames)) rnames <- as.character(1:nrow(ind.coord))
  rownames(ind.coord) <- rownames(ind.cos2) <- rownames(ind.contrib) <- rnames

  # Individuals coord, cos2 and contrib
  ind <- list(coord = ind.coord, cos2 = ind.cos2, contrib = ind.contrib)

  return(ind)
}

# Calculate cumulative contributions

sum_contrib_pca <- function(
    expr.PCA.list,
    pcas) {
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

  # expr.PCA.list <- do_pca(expression.matrix)
  expr <- as.data.frame(expr.PCA.list$x)
  data <- prcomp_reconstruct_data(expr.PCA.list)

  ind <- get_pca_results(expr, data, expr.PCA.list$sdev * expr.PCA.list$sdev)

  contrib <- ind$contrib[, pca.array]
  if (length(pca.array) > 1) {
    eig <- expr.PCA.list$sdev[pca.array] * expr.PCA.list$sdev[pca.array]
    # Adjust variable contributions by the Dimension eigenvalues
    contrib <- t(apply(
      contrib, 1,
      function(var.contrib, pc.eig) {
        var.contrib * pc.eig
      },
      eig
    ))
    contrib <- apply(contrib, 1, sum) / sum(eig)
  }

  return(contrib)
}

# Sort the result from PCA on X^T by the values of the loading matrix

sort_pca_genes_by_loading <- function(
    expression.matrix,
    pca = c("PC1"),
    n.abundant = 500) {
  if (pca == "PC1") {
    pca.chosen <- 1
  }

  if (pca == "PC2") {
    pca.chosen <- 1
  }

  if (pca == "PC3") {
    pca.chosen <- 1
  }

  matrix <- do_pca_reverse(expression.matrix)

  pca.contrib <- matrix$rotation

  ordered.pca.contrib <- pca.contrib[order(pca.contrib[, pca.chosen], decreasing = TRUE), pca.chosen]
  sort.pca.contrib <- as.data.frame(ordered.pca.contrib[1:n.abundant])

  return(sort.pca.contrib)
}

memo_do_pca <- memoise(do_pca)
