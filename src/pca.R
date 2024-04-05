library(dplyr)
library(memoise)

do_pca <- function(
    expression.matrix) {
  expr.PCA.list <- stats::prcomp(expression.matrix, center = TRUE, scale = TRUE)

  return(expr.PCA.list)
}

do_pca_reverse <- function(
    expression.matrix) {
  result <- expression.matrix %>%
    t()
  result <- result[, apply(result, 2, function(x) max(x) != min(x))] %>%
    stats::prcomp(center = TRUE, scale = TRUE)

  return(result)
}

prcomp_reconstruct_data <- function(
    expr) {
  result <- t(t(expr$x %*% t(expr$rotation)) * expr$scale + expr$center)

  return(result)
}

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

sum_contrib_pca <- function(
    expr.PCA.list,
    pcas) {
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

sort_pca_genes_by_scores <- function(
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

memo_do_pca <- memoise(do_pca)