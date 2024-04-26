## utils.R ##

library(org.Mm.eg.db)
library(memoise)
library(dplyr)

# Sort by the most abundant genes

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

get_genes_from_column <- function(
    table,
    row,
    column) {
  row.selected <- table()[row, column]
  return(row.selected)
}

convert_from_esembl_to_symbol <- function(
    expression.matrix) {
  convert.genes <- AnnotationDbi::mapIds(org.Mm.eg.db,
    keys = rownames(expression.matrix),
    column = "SYMBOL", keytype = "ENSEMBL"
  )

  ensembl <- names(convert.genes)
  names(convert.genes) <- NULL

  expression.matrix.final <- cbind(
    genes = convert.genes,
    ensembl, data.frame(expression.matrix, row.names = NULL)
  ) %>% dplyr::mutate(NAME = ifelse(is.na(.data$genes), .data$ensembl, .data$genes))

  return(expression.matrix.final)
}

select_genes_for_pathway <- function(
    expression.matrix,
    path.id,
    path.id.name,
    metadata,
    condition1,
    condition2,
    annotation.id) {
  pathways <- gmtPathways(gmt.file)

  find.index <- names(pathways)[grep(path.id.name, names(pathways))]
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

  sorted.data <- cbind(expr.new.refined1, expr.new.refined2)

  return(sorted.data)
}

select_meta_for_pathway <- function(
    metadata,
    condition1,
    condition2,
    annotation.id) {
  metadata.sorted <- metadata[, annotation.id]

  metadata1 <- metadata.sorted[grepl(as.character(condition1), metadata.sorted)]
  metadata2 <- metadata.sorted[grepl(as.character(condition2), metadata.sorted)]

  metadata.final <- c(metadata1, metadata2)

  return(metadata.final)
}
