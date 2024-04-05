## utils.R ##

library(org.Mm.eg.db)
library(memoise)

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
  expression.matrix
)
{
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


