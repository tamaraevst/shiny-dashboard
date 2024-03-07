## utils.R ##

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
