source(here::here("src/utils.R"))

expression.matrix.freeze <- as.matrix(read.csv(
  "example_data/expression_matrix_preprocessed.csv",
  sep = ",", row.names = 1
))

meta <- data.frame(
  srr = colnames(expression.matrix.freeze),
  timepoint = rep(c("0h", "12h", "36h"), each = 2)
)

expression.matrix.symbols <- convert_from_esembl_to_symbol(expression.matrix.freeze)

gmt.file <- file.path("example_data/20221221_kegg_mmu.gmt")
